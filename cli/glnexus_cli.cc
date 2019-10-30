// Basic GLnexus command-line interface for use on one compute node

#include <iostream>
#include <exception>
#include <fstream>
#include <getopt.h>
#include <sstream>
#include <cstdlib>
#include "vcf.h"
#include "hfile.h"
#include "service.h"
#include "unifier.h"
#include "BCFKeyValueData.h"
#include "RocksKeyValue.h"
#include "ctpl_stl.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_sinks.h"
#include "cli_utils.h"

using namespace std;

auto console = spdlog::stderr_logger_mt("GLnexus");
GLnexus::Status s;
#define H(desc,expr) \
    s = expr; \
    if (s.bad()) { \
        console->error("Failed to {}: {}", desc, s.str()); \
        if (getenv("DX_JOB_ID")) { \
            ofstream joberrorjson("/home/dnanexus/job_error.json"); \
            joberrorjson << "{\"error\": {\"type\": \"AppError\", \"message\": \""; \
            joberrorjson << "Failed to " << desc << ": " << s.str(); \
            joberrorjson << "\"}}"; \
            joberrorjson.close(); \
        } \
        return 1; \
    }

// Perform all the separate GLnexus operations in one go.
// return 0 on success, 1 on failure.
static int all_steps(const vector<string> &vcf_files,
                     const string &bedfilename,
                     const string &config_name,
                     bool more_PL, bool squeeze, bool trim_uncalled_alleles,
                     size_t mem_budget, size_t nr_threads,
                     bool debug,
                     bool iter_compare,
                     size_t bucket_size) {
    GLnexus::Status s;
    GLnexus::unifier_config unifier_cfg;
    GLnexus::genotyper_config genotyper_cfg;
    string cfg_txt, cfg_crc32c;

    if (vcf_files.empty()) {
        console->error("No source GVCF files specified");
        return 1;
    }

    H("load unifier/genotyper configuration",
        GLnexus::cli::utils::load_config(console, config_name, unifier_cfg, genotyper_cfg, cfg_txt, cfg_crc32c,
                                         more_PL, squeeze, trim_uncalled_alleles));

    // initilize empty database
    string dbpath("GLnexus.DB");
    vector<pair<string,size_t> > contigs;
    H("initialize database", GLnexus::cli::utils::db_init(console, dbpath, vcf_files[0], contigs,
                                                          bucket_size));

    {
        // sanity check, see that we can get the contigs back
        vector<pair<string,size_t> > contigs_dbg;
        H("read the contigs back from DB",
          GLnexus::cli::utils::db_get_contigs(console, dbpath, contigs_dbg));
        if (contigs_dbg != contigs)
            return GLnexus::Status::Invalid("error, contigs read from DB do not match originals");
    }

    // Load the GVCFs into the database
    {
        // use an empty range filter
        vector<GLnexus::range> ranges;
        H("bulk load into DB",
          GLnexus::cli::utils::db_bulk_load(console, mem_budget, nr_threads, vcf_files, dbpath, ranges, contigs, false));
    }

    if (iter_compare) {
        H("compare database iteration methods",
          GLnexus::cli::utils::compare_db_itertion_algorithms(console, dbpath, 50));
    }

    // discover alleles
    // TODO: overlap allele discovery with final compactions. have db_bulk_load output a RocksKeyValue pointer which we can reuse
    vector<GLnexus::range> ranges;
    if (bedfilename.empty()) {
        console->warn("Processing full length of {} contigs, as no --bed was provided. Providing a BED file with regions of interest, if applicable, can speed this up.", std::to_string(contigs.size()));
        for (int rid = 0; rid < contigs.size(); ++rid) {
            ranges.push_back(GLnexus::range(rid, 0, contigs[rid].second));
        }
    } else {
        H("parse the bed file", GLnexus::cli::utils::parse_bed_file(console, bedfilename, contigs, ranges));
    }
    GLnexus::discovered_alleles dsals;
    unsigned sample_count = 0;
    H("discover alleles",
      GLnexus::cli::utils::discover_alleles(console, mem_budget, nr_threads, dbpath, ranges, contigs, dsals, sample_count));
    if (debug) {
        string filename("/tmp/dsals.yml");
        console->info("Writing discovered alleles as YAML to {}", filename);
        H("serialize discovered alleles to a file",
          GLnexus::cli::utils::yaml_write_discovered_alleles_to_file(dsals, contigs, sample_count, filename));
    }

    // partition dsals by contig to reduce peak memory usage in the unifier
    std::vector<GLnexus::discovered_alleles> dsals_by_contig(contigs.size());
    for (auto p = dsals.begin(); p != dsals.end(); dsals.erase(p++)) {
        UNPAIR(*p, al, dai);
        assert(al.pos.rid >= 0 && al.pos.rid < contigs.size());
        dsals_by_contig[al.pos.rid][al] = dai;
    }

    // unify sites
    // we could parallelize over dsals_by_contig although this might increase memory usage.
    vector<GLnexus::unified_site> sites;
    GLnexus::unifier_stats stats;
    for (auto& dsals_i : dsals_by_contig) {
        GLnexus::unifier_stats stats1;
        H("unify sites",
          GLnexus::cli::utils::unify_sites(console, unifier_cfg, contigs, dsals_i, sample_count, sites, stats1));
        stats += stats1;
    }
    console->info("unified to {} sites cleanly with {} ALT alleles. {} ALT alleles were {} and {} were filtered out on quality thresholds.",
                  sites.size(), stats.unified_alleles, stats.lost_alleles,
                  (unifier_cfg.monoallelic_sites_for_lost_alleles ? "additionally included in monoallelic sites" : "lost due to failure to unify"),
                  stats.filtered_alleles);
    if (debug) {
        string filename("/tmp/sites.yml");
        console->info("Writing unified sites as YAML to {}", filename);
        H("write unified sites to file",
          GLnexus::cli::utils::write_unified_sites_to_file(sites, contigs, filename));
    }

    // genotype
    genotyper_cfg.output_residuals = debug;
    vector<string> hdr_lines = {
        ("##GLnexusConfigName="+config_name),
        ("##GLnexusConfigCRC32C="+cfg_crc32c),
        ("##GLnexusConfig="+cfg_txt)
    };
    auto DX_JOB_ID = std::getenv("DX_JOB_ID");
    if (DX_JOB_ID) {
        // if running in DNAnexus, record job ID in header
        hdr_lines.push_back(string("##DX_JOB_ID=")+DX_JOB_ID);
    }
    string outfile("-");
    H("genotype",
      GLnexus::cli::utils::genotype(console, mem_budget, nr_threads, dbpath, genotyper_cfg, sites, hdr_lines, outfile));

    return 0;
}


void help(const char* prog) {
    cout << "Usage: " << prog << " [options] /vcf/file/1 .. /vcf/file/N" << endl
         << "Merge and joint-call input gVCF files, emitting multi-sample BCF on standard output." << endl << endl
         << "Options:" << endl
         << "  --bed FILE, -b FILE            three-column BED file of ranges to analyze (if omitted, use full length of all contigs)" << endl
         << "  --config X, -c X               configuration preset name or .yml filename (default: gatk)" << endl
         << "  --more-PL, -P                  include PL from reference bands and other cases omitted by default" << endl
         << "  --squeeze, -S                  reduce pVCF size by suppressing detail in cells derived from reference bands" << endl
         << "  --trim-uncalled-alleles, -a    remove alleles with no output GT calls in postprocessing" << endl
         << "  --list, -l                     given files contain lists of gVCF filenames, one per line" << endl
         << "  --mem-gbytes X, -m X           memory budget, in gbytes (default: most of system memory)" << endl
         << "  --threads X, -t X              thread budget (default: all hardware threads)" << endl
         << "  --help, -h                     print this help message" << endl
         << endl << "Configuration presets:" << endl;
    cout << GLnexus::cli::utils::describe_config_presets() << endl;
}

// Expected usage:
//    glnexus [vcf files]
//
int main(int argc, char *argv[]) {
    GLnexus::Status s;
    spdlog::set_level(spdlog::level::info);
    spdlog::set_pattern("[%t] %+");
    #ifdef NDEBUG
    #define BUILD_CONFIG "release"
    #else
    #define BUILD_CONFIG "debug"
    #endif
    console->info("glnexus_cli {} {} {}", BUILD_CONFIG, GIT_REVISION, __TIMESTAMP__);
    GLnexus::cli::utils::detect_jemalloc(console);

    if (argc < 2) {
        help(argv[0]);
        return 1;
    }

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"bed", required_argument, 0, 'b'},
        {"config", required_argument, 0, 'c'},
        {"more-PL", no_argument, 0, 'P'},
        {"squeeze", no_argument, 0, 'S'},
        {"trim-uncalled-alleles", no_argument, 0, 'a'},
        {"list", no_argument, 0, 'l'},
        {"mem-gbytes", required_argument, 0, 'm'},
        {"threads", required_argument, 0, 't'},
        {"bucket_size", required_argument, 0, 'x'},
        {"debug", no_argument, 0, 'd'},
        {"iter_compare", no_argument, 0, 'i'},
        {0, 0, 0, 0}
    };

    int c;
    string config_name = "gatk";
    bool more_PL = false;
    bool squeeze = false;
    bool trim_uncalled_alleles = false;
    bool list_of_files = false;
    bool debug = false;
    bool iter_compare = false;
    string bedfilename;
    size_t mem_budget = 0, nr_threads = 0;
    size_t bucket_size = GLnexus::BCFKeyValueData::default_bucket_size;

    while (-1 != (c = getopt_long(argc, argv, "hPSadil:b:x:m:t:c:",
                                  long_options, nullptr))) {
        switch (c) {
            case 'b':
                bedfilename = string(optarg);
                if (bedfilename.size() == 0) {
                    cerr <<  "invalid BED filename" << endl;
                    return 1;
                }
                break;

            case 'l':
                list_of_files = true;
                break;

            case 'c':
                config_name = string(optarg);
                break;

            case 'P':
                more_PL = true;
                break;

            case 'S':
                squeeze = true;
                break;

            case 'a':
                trim_uncalled_alleles = true;
                break;

            case 'd':
                debug = true;
                break;

            case 'h':
            case '?':
                help(argv[0]);
                exit(0);
                break;

            case 'i':
                iter_compare = true;
                break;

            case 'x':
                bucket_size = strtoul(optarg, nullptr, 10);
                if (bucket_size == 0 || bucket_size > 1000000000) {
                    cerr << "bucket size should be in (1,1e9]" << endl;
                    return 1;
                }
                break;

            case 'm':
                mem_budget = strtoull(optarg, nullptr, 10);
                if (mem_budget == 0 || mem_budget > 16*1024) {
                    cerr << "invalid --mem-gbytes" << endl;
                    return 1;
                }
                mem_budget <<= 30;
                break;

            case 't':
                nr_threads = strtoull(optarg, nullptr, 10);
                if (nr_threads == 0 || nr_threads > 1024) {
                    cerr << "invalid --threads" << endl;
                    return 1;
                }
                break;

            default:
                abort ();
        }
    }

    if (optind > argc-1) {
        help(argv[0]);
        return 1;
    }

    vector<string> vcf_files, vcf_files_precursor;
    for (int i=optind; i < argc; i++) {
        vcf_files_precursor.push_back(string(argv[i]));
    }

    if (list_of_files) {
        for (const string& fn : vcf_files_precursor) {
            string gvcf;
            ifstream infile(fn);
            while (getline(infile, gvcf)) {
                vcf_files.push_back(gvcf);
            }
            if (infile.bad() || !infile.eof()) {
                H("read input file list", GLnexus::Status::IOError("reading", fn));
            }
        }
    } else {
        vcf_files = vcf_files_precursor;
    }

    return all_steps(vcf_files, bedfilename, config_name, more_PL, squeeze, trim_uncalled_alleles, mem_budget, nr_threads, debug, iter_compare, bucket_size);
}

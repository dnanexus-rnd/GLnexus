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
                     int nr_threads,
                     bool debug,
                     bool iter_compare,
                     size_t bucket_size) {
    GLnexus::Status s;
    GLnexus::unifier_config unifier_cfg;
    GLnexus::genotyper_config genotyper_cfg;
    string cfg_crc32c;

    if (vcf_files.empty()) {
        console->error("No source GVCF files specified");
        return 1;
    }

    H("load unifier/genotyper configuration",
        GLnexus::cli::utils::load_config(console, config_name, unifier_cfg, genotyper_cfg, cfg_crc32c));

    // initilize empty database
    string dbpath("GLnexus.DB");
    vector<pair<string,size_t> > contigs;
    H("initializing database", GLnexus::cli::utils::db_init(console, dbpath, vcf_files[0], contigs,
                                                            bucket_size));

    {
        // sanity check, see that we can get the contigs back
        vector<pair<string,size_t> > contigs_dbg;
        H("Reading the contigs back from DB",
          GLnexus::cli::utils::db_get_contigs(console, dbpath, contigs_dbg));
        if (contigs_dbg != contigs)
            return GLnexus::Status::Invalid("error, contigs read from DB do not match originals");
    }

    // Load the GVCFs into the database
    {
        // use an empty range filter
        vector<GLnexus::range> ranges;
        H("bulk load into DB",
          GLnexus::cli::utils::db_bulk_load(console, nr_threads, vcf_files, dbpath, ranges, contigs, false));
    }

    if (iter_compare) {
        H("compare database iteration methods",
          GLnexus::cli::utils::compare_db_itertion_algorithms(console, dbpath, 50));
    }

    // discover alleles
    // TODO: if bedfilename is empty, fill ranges with all contigs
    // TODO: overlap allele discovery with final compactions. have db_bulk_load output a RocksKeyValue pointer which we can reuse
    vector<GLnexus::range> ranges;
    H("parsing the bed file", GLnexus::cli::utils::parse_bed_file(console, bedfilename, contigs, ranges));
    GLnexus::discovered_alleles dsals;
    unsigned sample_count = 0;
    H("discover alleles",
      GLnexus::cli::utils::discover_alleles(console, nr_threads, dbpath, ranges, contigs, dsals, sample_count));
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
    vector<string> hdr_lines = { ("##GLnexusConfig="+config_name), ("##GLnexusConfigCRC32C="+cfg_crc32c) };
    auto DX_JOB_ID = std::getenv("DX_JOB_ID");
    if (DX_JOB_ID) {
        // if running in DNAnexus, record job ID in header
        hdr_lines.push_back(string("##DX_JOB_ID=")+DX_JOB_ID);
    }
    string outfile("-");
    H("Genotyping",
      GLnexus::cli::utils::genotype(console, nr_threads, dbpath, genotyper_cfg, sites, hdr_lines, outfile));

    return 0;
}


void help(const char* prog) {
    cout << "Usage: " << prog << " [options] /vcf/file/1 .. /vcf/file/N" << endl
         << "Merge and joint-call input gVCF files, emitting multi-sample BCF on" << endl
         << "standard output." << endl << endl
         << "Options:" << endl
         << "  --help, -h           print this help message" << endl
         << "  --bed FILE, -b FILE  three-column BED file of ranges to analyze (required)" << endl
         << "  --config X, -c X     configuration preset name or .yml filename (default: gatk)" << endl
         << "  --list, -l           given files contain lists of gVCF filenames, one per line" << endl
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
    console->info("glnexus_cli {} {}", GIT_REVISION, __TIMESTAMP__);

    if (argc < 2) {
        help(argv[0]);
        return 1;
    }

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"bed", required_argument, 0, 'b'},
        {"config", required_argument, 0, 'c'},
        {"list", no_argument, 0, 'l'},
        {"bucket_size", required_argument, 0, 'x'},
        {"debug", no_argument, 0, 'd'},
        {"iter_compare", no_argument, 0, 'i'},
        {0, 0, 0, 0}
    };

    int c;
    string config_name = "gatk";
    bool list_of_files = false;
    bool debug = false;
    bool iter_compare = false;
    string bedfilename;
    int nr_threads = std::thread::hardware_concurrency();
    size_t bucket_size = GLnexus::BCFKeyValueData::default_bucket_size;

    while (-1 != (c = getopt_long(argc, argv, "hb:dIx:",
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

            case 'd':
                debug = true;
                break;

            case 'h':
            case '?':
                help(argv[0]);
                exit(1);
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

    return all_steps(vcf_files, bedfilename, config_name, nr_threads, debug, iter_compare, bucket_size);
}

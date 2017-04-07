// GLnexus crude command-line interface. This is a temporary thing to get us
// bootstrapped with the core algorithms and storage engine before engineering
// an always-on "server"

#include <iostream>
#include <exception>
#include <fstream>
#include <getopt.h>
#include <sstream>
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
        console->error() << "Failed to " << desc << ": " << s.str(); \
        return 1; \
    }

// Perform all the separate GLnexus operations in one go.
// return 0 on success, 1 on failure.
static int all_steps(const vector<string> &vcf_files,
                     const string &bedfilename,
                     int nr_threads,
                     bool debug,
                     bool iter_compare,
                     size_t bucket_size) {
    GLnexus::Status s;
    GLnexus::unifier_config unifier_cfg;
    GLnexus::genotyper_config genotyper_cfg;

    if (vcf_files.empty()) {
        console->error() << "No source GVCF files specified";
        return 1;
    }

    string config_preset = "test";
    if (config_preset.size()) {
        H("load unifier/genotyper configuration",
          GLnexus::cli::utils::load_config_preset(console, config_preset, unifier_cfg, genotyper_cfg));
    }

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
          GLnexus::cli::utils::db_bulk_load(console, nr_threads, vcf_files, dbpath, ranges, contigs, true));
    }

    if (iter_compare) {
        H("compare database iteration methods",
          GLnexus::cli::utils::compare_db_itertion_algorithms(console, dbpath, 50));
    }

    // discover alleles
    vector<GLnexus::range> ranges;
    H("parsing the bed file", GLnexus::cli::utils::parse_bed_file(console, bedfilename, contigs, ranges));
    GLnexus::discovered_alleles dsals;
    unsigned sample_count = 0;
    H("discover alleles",
      GLnexus::cli::utils::discover_alleles(console, nr_threads, dbpath, ranges, contigs, dsals, sample_count));
    if (debug) {
        string filename("/tmp/dsals.yml");
        console->info() << "Writing discovered alleles as YAML to " << filename;
        H("serialize discovered alleles to a file",
          GLnexus::cli::utils::yaml_write_discovered_alleles_to_file(dsals, contigs, sample_count, filename));

        string filename_cflat("/tmp/dsals.cflat");
        console->info() << "Writing discovered alleles in serialized capnproto to " << filename
                        << ", verifying that it is readable";
        H("Verify that the cflat file is readable",
          GLnexus::capnp_discover_alleles_verify(sample_count, contigs, dsals, filename_cflat));
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
    for (auto& dsals_i : dsals_by_contig) {
        H("unify sites",
          GLnexus::cli::utils::unify_sites(console, unifier_cfg, ranges, contigs, dsals_i, sample_count, sites));
    }
    console->info() << "unified to " << sites.size() << " sites";
    if (debug) {
        string filename("/tmp/sites.yml");
        console->info() << "Writing unified sites as YAML to " << filename;
        H("write unified sites to file",
          GLnexus::cli::utils::write_unified_sites_to_file(sites, contigs, filename));
    }

    // genotype
    genotyper_cfg.output_residuals = debug;
    string outfile("-");
    H("Genotyping",
      GLnexus::cli::utils::genotype(console, nr_threads, dbpath, genotyper_cfg, sites, outfile));

    return 0;
}


void help(const char* prog) {
    cerr << "usage: " << prog << " [options] /vcf/file/1 .. /vcf/file/N" << endl
         << "Joint genotype all source VCF files, and generate a project VCF file" << endl
         << "on standard out. The source files must be in GVCF format." << endl
         << "Options:" << endl
         << "  --help, -h           print this help message" << endl
         << "  --bed FILE, -b FILE  path to three-column BED file" << endl
         << "  --bucket_size INT, -x INT  set the bucket size" << endl
         << "  --debug, -d          create additional file outputs for diagnostics/debugging" << endl
         << "  --iter_compare, -i   compare different implementations of database iteration" << endl;
}

// Expected usage:
//    glnexus [vcf files]
//
int main(int argc, char *argv[]) {
    GLnexus::Status s;
    spdlog::set_level(spdlog::level::info);
    spdlog::set_pattern("[%t] %+");
    console->info() << "glnexus_cli " << GIT_REVISION;

    if (argc <= 2) {
        help(argv[0]);
        return 1;
    }

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"bed", required_argument, 0, 'b'},
        {"bucket_size", required_argument, 0, 'x'},
        {"debug", no_argument, 0, 'd'},
        {"iter_compare", no_argument, 0, 'i'},
        {0, 0, 0, 0}
    };

    int c;
    bool debug = false;
    bool iter_compare = false;
    string bedfilename;
    int nr_threads = std::thread::hardware_concurrency();
    size_t bucket_size = GLnexus::BCFKeyValueData::default_bucket_size;

    optind = 1; // force optind past command positional argument
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

    vector<string> vcf_files;
    for (int i=optind; i < argc; i++)
        vcf_files.push_back(string(argv[i]));

    return all_steps(vcf_files, bedfilename, nr_threads, debug, iter_compare, bucket_size);
}

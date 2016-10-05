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
#include "compare_iter.h"
#include "cli_utils.h"

using namespace std;
using namespace GLnexus::cli;

auto console = spdlog::stderr_logger_mt("GLnexus");

GLnexus::Status parse_bed_file(const string &bedfilename,
                               const std::vector<std::pair<std::string,size_t> > &contigs,
                               vector<GLnexus::range> &ranges) {
    if (bedfilename.empty()) {
        return GLnexus::Status::Invalid("Empty bed file");
    }
    if (!utils::check_file_exists(bedfilename)) {
        return GLnexus::Status::IOError("bed file does not exist", bedfilename);
    }

    // read BED file
    string rname, beg_txt, end_txt;
    ifstream bedfile(bedfilename);
    while (bedfile >> rname >> beg_txt >> end_txt) {
        int rid = 0;
        for(; rid<contigs.size(); rid++)
            if (contigs[rid].first == rname)
                break;
        if (rid == contigs.size()) {
            return GLnexus::Status::Invalid("Unknown contig ", rname);
        }
        ranges.push_back(GLnexus::range(rid,
                                        strtol(beg_txt.c_str(), nullptr, 10),
                                        strtol(end_txt.c_str(), nullptr, 10)));
    }
    if (bedfile.bad() || !bedfile.eof()) {
        return GLnexus::Status::IOError( "Error reading ", bedfilename);
    }

    sort(ranges.begin(), ranges.end());
    for (auto& query : ranges) {
        if (query.beg < 0 || query.end < 1 || query.end <= query.beg) {
            return GLnexus::Status::Invalid("query range ", query.str(contigs));
        }
        if (query.end > contigs[query.rid].second) {
            query.end = contigs[query.rid].second;
            console->warn() << "Truncated query range at end of contig: " << query.str(contigs);
        }
    }

    // make an orderly pass on the ranges, and verify that there are no overlaps.
    if (ranges.size() > 1) {
        const GLnexus::range &prev = *(ranges.begin());
        for (auto it = ranges.begin() + 1; it != ranges.end(); ++it) {
            if (prev.overlaps(*it))
                return GLnexus::Status::Invalid("overlapping ranges ", prev.str(contigs) + " " + it->str(contigs));
        }
    }

    return GLnexus::Status::OK();
}




// Perform all the separate GLnexus operations in one go.
GLnexus::Status all_steps(const vector<string> &vcf_files,
                          const string &bedfilename,
                          int nr_threads,
                          bool residuals,
                          bool debug) {
    GLnexus::Status s;
    GLnexus::unifier_config unifier_cfg;
    GLnexus::genotyper_config genotyper_cfg;

    if (vcf_files.empty())
        return GLnexus::Status::Invalid("No source GVCF files specified");

    string config_preset = "test";
    if (config_preset.size()) {
        S(GLnexus::cli::utils::load_config_preset(console, config_preset, unifier_cfg, genotyper_cfg));
    }

    // initilize empty database
    string dbpath("GLnexus.DB");
    vector<pair<string,size_t> > contigs;
    S(utils::recursive_delete(dbpath));
    S(utils::db_init(console, dbpath, vcf_files[0], contigs));

    // Load the GVCFs into the database
    {
        // do not use a range filter
        vector<GLnexus::range> ranges;
        S(GLnexus::cli::utils::db_bulk_load(console, nr_threads, vcf_files, dbpath, ranges, contigs, true));
    }

    // discover alleles
    vector<GLnexus::range> ranges;
    S(parse_bed_file(bedfilename, contigs, ranges));
    GLnexus::discovered_alleles dsals;
    unsigned sample_count = 0;
    S(GLnexus::cli::utils::discover_alleles(console, nr_threads, dbpath, ranges, contigs, dsals, sample_count));
    if (debug) {
        string filename("/tmp/dsals.yml");
        S(GLnexus::cli::utils::yaml_write_discovered_alleles_to_file(dsals, contigs, sample_count, filename));
    }

    // unify sites
    vector<GLnexus::unified_site> sites;
    S(GLnexus::cli::utils::unify_sites(console, unifier_cfg, ranges, contigs, dsals, sample_count, sites));
    if (debug) {
        string filename("/tmp/sites.yml");
        S(GLnexus::cli::utils::write_unified_sites_to_file(sites, contigs, filename));
    }

    // genotype
    genotyper_cfg.output_residuals = residuals;
    string outfile("-");
    S(GLnexus::cli::utils::genotype(console, nr_threads, dbpath, genotyper_cfg, sites, outfile));

    return GLnexus::Status::OK();
}


void help(const char* prog) {
    cerr << "usage: " << prog << " [options] /vcf/file/1 .. /vcf/file/N" << endl
         << "Joint genotype all source VCF files, and generate a project VCF file" << endl
         << "on standard out. The source files must be in GVCF format." << endl
         << "Options:" << endl
         << "  --help, -h       print this help message" << endl
         << "  --bed FILE, -b FILE  path to three-column BED file" << endl
         << "  --debug, -d      create additional file outputs for diagnostics/debugging" << endl
         << "  --residuals, -r  generate detailed residuals output file" << endl;
}

// Expected usage:
//    glnexus_aio [vcf files]
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
        {"debug", no_argument, 0, 'd'},
        {"residuals", no_argument, 0, 'r'},
        {0, 0, 0, 0}
    };

    int c;
    bool debug = false;
    bool residuals = false;
    string bedfilename;
    int nr_threads = std::thread::hardware_concurrency();

    optind = 1; // force optind past command positional argument
    while (-1 != (c = getopt_long(argc, argv, "hb:dr",
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

            case 'r':
                residuals = true;
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

    s = all_steps(vcf_files, bedfilename, nr_threads, residuals, debug);
    if (s.bad()) {
        cerr << s.str() << endl;
        return 1;
    }
    return 0;
}

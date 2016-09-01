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
GLnexus::Status s;
#define H(desc,expr) \
    s = expr; \
    if (s.bad()) { \
        console->error() << "Failed to " << desc << ": " << s.str(); \
        return 1; \
    }

static GLnexus::RocksKeyValue::prefix_spec* GLnexus_prefix_spec() {
    static unique_ptr<GLnexus::RocksKeyValue::prefix_spec> p;
    if (!p) {
        p = make_unique<GLnexus::RocksKeyValue::prefix_spec>("bcf", GLnexus::BCFKeyValueDataPrefixLength());
    }
    return p.get();
}

// hard-coded configuration presets for unifier & genotyper. TODO: these
// should reside in some user-modifiable yml file
static const char* config_presets_yml = R"eof(
unifier_config:
  min_AQ1: 70
  min_AQ2: 40
  min_GQ: 70
genotyper_config:
  required_dp: 1
  liftover_fields:
    - orig_names: [GQ]
      name: GQ
      description: '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">'
      type: int
      number: basic
      combi_method: min
      count: 1
      ignore_non_variants: true
    - orig_names: [DP, MIN_DP]
      name: DP
      description: '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">'
      type: int
      combi_method: min
      number: basic
      count: 1
    - orig_names: [AD]
      name: AD
      description: '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">'
      type: int
      number: alleles
      combi_method: min
      default_type: zero
      count: 0
    - orig_names: [SB]
      name: SB
      description: '##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fishers Exact Test to detect strand bias.">'
      type: int
      combi_method: max
      number: basic
      count: 4
)eof";

static GLnexus::Status load_config_preset(const std::string& name,
                                          GLnexus::unifier_config& unifier_cfg,
                                          GLnexus::genotyper_config& genotyper_cfg) {
    GLnexus::Status s;
    cerr << "Loading config "<< endl << config_presets_yml;
    YAML::Node yaml = YAML::Load(config_presets_yml);
    if (!yaml) {
        return GLnexus::Status::NotFound("unknown configuration preset", name);
    }
    if (!yaml.IsMap()) {
        return GLnexus::Status::Invalid("configuration presets");
    }
    if (yaml["unifier_config"]) {
        S(GLnexus::unifier_config::of_yaml(yaml["unifier_config"], unifier_cfg));
    }
    if (yaml["genotyper_config"]) {
        S(GLnexus::genotyper_config::of_yaml(yaml["genotyper_config"], genotyper_cfg));
    }

    return GLnexus::Status::OK();
}


void help_load(const char* prog) {
    cerr << "usage: " << prog << " load [options] /db/path sample.gvcf[.gz] [sample2.gvcf[.gz] ...]" << endl
         << "Loads gVCF file(s) into an existing database. The data set name will be derived from" << endl
         << "the gVCF filename. It can be overridden with --dataset if loading only one gVCF." << endl
         << "If the final argument is - then gVCF filenames are read from standard input." << endl
         << "Options:" << endl
         << "  --range, -r chr1,chr2:1000-2000  load only records overlapping a given range" << endl
         << "  --delete, -X                     delete each gVCF file immediately after successful load" << endl
         << "  --threads N, -t N                override thread pool size (default: nproc)" << endl
         << endl;
}

int main_load(int argc, char *argv[]) {
    if (argc == 2) {
        help_load(argv[0]);
        return 1;
    }

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"dataset", required_argument, 0, 'd'},
        {"range", required_argument, 0, 'r'},
        {"and-delete", no_argument, 0, 'X'},
        {"threads", required_argument, 0, 't'},
        {0, 0, 0, 0}
    };

    string dataset, ranges_txt;
    bool and_delete = false;
    size_t threads = std::thread::hardware_concurrency();

    int c;
    optind = 2; // force optind past command positional argument
    while (-1 != (c = getopt_long(argc, argv, "hd:r:Xt:",
                                  long_options, nullptr))) {
        switch (c) {
            case 'd':
                dataset = string(optarg);
                if (dataset.size() == 0) {
                    cerr <<  "invalid --dataset" << endl;
                    return 1;
                }
                break;

            case 'r':
                ranges_txt = string(optarg);
                break;

            case 'X':
                and_delete = true;
                break;

            case 't':
                threads = strtoul(optarg, nullptr, 10);
                break;

            case 'h':
            case '?':
                help_load(argv[0]);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    if (argc-optind < 2) {
        help_load(argv[0]);
        return 1;
    }
    string dbpath(argv[optind]);
    vector<string> gvcfs;

    for (size_t i = optind+1; i < argc; i++) {
        gvcfs.push_back(string(argv[i]));
    }

    if (!dataset.empty() && gvcfs.size() > 1) {
        cerr << "--dataset applicable only when loading exactly one gVCF" << endl;
        return 1;
    }

    if (gvcfs[gvcfs.size()-1] == "-") {
        // read list of filenames from stdin
        gvcfs.erase(gvcfs.end()-1);
        for (string fn; std::getline(cin, fn);) {
            gvcfs.push_back(fn);
        }
    }

    if (gvcfs.empty()) {
        console->warn() << "Nothing to do";
        return 0;
    }

    // open the database
    unique_ptr<GLnexus::KeyValue::DB> db;
    H("open database", GLnexus::RocksKeyValue::Open(dbpath, db, GLnexus_prefix_spec(),
                                                    GLnexus::RocksKeyValue::OpenMode::BULK_LOAD));

    {
        unique_ptr<GLnexus::BCFKeyValueData> data;
        H("open database", GLnexus::BCFKeyValueData::Open(db.get(), data));

        {
            unique_ptr<GLnexus::MetadataCache> metadata;
            H("instantiate metadata cache", GLnexus::MetadataCache::Start(*data, metadata));
            const auto& contigs = metadata->contigs();

            set<GLnexus::range> ranges;
            if (!ranges_txt.empty()) {
                vector<GLnexus::range> vranges;
                if (!utils::parse_ranges(contigs, ranges_txt, vranges) || vranges.empty()) {
                    console->error() << "empty, invalid, or out-of-bound range(s): " << ranges_txt;
                    return 1;
                }
                ranges.insert(vranges.begin(), vranges.end());
                ostringstream ss;
                for (const auto& rng : ranges) {
                    ss << " " << rng.str(contigs);
                }
                console->info() << "Beginning bulk load of records overlapping:" << ss.str();
            } else {
                console->info() << "Beginning bulk load with no range filter.";
            }

            ctpl::thread_pool threadpool(threads);
            vector<future<GLnexus::Status>> statuses;
            set<string> datasets_loaded;
            GLnexus::BCFKeyValueData::import_result stats;
            mutex mu;

            // load the gVCFs on the thread pool
            for (const auto& gvcf : gvcfs) {
                // default dataset name (the gVCF filename)
                if (dataset.empty()) {
                    size_t p = gvcf.find_last_of('/');
                    if (p != string::npos && p < gvcf.size()-1) {
                        dataset = gvcf.substr(p+1);
                    } else {
                        dataset = gvcf;
                    }
                }


                auto fut = threadpool.push([&, gvcf, dataset](int tid){
                    GLnexus::BCFKeyValueData::import_result rslt;
                    GLnexus::Status ls = data->import_gvcf(*metadata, dataset, gvcf, ranges, rslt);
                    if (ls.ok()) {
                        if (and_delete && unlink(gvcf.c_str())) {
                            console->warn() << "Loaded " << gvcf << " successfully, but failed deleting it afterwards.";
                        }
                        lock_guard<mutex> lock(mu);
                        stats += rslt;
                        datasets_loaded.insert(dataset);
                        size_t n = datasets_loaded.size();
                        if (n % 100 == 0) {
                            console->info() << n << "...";
                        }
                    }
                    return ls;
                });
                statuses.push_back(move(fut));
                dataset.clear();
            }

            // collect results
            vector<pair<string,GLnexus::Status>> failures;
            for (size_t i = 0; i < gvcfs.size(); i++) {
                GLnexus::Status s_i(move(statuses[i].get()));
                if (!s_i.ok()) {
                    failures.push_back(make_pair(gvcfs[i],move(s_i)));
                }
            }

            // report results
            console->info() << "Loaded " << datasets_loaded.size() << " datasets with "
                            << stats.samples.size() << " samples; "
                            << stats.bytes << " bytes in "
                            << stats.records << " BCF records in "
                            << stats.buckets << " buckets. "
                            << "Bucket max " << stats.max_bytes << " bytes, max "
                            << stats.max_records << " records.";

            // call all_samples_sampleset to create the sample set including
            // the newly loaded ones. By doing this now we make it possible
            // for other CLI functions to open the database in purely read-
            // only mode (since the sample set has to get written into the
            // database to be used)
            string sampleset;
            H("update * sample set", data->all_samples_sampleset(sampleset));
            console->info() << "Created sample set " << sampleset;

            if (failures.size()) {
                console->error() << "FAILED to load " << failures.size() << " datasets:";
                for (const auto& p : failures) {
                    console->error() << p.first << " " << p.second.str();
                }
                return datasets_loaded.size() ? 2 : 1;
            }
        }
    }

    console->info() << "Flushing and compacting database...";
    H("final database flush", db->flush());
    db.reset();
    console->info() << "Bulk load complete!";

    return 0;
}



void help_discover_alleles(const char* prog) {
    cerr << "usage: " << prog << " discover [options] /db/path chrom:1234-2345" << endl
         << "Discover alleles in all samples in the database in the given interval. The positions" << endl
         << "are one-based, inclusive. As an alternative to providing one interval on the" << endl
         << "command line, you can provide a three-column BED file using --bed." << endl
         << "Options:" << endl
         << "  --bed FILE, -b FILE  path to three-column BED file" << endl
         << "  --config X, -c X     apply unifier/genotyper configuration preset X" << endl
         << "  --threads N, -t N    override thread pool size (default: nproc)" << endl
         << endl;
}

int main_discover_alleles(int argc, char *argv[]) {
    if (argc == 2) {
        help_discover_alleles(argv[0]);
        return 1;
    }

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"bed", required_argument, 0, 'b'},
        {"config", required_argument, 0, 'c'},
        {"threads", required_argument, 0, 't'},
        {0, 0, 0, 0}
    };

    string bedfilename;
    string config_preset;
    size_t threads = 0;

    int c;
    optind = 2; // force optind past command positional argument
    while (-1 != (c = getopt_long(argc, argv, "hb:c:t:",
                                  long_options, nullptr))) {
        switch (c) {
            case 'b':
                bedfilename = string(optarg);
                if (bedfilename.size() == 0) {
                    cerr <<  "invalid BED filename" << endl;
                    return 1;
                }
                break;
            case 'c':
                config_preset = string(optarg);
                break;
            case 't':
                threads = strtoul(optarg, nullptr, 10);
                console->info() << "pooling " << threads << " threads for genotyping";
                break;
            case 'h':
            case '?':
                help_discover_alleles(argv[0]);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    if (optind != argc-1) {
        help_discover_alleles(argv[0]);
        return 1;
    }
    string dbpath(argv[optind]);
    string range_txt;

    if (bedfilename.empty()) {
        if (optind != argc-2) {
            help_discover_alleles(argv[0]);
            return 1;
        }

        range_txt = argv[optind+1];
    } else {
        if (optind != argc-1) {
            help_discover_alleles(argv[0]);
            return 1;
        }
    }

    GLnexus::unifier_config unifier_cfg;
    GLnexus::genotyper_config genotyper_cfg;
    if (config_preset.size()) {
        H("load configuration preset", load_config_preset(config_preset, unifier_cfg, genotyper_cfg));
    }

    unique_ptr<GLnexus::KeyValue::DB> db;
    unique_ptr<GLnexus::BCFKeyValueData> data;

    // open the database in read-only mode
    H("open database", GLnexus::RocksKeyValue::Open(dbpath, db, GLnexus_prefix_spec(),
                                                    GLnexus::RocksKeyValue::OpenMode::READ_ONLY));
    H("open database", GLnexus::BCFKeyValueData::Open(db.get(), data));

    std::vector<std::pair<std::string,size_t> > contigs;
    H("read contig metadata", data->contigs(contigs));

    vector<GLnexus::range> ranges;
    H("Parsing range argument", parse_ranges_from_cmd_line(bedfilename, range_txt, contigs, ranges));

    // start service, discover alleles
    GLnexus::service_config svccfg;
    svccfg.threads = threads;
    unique_ptr<GLnexus::Service> svc;
    H("start GLnexus service", GLnexus::Service::Start(svccfg, *data, *data, svc));

    string sampleset;
    H("list all samples", data->all_samples_sampleset(sampleset));
    console->info() << "found sample set " << sampleset;

    console->info() << "discovering alleles in " << ranges.size() << " range(s)";
    vector<GLnexus::discovered_alleles> valleles;
    unsigned N;
    H("discover alleles", svc->discover_alleles(sampleset, ranges, N, valleles));

    GLnexus::discovered_alleles dsals;
    for (auto it = valleles.begin(); it != valleles.end(); ++it) {
        H("merge vectors", merge_discovered_alleles(*it, dsals));
    }
    console->info() << "discovered " << dsals.size() << " alleles";

    // Write the discovered alleles to stdout
    H("write results as yaml", utils::yaml_stream_of_discovered_alleles(N, contigs, dsals, cout));
    return 0;
}


void help_unify_sites(const char* prog) {
    cerr << "usage: " << prog << " unify_sites [options] /discovered_alleles/path/1 ... /discovered_alleles/path/N" << endl
         << "Unify the discovered allele file(s). There can be one or more files, format is YAML." << endl
         << "Options:" << endl
         << "  --bed FILE, -b FILE  path to three-column BED file, the ranges have to the same same as the discover-alleles phase" << endl
         << "  --config X, -c X     apply unifier/genotyper configuration preset X" << endl
         << "  --threads N, -t N    override thread pool size (default: nproc)" << endl
         << endl;
}


int main_unify_sites(int argc, char *argv[]) {
    if (argc == 2) {
        help_unify_sites(argv[0]);
        return 1;
    }

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"bed", required_argument, 0, 'b'},
        {"config", required_argument, 0, 'c'},
        {"threads", required_argument, 0, 't'},
        {0, 0, 0, 0}
    };

    string bedfilename;
    string config_preset;
    size_t threads = 0;

    int c;
    optind = 2; // force optind past command positional argument
    while (-1 != (c = getopt_long(argc, argv, "hb:c:",
                                  long_options, nullptr))) {
        switch (c) {
            case 'b':
                bedfilename = string(optarg);
                if (bedfilename.size() == 0) {
                    cerr <<  "invalid BED filename" << endl;
                    return 1;
                }
                break;
            case 'c':
                config_preset = string(optarg);
                break;
            case 't':
                threads = strtoul(optarg, nullptr, 10);
                console->info() << "pooling " << threads << " threads for site unification";
                break;
            case 'h':
            case '?':
                help_unify_sites(argv[0]);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    if (optind > argc-1) {
        help_unify_sites(argv[0]);
        return 1;
    }
    vector<string> discovered_allele_files;
    for (int i = optind; i < argc; i++) {
        discovered_allele_files.push_back(argv[i]);
    }
    if (discovered_allele_files.size() == 0) {
        help_unify_sites(argv[0]);
        return 1;
    }

    GLnexus::unifier_config unifier_cfg;
    GLnexus::genotyper_config genotyper_cfg;
    if (config_preset.size()) {
        H("load configuration preset", load_config_preset(config_preset, unifier_cfg, genotyper_cfg));
    }

    unsigned N;
    vector<pair<string,size_t> > contigs;
    GLnexus::discovered_alleles dsals;
    H("merge discovered allele files",
      utils::merge_discovered_allele_files(console, threads, discovered_allele_files, N, contigs, dsals));
    console->info() << "consolidated to " << dsals.size() << " alleles for site unification; N = " << N;

    set<GLnexus::range> ranges;
    if (!bedfilename.empty()) {
        string range_txt; // empty string, the plan is to get rid of the argument
        vector<GLnexus::range> i_ranges;
        H("Parsing range argument", parse_ranges_from_cmd_line(bedfilename, range_txt, contigs, i_ranges));
        for (auto &r : i_ranges)
            ranges.insert(r);
    }

    vector<GLnexus::unified_site> sites;
    H("unify sites", GLnexus::unified_sites(unifier_cfg, N, dsals, sites));
    // sanity check, sites are in-order and non-overlapping
    if (sites.size() > 1) {
        auto p = sites.begin();
        for (auto q = p+1; q != sites.end(); ++p, ++q) {
            if (!(p->pos < q->pos) || p->pos.overlaps(q->pos)) {
                console->error() << "BUG: unified sites failed sanity check -- sites are out of order or overlapping. "
                                 << p->pos.str(contigs) << " & " << q->pos.str(contigs);
                return 1;
            }
        }
    }
    console->info() << "unified to " << sites.size() << " sites";

    if (!ranges.empty()) {
        // set the containing ranges for each site
        for (auto &us : sites) {
            H("find target range containing site " + us.pos.str(contigs),
              utils::find_containing_range(ranges, us.pos, us.containing_target));
        }
    }

    // Write the unified sites to stdout
    H("write results as yaml", utils::yaml_stream_of_unified_sites(sites, contigs, cout));

    return 0;
}

void help_genotype(const char* prog) {
    cerr << "usage: " << prog << " genotype [options] /db/path /unified_sites/path" << endl
         << "Genotype all samples in the database. The unified sites" << endl
         << "should be provided as a file in yaml format. The positions are" << endl
         << "one-based, inclusive." << endl
         << "Options:" << endl
         << "  --residuals, -r      generate detailed residuals output file" << endl
         << "  --config X, -c X     apply unifier/genotyper configuration preset X" << endl
         << "  --threads N, -t N    override thread pool size (default: nproc)" << endl
         << endl;
}

int main_genotype(int argc, char *argv[]) {
    if (argc == 2) {
        help_genotype(argv[0]);
        return 1;
    }

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"residuals", no_argument, 0, 'r'},
        {"config", required_argument, 0, 'c'},
        {"threads", required_argument, 0, 't'},
        {0, 0, 0, 0}
    };

    string config_preset;
    bool residuals = false;
    size_t threads = 0;

    int c;
    optind = 2; // force optind past command positional argument
    while (-1 != (c = getopt_long(argc, argv, "hrc:t:",
                                  long_options, nullptr))) {
        switch (c) {
            case 'r':
                residuals = true;
                break;
            case 'c':
                config_preset = string(optarg);
                break;
            case 't':
                threads = strtoul(optarg, nullptr, 10);
                console->info() << "pooling " << threads << " threads for genotyping";
                break;
            case 'h':
            case '?':
                help_genotype(argv[0]);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    if (optind > argc-1) {
        help_genotype(argv[0]);
        return 1;
    }
    string dbpath(argv[optind]);
    string unified_sites_file(argv[optind + 1]);
    if (optind != argc-2) {
        help_genotype(argv[0]);
        return 1;
    }

    GLnexus::unifier_config unifier_cfg;
    GLnexus::genotyper_config genotyper_cfg;
    if (config_preset.size()) {
        H("load configuration preset", load_config_preset(config_preset, unifier_cfg, genotyper_cfg));
    }

    genotyper_cfg.output_residuals = residuals;
    console->info() << "Lifting over " << genotyper_cfg.liftover_fields.size() << " fields.";
    unique_ptr<GLnexus::KeyValue::DB> db;

    // open the database in read-only mode
    H("open database", GLnexus::RocksKeyValue::Open(dbpath, db, GLnexus_prefix_spec(),
                                                    GLnexus::RocksKeyValue::OpenMode::READ_ONLY));
    {
        unique_ptr<GLnexus::BCFKeyValueData> data;
        H("open database", GLnexus::BCFKeyValueData::Open(db.get(), data));

        std::vector<std::pair<std::string,size_t> > contigs;
        H("read contig metadata", data->contigs(contigs));

        vector<GLnexus::unified_site> sites;
        {
            std::ifstream ifs;
            ifs.open(unified_sites_file, std::ifstream::in);
            if (!ifs.good()) {
                console->info() << "Error opening file " << unified_sites_file;
                return 1;
            }
            H("load unified sites",
              utils::unified_sites_of_yaml_stream(ifs, contigs, sites));
            ifs.close();
        }
        console->info() << "loaded " << sites.size() << " sites from " << unified_sites_file;

        // start service, discover alleles, unify sites, genotype sites
        GLnexus::service_config svccfg;
        svccfg.threads = threads;
        unique_ptr<GLnexus::Service> svc;
        H("start GLnexus service", GLnexus::Service::Start(svccfg, *data, *data, svc));

        string sampleset;
        H("list all samples", data->all_samples_sampleset(sampleset));
        console->info() << "found sample set " << sampleset;

        H("genotype sites",
          svc->genotype_sites(genotyper_cfg, sampleset, sites, string("-")));
        console->info("genotyping complete!");

        auto stalls_ms = svc->threads_stalled_ms();
        if (stalls_ms) {
            console->info() << "worker threads were cumulatively stalled for " << stalls_ms << "ms";
        }

        std::shared_ptr<GLnexus::StatsRangeQuery> statsRq = data->getRangeStats();
        cerr << statsRq->str() << endl;
    }

    return 0;
}

void help_iter_compare(const char* prog) {
    cerr << "usage: " << prog << " iter_compare /db/path" << endl
         << "Run tests comparing the two BCF iterators" << endl
         << endl;
}

int main_iter_compare(int argc, char *argv[]) {
    if (argc != 3) {
        help_iter_compare(argv[0]);
        return 1;
    }
    string dbpath(argv[2]);

    unique_ptr<GLnexus::KeyValue::DB> db;
    unique_ptr<GLnexus::BCFKeyValueData> data;
    string sampleset;
    H("open database", GLnexus::RocksKeyValue::Open(dbpath, db, GLnexus_prefix_spec(),
                                                    GLnexus::RocksKeyValue::OpenMode::READ_ONLY));
    H("open database", GLnexus::BCFKeyValueData::Open(db.get(), data));

    unique_ptr<GLnexus::MetadataCache> metadata;
    H("instantiate metadata cache", GLnexus::MetadataCache::Start(*data, metadata));

    H("all_samples_sampleset", data->all_samples_sampleset(sampleset));
    console->info() << "using sample set " << sampleset;

    const auto& contigs = metadata->contigs();

    // get samples and datasets
    shared_ptr<const set<string>> samples, datasets;
    H("sampleset_datasets", metadata->sampleset_datasets(sampleset, samples, datasets));

    int nChroms = min((size_t)22, contigs.size());
    int nIter = 50;
    int maxRangeLen = 1000000;
    int minLen = 10; // ensure that the the range is of some minimal size

    for (int i = 0; i < nIter; i++) {
        int rid = genRandNumber(nChroms);
        int lenChrom = (int)contigs[rid].second;
        assert(lenChrom > minLen);

        // bound the range to be no larger than the chromosome
        int rangeLen = min(maxRangeLen, lenChrom);
        assert(rangeLen > minLen);

        int beg = genRandNumber(lenChrom - rangeLen);
        int rlen = genRandNumber(rangeLen - minLen);
        GLnexus::range rng(rid, beg, beg + minLen + rlen);

        int rc = compare_query(*data, *metadata, sampleset, rng);
        switch (rc) {
        case 1: break;
        case 0: return 1; // ERROR Status
        case -1: break;  // Query used too much memory, aborted
        }
    }

    cout << "Passed " << nIter << " iterator comparison tests" << endl;
    return 0; // GOOD STATUS
}

// Initialize the database
static Status init_db(const string &dbpath, const string &exemplar_gvcf) {
    // load exemplar contigs
    unique_ptr<vcfFile, void(*)(vcfFile*)> vcf(bcf_open(exemplar_gvcf.c_str(), "r"),
                                               [](vcfFile* f) { bcf_close(f); });
    if (!vcf) {
        return GLnexus::Status::IOError("Failed to open exemplar gVCF file at ", exemplar_gvcf);
    }
    unique_ptr<bcf_hdr_t, void(*)(bcf_hdr_t*)> hdr(bcf_hdr_read(vcf.get()), &bcf_hdr_destroy);
    if (!hdr) {
        cerr << "Failed to read gVCF file header from " << exemplar_gvcf << endl;
        return 1;
    }
    int ncontigs = 0;
    const char **contignames = bcf_hdr_seqnames(hdr.get(), &ncontigs);
    vector<pair<string,size_t>> contigs;
    for (int i = 0; i < ncontigs; i++) {
        if (hdr->id[BCF_DT_CTG][i].val == nullptr) {
            return GLnexus::Status::Invalid("Invalid gVCF header in ", exemplar_gvcf);
        }
        contigs.push_back(make_pair(string(contignames[i]),
                                    hdr->id[BCF_DT_CTG][i].val->info[0]));
    }

    // create and initialize the database
    size_t bucket_size = 30000;
    unique_ptr<GLnexus::KeyValue::DB> db;
    S(GLnexus::RocksKeyValue::Initialize(dbpath, db, GLnexus_prefix_spec()));
    S(GLnexus::BCFKeyValueData::InitializeDB(db.get(), contigs, bucket_size));

    // report success
    console->info() << "Initialized GLnexus database in " << dbpath;
    console->info() << "bucket size: " << bucket_size;

    stringstream ss;
    ss << "contigs:";
    for (const auto& contig : contigs) {
        ss << " " << contig.first;
    }
    console->info() << ss.str() << endl;

    return 0;
}

// It is equivalent to the following subcommands, issued from glnexus_cli:
//    glnexus_cli init $bucket_size_arg GLnexus.db $(find in/gvcf -type f | head -n 1)
//    cat all_gvcfs.txt | time glnexus_cli load $ranges_to_load_arg --and-delete GLnexus.db -
//
//    glnexus_cli discover_alleles GLnexus.db --bed ranges.bed $config_flag > discovered_alleles.yml
//    glnexus_cli unify_sites discovered_alleles.yml --bed ranges.bed $config_flag > unified_sites.yml
//    glnexus_cli genotype GLnexus.db unified_sites.yml $residuals_flag $config_flag | bcftools view - | $vcf_compressor -c > "out/vcf/${output_name}.vcf.${compress_ext}"
//
Status all_in_one(const vector<string> &vcf_files, const string &proj_vcf) {
    // initilize DB
    string dbpath(argv[optind]);
    string exemplar_gvcf(argv[optind+1]);
    S(init_db(dbpath, exemplar_gvcf));

    return GLnexus::Status::OK();
}


void help(const char* prog) {
    cerr << "usage: " << prog << "/project/vcf/path [options] /vcf/file/1 .. /vcf/file/N" << endl
         << "Joint genotype all source VCF files, and generate a project VCF file. The" << endl
         << "source files must be in GVCF format." << endl
         << "Options:" << endl
         << "  --help, -h      print this help message" << endl;
}

// Expected usage:
//    glnexus_cli_aio proj.vcf.gz [vcf files]
//
int main(int argc, char *argv[]) {
    spdlog::set_level(spdlog::level::info);
    spdlog::set_pattern("[%t] %+");
    console->info() << "glnexus_cli_aio " << GIT_REVISION;

    if (argc <= 2) {
        help(argv[0]);
        return 1;
    }

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    optind = 1; // force optind past command positional argument
    while (-1 != (c = getopt_long(argc, argv, "h",
                                  long_options, nullptr))) {
        switch (c) {
            case 'h':
            case '?':
                help(argv[0]);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    if (optind > argc-2) {
        help_genotype(argv[0]);
        return 1;
    }

    string proj_vcf(argv[optind]);
    vector<string> vcf_files;
    for (int i=optind+1; i < argc; i++)
        vcf_files.push_back(string(argv[i]));
    assert(!vcf_files.empty());

    S(all_in_one(vcf_files, proj_vcf));
    return 0;
}

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
      description: '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">'
      type: int
      number: basic
      combi_method: min
      count: 1
      ignore_non_variants: true
    - orig_names: [DP, MIN_DP]
      name: DP
      description: '##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">'
      type: int
      combi_method: min
      number: basic
      count: 1
    - orig_names: [AD]
      name: AD
      description: '##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">'
      type: int
      number: alleles
      combi_method: min
      default_type: zero
      count: 0
    - orig_names: [SB]
      name: SB
      description: '##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Per-sample component statistics which comprise the Fishers Exact Test to detect strand bias.\">'
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


void help_init(const char* prog) {
    cerr << "usage: " << prog << " init [options] /desired/db/path exemplar.gvcf[.gz]" << endl
         << "Initializes a new GLnexus database in the specified directory; the parent directory" << endl
         << "must exist but the directory itself must not. The exemplar gVCF file is used to set" << endl
         << "the contigs in the database configuration; it is NOT loaded into the database."
         << "Options:" << endl
         << "  --bucket-size, -b N  gVCF bucket size (default 30,000bp)" << endl
         << endl;
}

int main_init(int argc, char *argv[]) {
    if (argc == 2) {
        help_init(argv[0]);
        return 1;
    }

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"bucket-size", required_argument, 0, 'b'},
        {0, 0, 0, 0}
    };

    size_t bucket_size = 30000;

    int c;
    optind = 2; // force optind past command positional argument
    while (-1 != (c = getopt_long(argc, argv, "hb:",
                                  long_options, nullptr))) {
        switch (c) {
            case 'b':
                bucket_size = strtoul(optarg, nullptr, 10);
                if (bucket_size == 0 || bucket_size > 1000000000) {
                    cerr << "bucket size should be in (1,1e9]" << endl;
                    return 1;
                }
                break;

            case 'h':
            case '?':
                help_init(argv[0]);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    if (optind != argc-2) {
        help_init(argv[0]);
        return 1;
    }
    string dbpath(argv[optind]);
    string exemplar_gvcf(argv[optind+1]);

    // load exemplar contigs
    unique_ptr<vcfFile, void(*)(vcfFile*)> vcf(bcf_open(exemplar_gvcf.c_str(), "r"),
                                               [](vcfFile* f) { bcf_close(f); });
    if (!vcf) {
        cerr << "Failed to open exemplar gVCF file at " << exemplar_gvcf << endl;
        return 1;
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
            cerr << "Invalid gVCF header in " << exemplar_gvcf << endl;
            return 1;
        }
        contigs.push_back(make_pair(string(contignames[i]),
                                    hdr->id[BCF_DT_CTG][i].val->info[0]));
    }

    // create and initialize the database
    unique_ptr<GLnexus::KeyValue::DB> db;
    H("create database", GLnexus::RocksKeyValue::Initialize(dbpath, db, GLnexus_prefix_spec()));
    H("initialize database", GLnexus::BCFKeyValueData::InitializeDB(db.get(), contigs, bucket_size));

    // report success
    cout << "Initialized GLnexus database in " << dbpath << endl
         << "bucket size: " << bucket_size << endl
         << "contigs:";
    for (const auto& contig : contigs) {
        cout << " " << contig.first;
    }
    cout << endl;

    return 0;
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

void help_dump(const char* prog) {
    cerr << "usage: " << prog << " dump [options] /db/path chrom:1234-2345" << endl
         << "Dump all gVCF records in the database overlapping the given range. The positions"
         << "are one-based, inclusive."
         << endl;
}

int main_dump(int argc, char *argv[]) {
    if (argc == 2) {
        help_dump(argv[0]);
        return 1;
    }

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    optind = 2; // force optind past command positional argument
    while (-1 != (c = getopt_long(argc, argv, "h",
                                  long_options, nullptr))) {
        switch (c) {
            case 'h':
            case '?':
                help_dump(argv[0]);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    if (optind != argc-2) {
        help_dump(argv[0]);
        return 1;
    }
    string dbpath(argv[optind]);
    string range_txt(argv[optind+1]);

    // open the database
    unique_ptr<GLnexus::KeyValue::DB> db;
    H("open database", GLnexus::RocksKeyValue::Open(dbpath, db, GLnexus_prefix_spec(),
                                                    GLnexus::RocksKeyValue::OpenMode::READ_ONLY));

    {
        unique_ptr<GLnexus::BCFKeyValueData> data;
        H("open database", GLnexus::BCFKeyValueData::Open(db.get(), data));

        {
            unique_ptr<GLnexus::MetadataCache> metadata;
            H("instantiate metadata cache", GLnexus::MetadataCache::Start(*data, metadata));

            // parse desired range
            const auto& contigs = metadata->contigs();
            GLnexus::range query(-1,-1,-1);
            if (!utils::parse_range(contigs, range_txt, query)) {
                cerr << "Invalid range: " << range_txt << endl;
                return 1;
            }

            // query and output records
            string sampleset;
            H("list all samples", metadata->all_samples_sampleset(sampleset));
            shared_ptr<const set<string>> samples, datasets;
            H("sampleset_datasets", metadata->sampleset_datasets(sampleset, samples, datasets));
            vcfFile *vcfout = vcf_open("-", "w");
            for (const auto& dataset : *datasets) {
                shared_ptr<const bcf_hdr_t> hdr;
                vector<shared_ptr<bcf1_t>> records;
                H("dataset_range", data->dataset_range_and_header(dataset, query, 0, hdr, records));

                int bcf_nsamples = bcf_hdr_nsamples(hdr.get());
                for (const auto& record : records) {
                    cout << dataset << "\t";
                    for (int i = 0; i < bcf_nsamples; i++) {
                        if (i) cout << ",";
                        cout << bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, i);
                    }
                    cout << "\t" << flush;
                    vcf_write(vcfout, hdr.get(), record.get());
                    ignore_retval(hflush(vcfout->fp.hfile));
                }
            }
            // intentionally leaking vcfout because closing it closes stdout
        }

        std::shared_ptr<GLnexus::StatsRangeQuery> statsRq = data->getRangeStats();
        cout << statsRq->str() << endl;
    }

    return 0;
}


GLnexus::Status parse_ranges_from_cmd_line(const string &bedfilename,
                                           const string &range_txt,
                                           const std::vector<std::pair<std::string,size_t> > &contigs,
                                           vector<GLnexus::range> &ranges) {
    if (bedfilename.empty()) {
        // single range from the command line
        GLnexus::range range(-1,-1,-1);
        if (!utils::parse_range(contigs, range_txt, range)) {
            return GLnexus::Status::Invalid("range: ", range_txt);
        }
        ranges.push_back(range);
    } else {
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
    H("discover alleles", svc->discover_alleles(sampleset, ranges, valleles));

    GLnexus::discovered_alleles dsals;
    for (auto it = valleles.begin(); it != valleles.end(); ++it) {
        H("merge vectors", merge_discovered_alleles(*it, dsals));
    }
    console->info() << "discovered " << dsals.size() << " alleles";

    // Write the discovered alleles to stdout
    H("write results as yaml", utils::yaml_stream_of_discovered_alleles(contigs, dsals, cout));
    return 0;
}


void help_unify_sites(const char* prog) {
    cerr << "usage: " << prog << " unify_sites [options] /discovered_alleles/path/1 ... /discovered_alleles/path/N" << endl
         << "Unify the discovered allele file(s). There can be one or more files, format is YAML." << endl
         << "Options:" << endl
         << "  --bed FILE, -b FILE  path to three-column BED file, the ranges have to the same same as the discover-alleles phase" << endl
         << "  --config X, -c X     apply unifier/genotyper configuration preset X" << endl
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
        {0, 0, 0, 0}
    };

    string bedfilename;
    string config_preset;

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

    vector<pair<string,size_t> > contigs;
    GLnexus::discovered_alleles dsals;
    H("merge discovered allele files",
      utils::merge_discovered_allele_files(console, discovered_allele_files, contigs, dsals));
    console->info() << "consolidated to " << dsals.size() << " alleles for site unification";

    set<GLnexus::range> ranges;
    if (!bedfilename.empty()) {
        string range_txt; // empty string, the plan is to get rid of the argument
        vector<GLnexus::range> i_ranges;
        H("Parsing range argument", parse_ranges_from_cmd_line(bedfilename, range_txt, contigs, i_ranges));
        for (auto &r : i_ranges)
            ranges.insert(r);
    }

    vector<GLnexus::unified_site> sites;
    H("unify sites", GLnexus::unified_sites(unifier_cfg, dsals, sites));
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


void help(const char* prog) {
    cerr << "usage: " << prog << " <command> [options]" << endl
             << endl
             << "commands:" << endl
             << "  init             initialize new database" << endl
             << "  load             load a gVCF file into an existing database" << endl
             << "  discover_alleles discover alleles in the database  " << endl
             << "  unify_sites      unify the sites where there are alleles" << endl
             << "  genotype         genotype samples in the database" << endl
             << "  iter_compare     compare the two BCF iterator impelementations" << endl
             << endl;
}

int main(int argc, char *argv[]) {
    spdlog::set_level(spdlog::level::info);
    spdlog::set_pattern("[%t] %+");
    console->info() << "glnexus_cli " << GIT_REVISION;

    if (argc == 1) {
        help(argv[0]);
        return 1;
    }

    string command = argv[1];
    if (command == "init") {
        return main_init(argc, argv);
    } else if (command == "load") {
        return main_load(argc, argv);
    } else if (command == "dump") {
        return main_dump(argc, argv);
    } else if (command == "discover_alleles") {
        return main_discover_alleles(argc, argv);
    } else if (command == "unify_sites") {
        return main_unify_sites(argc, argv);
    } else if (command == "genotype") {
        return main_genotype(argc, argv);
    } else if (command == "iter_compare") {
        return main_iter_compare(argc, argv);
    } else {
        cerr << "unknown command " << command << endl;
        help(argv[0]);
        return 1;
    }

    return 0;
}

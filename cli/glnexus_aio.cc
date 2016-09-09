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
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

using namespace std;
using namespace GLnexus::cli;

auto console = spdlog::stderr_logger_mt("GLnexus");

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
    console->info() << "Loading config " << config_presets_yml;
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


// Initialize the database
static GLnexus::Status init(const string &dbpath,
                            const string &exemplar_gvcf,
                            vector<pair<string,size_t>> &contigs) {
    GLnexus::Status s;

    // load exemplar contigs
    unique_ptr<vcfFile, void(*)(vcfFile*)> vcf(bcf_open(exemplar_gvcf.c_str(), "r"),
                                               [](vcfFile* f) { bcf_close(f); });
    if (!vcf) {
        return GLnexus::Status::IOError("Failed to open exemplar gVCF file at ", exemplar_gvcf);
    }
    unique_ptr<bcf_hdr_t, void(*)(bcf_hdr_t*)> hdr(bcf_hdr_read(vcf.get()), &bcf_hdr_destroy);
    if (!hdr) {
        return GLnexus::Status::IOError("Failed to read gVCF file header from", exemplar_gvcf);
    }
    int ncontigs = 0;
    const char **contignames = bcf_hdr_seqnames(hdr.get(), &ncontigs);
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
    console->info() << ss.str();

    return GLnexus::Status::OK();
}


static GLnexus::Status load(const vector<string> &gvcfs,
                            const string &dbpath,
                            int nr_threads,
                            std::vector<std::pair<std::string,size_t> > &contigs) {
    GLnexus::Status s;

    // open the database
    unique_ptr<GLnexus::KeyValue::DB> db;
    S(GLnexus::RocksKeyValue::Open(dbpath, db, GLnexus_prefix_spec(),
                                   GLnexus::RocksKeyValue::OpenMode::BULK_LOAD));
    unique_ptr<GLnexus::BCFKeyValueData> data;
    S(GLnexus::BCFKeyValueData::Open(db.get(), data));

    unique_ptr<GLnexus::MetadataCache> metadata;
    S(GLnexus::MetadataCache::Start(*data, metadata));
    contigs = metadata->contigs();

    console->info() << "Beginning bulk load with no range filter.";

    ctpl::thread_pool threadpool(nr_threads);
    vector<future<GLnexus::Status>> statuses;
    set<string> datasets_loaded;
    GLnexus::BCFKeyValueData::import_result stats;
    mutex mu;
    string dataset;

    // load the gVCFs on the thread pool
    for (const auto& gvcf : gvcfs) {
        // default dataset name (the gVCF filename)
        size_t p = gvcf.find_last_of('/');
        if (p != string::npos && p < gvcf.size()-1) {
            dataset = gvcf.substr(p+1);
        } else {
            dataset = gvcf;
        }

        auto fut = threadpool.push([&, gvcf, dataset](int tid){
                set<GLnexus::range> ranges;
                GLnexus::BCFKeyValueData::import_result rslt;
                GLnexus::Status ls = data->import_gvcf(*metadata, dataset, gvcf, ranges, rslt);
                if (ls.ok()) {
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
    S(data->all_samples_sampleset(sampleset));
    console->info() << "Created sample set " << sampleset;

    if (failures.size()) {
        for (const auto& p : failures) {
            console->error() << p.first << " " << p.second.str();
        }
        return GLnexus::Status::Failure("FAILED to load ", failures.size() + " datasets:");
    }

    console->info() << "Flushing and compacting database...";
    S(db->flush());
    db.reset();
    console->info() << "Bulk load complete!";
    return GLnexus::Status::OK();
}


static GLnexus::Status discover_alleles(const string &dbpath,
                                        const vector<GLnexus::range> &ranges,
                                        const std::vector<std::pair<std::string,size_t> > &contigs,
                                        int nr_threads,
                                        GLnexus::discovered_alleles &dsals,
                                        unsigned &sample_count) {
    GLnexus::Status s;
    unique_ptr<GLnexus::KeyValue::DB> db;
    unique_ptr<GLnexus::BCFKeyValueData> data;

    // open the database in read-only mode
    S(GLnexus::RocksKeyValue::Open(dbpath, db, GLnexus_prefix_spec(),
                                   GLnexus::RocksKeyValue::OpenMode::READ_ONLY));
    S(GLnexus::BCFKeyValueData::Open(db.get(), data));

    // start service, discover alleles
    GLnexus::service_config svccfg;
    svccfg.threads = nr_threads;
    unique_ptr<GLnexus::Service> svc;
    S(GLnexus::Service::Start(svccfg, *data, *data, svc));

    string sampleset;
    S(data->all_samples_sampleset(sampleset));
    console->info() << "found sample set " << sampleset;

    console->info() << "discovering alleles in " << ranges.size() << " range(s)";
    vector<GLnexus::discovered_alleles> valleles;
    S(svc->discover_alleles(sampleset, ranges, sample_count, valleles));

    for (auto it = valleles.begin(); it != valleles.end(); ++it) {
        S(merge_discovered_alleles(*it, dsals));
    }
    console->info() << "discovered " << dsals.size() << " alleles";
    return GLnexus::Status::OK();
}


static GLnexus::Status unify_sites(const GLnexus::unifier_config &unifier_cfg,
                                   const string &dbpath,
                                   const vector<GLnexus::range> &ranges,
                                   const std::vector<std::pair<std::string,size_t> > &contigs,
                                   int nr_threads,
                                   const GLnexus::discovered_alleles &dsals,
                                   unsigned sample_count,
                                   vector<GLnexus::unified_site> &sites) {
    GLnexus::Status s;
    sites.clear();
    S(GLnexus::unified_sites(unifier_cfg, sample_count, dsals, sites));

    // sanity check, sites are in-order and non-overlapping
    if (sites.size() > 1) {
        auto p = sites.begin();
        for (auto q = p+1; q != sites.end(); ++p, ++q) {
            if (!(p->pos < q->pos) || p->pos.overlaps(q->pos)) {
                return GLnexus::Status::Failure(
                    "BUG: unified sites failed sanity check -- sites are out of order or overlapping.",
                    p->pos.str(contigs)  + " " + q->pos.str(contigs));
            }
        }
    }
    console->info() << "unified to " << sites.size() << " sites";

    if (!ranges.empty()) {
        // convert the ranges vector to a set
        set<GLnexus::range> ranges_set;
        for (auto &r : ranges)
            ranges_set.insert(r);

        // set the containing ranges for each site
        for (auto &us : sites) {
            S(utils::find_containing_range(ranges_set, us.pos, us.containing_target));
        }
    }

    return GLnexus::Status::OK();
}

GLnexus::Status genotype(const string &dbpath,
                         const GLnexus::genotyper_config &genotyper_cfg,
                         const vector<GLnexus::unified_site> &sites,
                         int nr_threads) {
    GLnexus::Status s;
    console->info() << "Lifting over " << genotyper_cfg.liftover_fields.size() << " fields.";

    // open the database in read-only mode
    unique_ptr<GLnexus::KeyValue::DB> db;
    S(GLnexus::RocksKeyValue::Open(dbpath, db, GLnexus_prefix_spec(),
                                   GLnexus::RocksKeyValue::OpenMode::READ_ONLY));
    unique_ptr<GLnexus::BCFKeyValueData> data;
    S(GLnexus::BCFKeyValueData::Open(db.get(), data));

    std::vector<std::pair<std::string,size_t> > contigs;
    S(data->contigs(contigs));

    // start service, discover alleles, unify sites, genotype sites
    GLnexus::service_config svccfg;
    svccfg.threads = nr_threads;
    unique_ptr<GLnexus::Service> svc;
    S(GLnexus::Service::Start(svccfg, *data, *data, svc));

    string sampleset;
    S(data->all_samples_sampleset(sampleset));
    console->info() << "found sample set " << sampleset;

    S(svc->genotype_sites(genotyper_cfg, sampleset, sites, string("-")));
    console->info("genotyping complete!");

    auto stalls_ms = svc->threads_stalled_ms();
    if (stalls_ms) {
        console->info() << "worker threads were cumulatively stalled for " << stalls_ms << "ms";
    }

    std::shared_ptr<GLnexus::StatsRangeQuery> statsRq = data->getRangeStats();
    console->info() << statsRq->str();

    return GLnexus::Status::OK();
}


// Perform all the separate GLnexus operations in one go.
GLnexus::Status all_in_one(const vector<string> &vcf_files,
                  const string &bedfilename,
                  int nr_threads,
                  bool residuals) {
    GLnexus::Status s;
    GLnexus::unifier_config unifier_cfg;
    GLnexus::genotyper_config genotyper_cfg;
    string config_preset = "test";
    if (config_preset.size()) {
        S(load_config_preset(config_preset, unifier_cfg, genotyper_cfg));
    }

    if (vcf_files.empty())
        return GLnexus::Status::Invalid("No source GVCF files specified");


    // initilize empty database
    string dbpath("GLnexus.DB");
    vector<pair<string,size_t> > contigs;
    boost::filesystem::remove_all(dbpath);
    S(init(dbpath, vcf_files[0], contigs));

    vector<GLnexus::range> ranges;
    if (!bedfilename.empty()) {
        string range_txt; // empty string, the plan is to get rid of the argument
        S(parse_ranges_from_cmd_line(bedfilename, range_txt, contigs, ranges));
    }

    // Load the GVCFs into the database
    S(load(vcf_files, dbpath, nr_threads, contigs));

    // discover alleles
    GLnexus::discovered_alleles dsals;
    unsigned sample_count = 0;
    S(discover_alleles(dbpath, ranges, contigs, nr_threads, dsals, sample_count));

    // unify sites
    vector<GLnexus::unified_site> sites;
    S(unify_sites(unifier_cfg, dbpath, ranges, contigs, nr_threads, dsals, sample_count, sites));

    // genotype
    genotyper_cfg.output_residuals = residuals;
    S(genotype(dbpath, genotyper_cfg, sites, nr_threads));

    return GLnexus::Status::OK();
}


void help(const char* prog) {
    cerr << "usage: " << prog << " [options] /vcf/file/1 .. /vcf/file/N" << endl
         << "Joint genotype all source VCF files, and generate a project VCF file" << endl
         << "on standard out. The source files must be in GVCF format." << endl
         << "Options:" << endl
         << "  --help, -h      print this help message" << endl
         << "  --bed FILE, -b FILE  path to three-column BED file" << endl
         << "  --residuals, -r      generate detailed residuals output file" << endl;
}

// Expected usage:
//    glnexus_cli_aio proj.vcf.gz [vcf files]
//
int main(int argc, char *argv[]) {
    GLnexus::Status s;
    spdlog::set_level(spdlog::level::info);
    spdlog::set_pattern("[%t] %+");
    console->info() << "glnexus_cli_aio " << GIT_REVISION;

    if (argc <= 2) {
        help(argv[0]);
        return 1;
    }

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"bed", required_argument, 0, 'b'},
        {"residuals", no_argument, 0, 'r'},
        {0, 0, 0, 0}
    };

    int c;
    bool residuals = false;
    string bedfilename;
    int nr_threads = std::thread::hardware_concurrency();

    optind = 1; // force optind past command positional argument
    while (-1 != (c = getopt_long(argc, argv, "hb:",
                                  long_options, nullptr))) {
        switch (c) {
            case 'b':
                bedfilename = string(optarg);
                if (bedfilename.size() == 0) {
                    cerr <<  "invalid BED filename" << endl;
                    return 1;
                }
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

    S(all_in_one(vcf_files, bedfilename, nr_threads, residuals));
    return 0;
}

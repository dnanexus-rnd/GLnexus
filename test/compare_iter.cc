#include <iostream>
#include <fstream>
#include <map>
#include "spdlog/spdlog.h"
#include "BCFKeyValueData.h"
#include "BCFSerialize.h"
#include "cli_utils.h"
#include "catch.hpp"

using namespace std;
using namespace GLnexus;

using IterResults = map<string, shared_ptr<vector<shared_ptr<bcf1_t>>>>;
using T = BCFKeyValueData;

// This module provides utility functions for comparing BCF iterators

// Seed for pseudo random number generator
#define PSEUDO_RAND_SEED (91451103L)

// maximal number of entries that we keep in memory
#define MAX_NUM_ENTRIES  (2000000)

// generate a random number in the range [0 .. n-1]
int gen_rand_number(int n){
    static bool firstTime = true;

    // initialization
    if (firstTime) {
        firstTime = false;
        srand(PSEUDO_RAND_SEED);
    }

    int i = rand() % n;
    return i;
}

double gen_rand_double(int n) {
    return double(gen_rand_number(n)) / double(n);
}

static void print_bcf_record(bcf1_t *x) {
    cout << "rid=" << x->rid << endl;
    cout << "pos=" << x->pos << endl;
    cout << "rlen=" << x->rlen << endl;
    cout << "qual=" << x->qual << endl;
    cout << "n_info=" << x->n_info << endl;
    cout << "n_allele=" << x->n_allele << endl;
    cout << "n_sample=" << x->n_sample << endl;
    cout << "shared.l=" << x->shared.l << endl;
    cout << "indiv.l=" << x->indiv.l << endl;
}

static int calc_tot_num_entries(shared_ptr<IterResults> results) {
    int nEntries = 0;
    for(IterResults::iterator iter = results->begin();
        iter != results->end(); ++iter) {
        shared_ptr<vector<shared_ptr<bcf1_t>>> v = iter->second;
        nEntries += v->size();
    }
    return nEntries;
}

// Read an entire iterator, and put into memory.
static void read_entire_iter(RangeBCFIterator *iter, shared_ptr<IterResults> results) {
    Status s;

    string dataset;
    shared_ptr<const bcf_hdr_t> hdr;
    do {
        vector<shared_ptr<bcf1_t>> v;
        v.clear();
        s = iter->next(dataset, hdr, v);
        for (auto rec : v) {
            if (results->find(dataset) == results->end()) {
                auto w = make_shared<vector<shared_ptr<bcf1_t>>>();
                (*results)[dataset] = w;
            }
            results->at(dataset)->push_back(rec);
        }
    } while (s.ok());
}

// Return true if the results are equivalent, false otherwise
static int compare_results(shared_ptr<IterResults> resultsBase,
                            shared_ptr<IterResults> resultsSoph) {
    for(IterResults::iterator iter = resultsBase->begin();
        iter != resultsBase->end();
        ++iter) {
        string dataset = iter->first;
        shared_ptr<vector<shared_ptr<bcf1_t>>> v = iter->second;
        shared_ptr<vector<shared_ptr<bcf1_t>>> w = resultsSoph->at(dataset);

        assert(v->size() == w->size());
        int len = v->size();
        //cout << dataset << " len=" << len << endl;
        for (int i=0; i < len; i++) {
            if (bcf_shallow_compare((*v)[i].get(), (*w)[i].get()) == 0) {
                cout << "Error comparing two records"  << endl;
                print_bcf_record((*v)[i].get());
                print_bcf_record((*w)[i].get());
                return 0;
            }
        }
    }

    return 1;
}

// Read all the records in the range, and return as a vector of BCF records.
// This is not a scalable function, because the return vector could be very
// large. For debugging/testing use only.
int compare_query(T &data, MetadataCache &cache,
                  const std::string& sampleset, const range& rng) {
    Status s;
    vector<unique_ptr<RangeBCFIterator>> iterators;
    shared_ptr<const set<string>> samples, datasets;

    // simple iterator
    s = data.sampleset_range_base(cache, sampleset, rng, 0,
                                  samples, datasets, iterators);
    assert(s.ok());

    auto resultsBase = make_shared<IterResults>();
    for (int i=0; i < iterators.size(); i++) {
        read_entire_iter(iterators[i].get(), resultsBase);

        if (calc_tot_num_entries(resultsBase) > MAX_NUM_ENTRIES) {
            // We are overrunning the memory limits, abort
            cout << "compare_query " << rng.str() << " ran over memory limit, aborting" << endl;
            return -1;
        }
    }
    iterators.clear();

    // sophisticated iterator
    auto resultsSoph = make_shared<IterResults>();
    s = data.sampleset_range(cache, sampleset, rng, 0,
                             samples, datasets, iterators);
    assert(s.ok());
    for (int i=0; i < iterators.size(); i++) {
        read_entire_iter(iterators[i].get(), resultsSoph);

        if (calc_tot_num_entries(resultsBase) > MAX_NUM_ENTRIES) {
            // We are overrunning the memory limits, abort
            cout << "compare_query " << rng.str() << " ran over memory limit, aborting" << endl;
            return -1;
        }
    }

    cout << "compare_query " << rng.str()
         << " num_entries=" << calc_tot_num_entries(resultsSoph) << endl;

    // compare all the elements between the two results
    return compare_results(resultsBase, resultsSoph);
}

static auto console = spdlog::stderr_logger_mt("GLnexus");

TEST_CASE("iter_compare") {
    Status s;

    // setup database directory
    string dbdir = "/tmp/iter_compare";
    REQUIRE(system(("rm -rf " + dbdir).c_str()) == 0);
    REQUIRE(system(("mkdir -p " + dbdir).c_str()) == 0);
    string dbpath = dbdir + "/DB";

    string basedir = "test/data/cli";
    string exemplar_gvcf = basedir + "/" + "F1.gvcf.gz";
    vector<pair<string,size_t>> contigs;
    s = cli::utils::db_init(console, dbpath, exemplar_gvcf, contigs);
    REQUIRE(s.ok());
    REQUIRE(contigs.size() >= 1);

    vector<string> gvcfs;
    for (auto fname : {"F1.gvcf.gz", "F2.gvcf.gz"}) {
         gvcfs.push_back(basedir + "/" + fname);
    }
    vector<range> ranges;
    s = cli::utils::db_bulk_load(console, 8, gvcfs, dbpath, ranges, contigs);
    REQUIRE(s.ok());
    REQUIRE(contigs.size() >= 1);

    unique_ptr<KeyValue::DB> db;
    unique_ptr<BCFKeyValueData> data;
    string sampleset;
    s = RocksKeyValue::Open(dbpath, db, cli::utils::GLnexus_prefix_spec(),
                                     RocksKeyValue::OpenMode::READ_ONLY);
    REQUIRE(s.ok());
    s = BCFKeyValueData::Open(db.get(), data);
    REQUIRE(s.ok());

    unique_ptr<MetadataCache> metadata;
    s = MetadataCache::Start(*data, metadata);
    REQUIRE(s.ok());

    s = data->all_samples_sampleset(sampleset);
    REQUIRE(s.ok());
    console->info() << "using sample set " << sampleset;

    // get samples and datasets
    shared_ptr<const set<string>> samples, datasets;
    s = metadata->sampleset_datasets(sampleset, samples, datasets);
    REQUIRE(s.ok());

    int nChroms = (int) (min((size_t)22, contigs.size()));
    int nIter = 50;
    int maxRangeLen = 1000000;
    int minLen = 10; // ensure that the the range is of some minimal size

    for (int i = 0; i < nIter; i++) {
        int rid = gen_rand_number(nChroms);
        int lenChrom = (int)contigs[rid].second;
        assert(lenChrom > minLen);

        // bound the range to be no larger than the chromosome
        int rangeLen = min(maxRangeLen, (lenChrom/10));
        assert(rangeLen > minLen);

        int beg = gen_rand_number(lenChrom - rangeLen);
        int rlen = gen_rand_number(rangeLen - minLen);
        range rng(rid, beg, beg + minLen + rlen);

        int rc = compare_query(*data, *metadata, sampleset, rng);
        REQUIRE(rc != 0);
        // 1: success
        // 0: error
        // -1: query used too much memory, aborted
    }

    cout << "Passed " << nIter << " iterator comparison tests" << endl;
}

#include <iostream>
#include <fstream>
#include <map>
#include "compare_queries.h"

namespace GLnexus {
namespace compare_queries {

using namespace std;
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
    cerr << "rid=" << x->rid << endl;
    cerr << "pos=" << x->pos << endl;
    cerr << "rlen=" << x->rlen << endl;
    cerr << "qual=" << x->qual << endl;
    cerr << "n_info=" << x->n_info << endl;
    cerr << "n_allele=" << x->n_allele << endl;
    cerr << "n_sample=" << x->n_sample << endl;
    cerr << "shared.l=" << x->shared.l << endl;
    cerr << "indiv.l=" << x->indiv.l << endl;
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
                cerr << "Error comparing two records"  << endl;
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
            cerr << "compare_query " << rng.str() << " ran over memory limit, aborting" << endl;
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
            cerr << "compare_query " << rng.str() << " ran over memory limit, aborting" << endl;
            return -1;
        }
    }

//    cerr << "compare_query " << rng.str()
//         << " num_entries=" << calc_tot_num_entries(resultsSoph) << endl;

    // compare all the elements between the two results
    return compare_results(resultsBase, resultsSoph);
}


Status compare_n_queries(int n_iter,
                         BCFKeyValueData &data,
                         MetadataCache &metadata,
                         const std::string& sampleset) {
    // get the contigs
    Status s;
    std::vector<std::pair<std::string,size_t> > contigs;
    S(data.contigs(contigs));

    int n_chroms = (int) (min((size_t)22, contigs.size()));
    int max_range_len = 1000000;
    int min_len = 10; // ensure that the the range is of some minimal size

    for (int i = 0; i < n_iter; i++) {
        int rid = gen_rand_number(n_chroms);
        int len_chrom = (int)contigs[rid].second;
        assert(len_chrom > min_len);

        // bound the range to be no larger than the chromosome
        int range_len = min(max_range_len, (len_chrom/10));
        assert(range_len > min_len);

        int beg = gen_rand_number(len_chrom - range_len);
        int rlen = gen_rand_number(range_len - min_len);
        range rng(rid, beg, beg + min_len + rlen);

        int rc = compare_query(data, metadata, sampleset, rng);

        // 1: success
        // 0: error
        // -1: query used too much memory, aborted
        if (rc == 0) {
            return Status::Failure("query comparison failed");
        }
    }

    return Status::OK();
}

}}

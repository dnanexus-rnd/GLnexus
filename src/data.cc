#include <assert.h>
#include "data.h"
#include "fcmm.hpp"
#include "khash.h"

using namespace std;

namespace GLnexus {

// hash<string> using the string hash function from htslib
class KStringHash {
public:
    size_t operator()(string const& s) const  {
        return (size_t) kh_str_hash_func(s.c_str());
    }
};
using StringCache = fcmm::Fcmm<string,string,hash<string>,KStringHash>;
using StringSetCache = fcmm::Fcmm<string,shared_ptr<const set<string>>,hash<string>,KStringHash>;
// this is not a hard limit but the FCMM performance degrades if it's too low
const size_t CACHE_SIZE = 4096;

struct MetadataCache::body {
    Metadata* inner;
    vector<pair<string,size_t> > contigs;
    unique_ptr<StringSetCache> sampleset_samples_cache;
    unique_ptr<StringCache> sample_dataset_cache;
    unique_ptr<StringSetCache> sampleset_datasets_cache;
};

MetadataCache::MetadataCache() = default;
MetadataCache::~MetadataCache() = default;

Status MetadataCache::Start(Metadata& inner, unique_ptr<MetadataCache>& ptr) {
    ptr.reset(new MetadataCache);
    ptr->body_.reset(new MetadataCache::body);
    ptr->body_->inner = &inner;
    ptr->body_->sampleset_samples_cache = make_unique<StringSetCache>(CACHE_SIZE);
    ptr->body_->sample_dataset_cache = make_unique<StringCache>(16 * CACHE_SIZE);
    ptr->body_->sampleset_datasets_cache = make_unique<StringSetCache>(CACHE_SIZE);
    return ptr->body_->inner->contigs(ptr->body_->contigs);
}

Status MetadataCache::contigs(vector<pair<string,size_t> >& ans) const {
    ans = body_->contigs;
    return Status::OK();
}

Status MetadataCache::sampleset_samples(const string& sampleset, shared_ptr<const set<string> >& ans) const {
    auto cached = body_->sampleset_samples_cache->end();
    if ((cached = body_->sampleset_samples_cache->find(sampleset))
            != body_->sampleset_samples_cache->end()) {
        ans = cached->second;
        assert(ans);
        return Status::OK();
    }

    Status s;
    S(body_->inner->sampleset_samples(sampleset, ans));
    body_->sampleset_samples_cache->insert(make_pair(sampleset,ans));
    return Status::OK();
}

Status MetadataCache::sample_dataset(const string& sample, string& ans) const {
    auto cached = body_->sample_dataset_cache->end();
    if ((cached = body_->sample_dataset_cache->find(sample))
            != body_->sample_dataset_cache->end()) {
        ans = cached->second;
        return Status::OK();
    }

    Status s;
    S(body_->inner->sample_dataset(sample, ans));
    body_->sample_dataset_cache->insert(make_pair(sample,ans));
    return Status::OK();
}

Status MetadataCache::all_samples_sampleset(string& ans) {
    // not safe to cache this as it's not immutable
    return body_->inner->all_samples_sampleset(ans);
}

Status MetadataCache::sample_count(size_t& ans) const {
    return body_->inner->sample_count(ans);
}

const vector<pair<string,size_t> >& MetadataCache::contigs() const {
    return body_->contigs;
}

Status MetadataCache::sampleset_datasets(const string& sampleset,
                                         shared_ptr<const set<string>>& samples,
                                         shared_ptr<const set<string>>& datasets_out) const {
    Status s;
    S(sampleset_samples(sampleset, samples));

    auto cached = body_->sampleset_datasets_cache->end();
    if ((cached = body_->sampleset_datasets_cache->find(sampleset))
            != body_->sampleset_datasets_cache->end()) {
        datasets_out = cached->second;
        assert(datasets_out);
        return Status::OK();
    }

    auto datasets = make_shared<set<string>>();
    for (const auto& it : *samples) {
        string dataset;
        S(sample_dataset(it, dataset));
        datasets->insert(dataset);
    }
    datasets_out = datasets;
    body_->sampleset_datasets_cache->insert(make_pair(sampleset,datasets_out));
    return Status::OK();
}

Status BCFData::dataset_range_and_header(const string& dataset, const range& pos, bcf_predicate predicate,
                                         shared_ptr<const bcf_hdr_t>* hdr,
                                         vector<shared_ptr<bcf1_t>>* records) {
    Status s;
    S(dataset_header(dataset, hdr));
    return dataset_range(dataset, hdr->get(), pos, predicate, records);
}

// default sampleset_range implementation:

// Return one iterator per 100kbp of the requested range. Each iterator simply
// calls dataset_range_in_header on each dataset in the 100kb block.

class DefaultRangeBCFIteratorImpl : public RangeBCFIterator {
    BCFData& data_;
    range range_;
    bool first_range_;
    shared_ptr<const set<string>> datasets_;
    set<string>::const_iterator it_;
    bcf_predicate predicate_;

public:
    DefaultRangeBCFIteratorImpl(BCFData& data, range range, bool first_range, bcf_predicate predicate,
                                shared_ptr<const set<string>>& datasets)
        : data_(data), range_(range), first_range_(first_range),
          datasets_(datasets), it_(datasets->begin()), predicate_(predicate) {}

    Status next(string& dataset, shared_ptr<const bcf_hdr_t>& hdr,
                vector<shared_ptr<bcf1_t>>& records) override {
        if (it_ == datasets_->end()) {
            return Status::NotFound();
        }
        dataset = *it_++;

        vector<shared_ptr<bcf1_t>> all_records;
        Status s = data_.dataset_range_and_header(dataset, range_, predicate_, &hdr, &all_records);
        if (s.bad()) {
             if (s == StatusCode::NOT_FOUND) {
                // censor this error so caller doesn't think this is the normal
                // end of iteration
                return Status::Failure("RangeBCFIterator::next(): unexpected NotFound error", dataset);
            }
            return s;
        }

        // subtlety: to ensure the iterators together return each record
        // exactly once, we may need to skip records that don't begin within
        // range_ (unless this is the first or lowest-range iterator)
        assert(std::is_sorted(all_records.begin(), all_records.end(),
                              [] (shared_ptr<bcf1_t>& p1, shared_ptr<bcf1_t>& p2) {
                                return range(p1) < range(p2);
                              }));
        int skip=0;
        if (!first_range_) {
            for (; skip<all_records.size(); skip++) {
                range rng(all_records[skip]);
                assert(range_.overlaps(rng));
                if (rng.beg >= range_.beg) {
                    break;
                }
            }
        }
        records.clear();
        records.insert(records.begin(), all_records.begin()+skip, all_records.end());

        return Status::OK();
    }
};

/* A naive implementation that splits the range into fixed sized
 * sub-ranges. Inside a sub-range, it iterates over the datasets, and
 * reads their records with point lookups.
 *
 * This is an inefficient DB access pattern, because it ignores the fact
 * that the data is ordered lexicographically by bucket, and then dataset. In other words,
 * data for one bucket for all datasets lies adjacently on disk.
 */
Status BCFData::sampleset_range(const MetadataCache& metadata, const string& sampleset,
                                const range& pos, bcf_predicate predicate,
                                shared_ptr<const set<string>>& samples,
                                shared_ptr<const set<string>>& datasets,
                                vector<unique_ptr<RangeBCFIterator>>& iterators) {
    const int RANGE_STEP = 100000;
    Status s;
    S(metadata.sampleset_datasets(sampleset, samples, datasets));

    iterators.clear();
    bool first = true;
    for (int beg = pos.beg; beg < pos.end; beg += RANGE_STEP) {
        range sub(pos.rid, beg, min(pos.end,beg+RANGE_STEP));
        iterators.push_back(make_unique<DefaultRangeBCFIteratorImpl>(*this, sub, first, predicate, datasets));
        first = false;
    }

    return Status::OK();
}


}

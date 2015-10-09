#include <assert.h>
#include "data.h"
#include "fcmm.hpp"
#include "khash.h"

using namespace std;

namespace GLnexus {

// std::hash<string> using the string hash function from htslib
class KStringHash {
public:
    std::size_t operator()(string const& s) const  {
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

Status BCFData::dataset_range_and_header(const string& dataset, const range& pos,
                                         shared_ptr<const bcf_hdr_t>& hdr,
                                         vector<shared_ptr<bcf1_t> >& records) {
    Status s;
    S(dataset_header(dataset, hdr));
    return dataset_range(dataset, hdr.get(), pos, records);
}

}

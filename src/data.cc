#include <assert.h>
#include "data.h"

using namespace std;

namespace GLnexus {

struct MetadataCache::body {
    Metadata* inner;
    vector<pair<string,size_t> > contigs;
};

MetadataCache::MetadataCache() = default;
MetadataCache::~MetadataCache() = default;

Status MetadataCache::Start(Metadata& inner, unique_ptr<MetadataCache>& ptr) {
    ptr.reset(new MetadataCache);
    ptr->body_.reset(new MetadataCache::body);
    ptr->body_->inner = &inner;
    return ptr->body_->inner->contigs(ptr->body_->contigs);
}

Status MetadataCache::contigs(vector<pair<string,size_t> >& ans) const {
    ans = body_->contigs;
    return Status::OK();
}

Status MetadataCache::sampleset_samples(const string& sampleset, shared_ptr<const set<string> >& ans) const {
    // TODO cache
    return body_->inner->sampleset_samples(sampleset, ans);
}

Status MetadataCache::sample_dataset(const string& sample, string& ans) const {
    // TODO cache
    return body_->inner->sample_dataset(sample, ans);
}

const vector<pair<string,size_t> >& MetadataCache::contigs() const {
    return body_->contigs;
}

Status MetadataCache::sampleset_datasets(const string& sampleset,
                                         shared_ptr<const set<string>>& samples,
                                         shared_ptr<const set<string>>& datasets_out) const {
    // TODO cache
    Status s;
    auto datasets = make_shared<set<string>>();
    S(sampleset_samples(sampleset, samples));
    for (const auto& it : *samples) {
        string dataset;
        S(sample_dataset(it, dataset));
        datasets->insert(dataset);
    }
    datasets_out = datasets;
    return Status::OK();
}

Status BCFData::dataset_range_and_header(const string& dataset, const range& pos, shared_ptr<const bcf_hdr_t>& hdr, vector<shared_ptr<bcf1_t> >& records) {
    Status s;
    S(dataset_header(dataset, hdr));
    return dataset_range(dataset, hdr.get(), pos, records);
}

}

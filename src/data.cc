#include <assert.h>
#include "data.h"

using namespace std;

namespace GLnexus {

struct DataCache::body {
    Data* data;
    vector<pair<string,size_t> > contigs;
};

Status DataCache::Init(Data* data) {
    body_ = new body;
    body_->data = data;
    return data->contigs(body_->contigs);
}

DataCache::~DataCache() {
    if (body_) {
        delete body_;
    }
}

Status DataCache::contigs(vector<pair<string,size_t> >& ans) const {
    return body_->data->contigs(ans);
}

Status DataCache::sampleset_samples(const string& sampleset, shared_ptr<const set<string> >& ans) const {
    // TODO cache
    return body_->data->sampleset_samples(sampleset, ans);
}

Status DataCache::sample_dataset(const string& sample, string& ans) const {
    // TODO cache
    return body_->data->sample_dataset(sample, ans);
}

Status DataCache::dataset_bcf_header(const string& dataset, shared_ptr<const bcf_hdr_t>& hdr) const {
    // TODO cache
    return body_->data->dataset_bcf_header(dataset, hdr);
}

Status DataCache::dataset_bcf(const string& dataset, const range& pos, shared_ptr<const bcf_hdr_t>& hdr, vector<shared_ptr<bcf1_t> >& records) const {
    return body_->data->dataset_bcf(dataset, pos, hdr, records);
}

const vector<pair<string,size_t> >& DataCache::contigs() const {
    return body_->contigs;
}

Status DataCache::sampleset_datasets(const string& sampleset, shared_ptr<const set<string>>& ans) const {
    // TODO cache
    shared_ptr<const set<string> > samples;
    auto datasets = make_shared<set<string>>();
    Status s;
    S(sampleset_samples(sampleset, samples));
    for (const auto& it : *samples) {
        string dataset;
        S(sample_dataset(it, dataset));
        datasets->insert(dataset);
    }
    ans = datasets;
    return Status::OK();
}

}

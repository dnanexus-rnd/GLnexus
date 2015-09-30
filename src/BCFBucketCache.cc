#include <thread>
#include <mutex>
#include "BCFSerialize.h"
#include "BCFBucketCache.h"

using namespace std;

namespace GLnexus {

// pImpl idiom
struct BCFBucketCache_body {
    KeyValue::DB* db;
};

BCFBucketCache::BCFBucketCache() = default;
BCFBucketCache::~BCFBucketCache() = default;


Status BCFBucketCache::Open(KeyValue::DB* db,
                            std::unique_ptr<BCFBucketCache>& ans) {
    assert(db != nullptr);

    ans.reset(new BCFBucketCache());
    ans->body_.reset(new BCFBucketCache_body);
    ans->body_->db = db;

    return Status::OK();
}

// Get a shared read-only pointer to a bucket.
Status BCFBucketCache::get_bucket(const string& key,
                                  const string& dataset,
                                  const range& r,
                                  shared_ptr<vector<shared_ptr<bcf1_t> > >& ans) {
    // Retrieve the pertinent DB entries
    ans = make_shared<vector<shared_ptr<bcf1_t> > >();
    Status s;
    KeyValue::CollectionHandle coll;
    S(body_->db->collection("bcf",coll));

    string data;
    s = body_->db->get(coll, key, data);
    if (s == StatusCode::NOT_FOUND) {
        return Status::NotFound();
    }

    unique_ptr<BCFReader> reader;
    S(BCFReader::Open(data.c_str(), data.size(), reader));

    shared_ptr<bcf1_t> vt;
    while ((s = reader->read(vt)).ok()) {
        assert(vt);
        if (bcf_unpack(vt.get(), BCF_UN_ALL) != 0) {
            return Status::IOError("BCFKeyValueData::dataset_bcf bcf_unpack",
                                   dataset + "@" + r.str());
        }
        ans->push_back(vt);
        vt.reset(); // important! otherwise reader overwrites the stored copy.
    }
    return Status::OK();
}

} // namespace GLnexus

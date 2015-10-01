#include <thread>
#include <mutex>
#include "BCFSerialize.h"
#include "BCFBucketCache.h"
#include "rocksdb/cache.h"

using namespace std;

namespace GLnexus {

// pImpl idiom
struct BCFBucketCache_body {
    KeyValue::DB* db;
    KeyValue::CollectionHandle coll;
    CacheMode mode;
    shared_ptr<rocksdb::Cache> cache;
};

// Type of values stored in a buckets
using BktT = shared_ptr<vector<shared_ptr<bcf1_t> > >;

BCFBucketCache::BCFBucketCache() = default;
BCFBucketCache::~BCFBucketCache() = default;


Status BCFBucketCache::Open(KeyValue::DB* db,
                            CacheMode mode,
                            int nBCFcapacity,
                            std::unique_ptr<BCFBucketCache>& ans) {
    assert(db != nullptr);

    if (mode == CacheMode::DISABLE) {
        // Don't waste memory if there if the cache is disabled.
        nBCFcapacity = 0;
    }

    ans.reset(new BCFBucketCache());
    ans->body_.reset(new BCFBucketCache_body);
    ans->body_->db = db;
    ans->body_->mode = mode;
    ans->body_->cache = rocksdb::NewLRUCache(nBCFcapacity);

    // calculate once the BCF collection handle
    Status s;
    S(ans->body_->db->collection("bcf", ans->body_->coll));

    return Status::OK();
}

// Get a shared read-only pointer to a bucket.
static Status get_bucket_from_db(BCFBucketCache_body *body_,
                                 const string& key,
                                 const string& dataset,
                                 const range& r,
                                 StatsRangeQuery &accu,
                                 shared_ptr<vector<shared_ptr<bcf1_t> > >& ans) {
    // Retrieve the pertinent DB entries
    string data;
    Status s = body_->db->get(body_->coll, key, data);
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
    accu.nBCFRecordsReadFromDB += ans->size();
    return Status::OK();
}

// Release the memory held by a bucket
static void bucket_mem_free(const rocksdb::Slice& key, void* val) {
//    cout << "---  bucket_mem_free ---" << endl;
    BktT *bucketPtr = reinterpret_cast<BktT*>(val);
    bucketPtr->reset();
    free(bucketPtr);
}

Status BCFBucketCache::get_bucket(const string& key,
                                  const string& dataset,
                                  const range& r,
                                  StatsRangeQuery &accu,
                                  shared_ptr<vector<shared_ptr<bcf1_t> > >& ans) {
    if (body_->mode == CacheMode::DISABLE ) {
        // no real caching
        if (ans == nullptr) {
            ans =  make_shared<vector<shared_ptr<bcf1_t> > >();
        }
        return get_bucket_from_db(body_.get(), key, dataset, r, accu, ans);
    }

    // Check if the bucket is in memory. If so, hand a shared
    // pointer to the caller.
    rocksdb::Cache::Handle *hndl = body_->cache->Lookup(key);
    if (hndl != nullptr) {
        void *val = body_->cache->Value(hndl);
        BktT bucket = *reinterpret_cast<BktT*>(val);
        ans = bucket;
        body_->cache->Release(hndl);
        return Status::OK();
    }

    // The bucket is not in memory. Read it from the DB and insert into
    // the cache.
    BktT* bucketPtr = new BktT;
    bucketPtr->reset(new vector<shared_ptr<bcf1_t> > );
    Status s = get_bucket_from_db(body_.get(), key, dataset, r, accu, *bucketPtr);
    int nRecordsInBucket = (*bucketPtr)->size();
    if (s.bad()) {
        return s;
    }
    ans = *bucketPtr;

//    cout << "Insert into cache capacity=" << nRecordsInBucket
//         << "/" << body_->cache->GetUsage() << endl;
    body_->cache->Insert(key, reinterpret_cast<void*>(bucketPtr),
                         nRecordsInBucket, bucket_mem_free);

    return Status::OK();
}

} // namespace GLnexus

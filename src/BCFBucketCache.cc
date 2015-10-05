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
    int64_t capacityRAM;
    shared_ptr<rocksdb::Cache> cache;
};

// Type of values stored in a buckets
using BktT = vector<shared_ptr<bcf1_t> >;

BCFBucketCache::BCFBucketCache() = default;

BCFBucketCache::~BCFBucketCache() {
    if (body_->cache != nullptr) {
        body_->cache.reset();
    }
}


Status BCFBucketCache::Open(KeyValue::DB* db,
                            int capacityRAM,
                            std::unique_ptr<BCFBucketCache>& ans) {
    assert(db != nullptr);

    ans.reset(new BCFBucketCache());
    ans->body_.reset(new BCFBucketCache_body);
    ans->body_->db = db;
    ans->body_->capacityRAM = capacityRAM;
    if (capacityRAM > 0) {
        ans->body_->cache = rocksdb::NewLRUCache(capacityRAM);
    }

    // calculate once the BCF collection handle
    Status s;
    S(ans->body_->db->collection("bcf", ans->body_->coll));

    return Status::OK();
}

// Get a shared read-only pointer to a bucket.
// [memCost] is an estimate for how much memory is used by this bucket.
static Status get_bucket_from_db(BCFBucketCache_body *body_,
                                 const string& key,
                                 StatsRangeQuery &accu,
                                 int &memCost,
                                 vector<shared_ptr<bcf1_t> > &ans) {
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
            return Status::IOError("BCFKeyValueData::dataset_bcf bcf_unpack", key);
        }
        ans.push_back(vt);
        vt.reset(); // important! otherwise reader overwrites the stored copy.
    }
    memCost = data.size() * 2;
    accu.nBCFRecordsReadFromDB += ans.size();
    return Status::OK();
}

// Release the memory held by a bucket
static void delete_cached_bucket(const rocksdb::Slice& key, void* val) {
    if (val == NULL)
        return;
    BktT *bucketPtr = static_cast<BktT*>(val);
    delete bucketPtr;
}

Status BCFBucketCache::get_bucket(const string& key,
                                  StatsRangeQuery &accu,
                                  void *&bucket_hndl,
                                  vector<shared_ptr<bcf1_t> > & ans) {
    bucket_hndl = NULL;
    int memCost = 0;
    if (body_->capacityRAM == 0) {
        // no real caching
        return get_bucket_from_db(body_.get(), key, accu, memCost, ans);
    }

    // Check if the bucket is in memory. If so, hand a shared
    // pointer to the caller.
    rocksdb::Slice sliceKey(key);
    rocksdb::Cache::Handle *hndl = body_->cache->Lookup(sliceKey);
    if (hndl == nullptr) {
        // The bucket is not in memory. Read it from the DB and insert into
        // the cache.
        BktT* bucketPtr = new BktT;
        Status s = get_bucket_from_db(body_.get(), key, accu, memCost, *bucketPtr);
        if (s.bad()) {
            delete bucketPtr;
            return s;
        }

        //cout << "Insert into cache " << memCost << "/" << body_->cache->GetUsage() << endl;
        body_->cache->Insert(sliceKey, static_cast<void*>(bucketPtr),
                             memCost, &delete_cached_bucket);

        // This should now succeed
        hndl = body_->cache->Lookup(sliceKey);
    }

    void *val = body_->cache->Value(hndl);
    ans = *static_cast<BktT*>(val);
    bucket_hndl = static_cast<void*>(hndl);
    return Status::OK();
}

void BCFBucketCache::release_bucket(void *bucket_hndl) {
    rocksdb::Cache::Handle *hndl = static_cast<rocksdb::Cache::Handle*>(bucket_hndl);
    body_->cache->Release(hndl);
}

} // namespace GLnexus

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
    size_t capacityRAM;
    shared_ptr<rocksdb::Cache> cache;
};

// Type of values stored in a buckets
using BktT = vector<shared_ptr<bcf1_t> >;

BCFBucketCache::BCFBucketCache() = default;
BCFBucketCache::~BCFBucketCache() = default;


Status BCFBucketCache::Open(KeyValue::DB* db,
                            size_t capacityRAM,
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
                                 BktT **ans) {
    // Retrieve the pertinent DB entries
    string data;
    Status s = body_->db->get(body_->coll, key, data);
    if (s == StatusCode::NOT_FOUND) {
        return Status::NotFound();
    }

    unique_ptr<BCFReader> reader;
    S(BCFReader::Open(data.c_str(), data.size(), reader));

    *ans = new BktT;
    shared_ptr<bcf1_t> vt;
    while ((s = reader->read(vt)).ok()) {
        assert(vt);
        if (bcf_unpack(vt.get(), BCF_UN_ALL) != 0) {
            delete *ans;
            return Status::IOError("BCFKeyValueData::dataset_bcf bcf_unpack", key);
        }
        (*ans)->push_back(vt);
        vt.reset(); // important! otherwise reader overwrites the stored copy.
    }
    memCost = data.size() * 2;
    accu.nBCFRecordsReadFromDB += (*ans)->size();
    return Status::OK();
}

// Release the memory held by a bucket (RocksDB cache deleter function)
static void delete_cached_bucket(const rocksdb::Slice& key, void* val) {
    assert(val != nullptr);
    BktT *bucketPtr = static_cast<BktT*>(val);
    delete bucketPtr;
}

Status BCFBucketCache::get_bucket(const string& key,
                                  StatsRangeQuery &accu,
                                  shared_ptr<BktT>& ans) {
    Status s;
    int memCost = 0;
    if (body_->capacityRAM == 0) {
        // no real caching
        BktT *bucketPtr;
        s = get_bucket_from_db(body_.get(), key, accu, memCost, &bucketPtr);
        if (s.ok()) {
            ans.reset(bucketPtr);
        }
        return s;
    }

    // Check if the bucket is in memory. If so, hand a shared
    // pointer to the caller.
    rocksdb::Slice sliceKey(key);
    rocksdb::Cache::Handle *hndl = body_->cache->Lookup(sliceKey);
    if (hndl == nullptr) {
        // The bucket is not in memory. Read it from the DB and insert into
        // the cache.
        BktT *bucketPtr;
        S(get_bucket_from_db(body_.get(), key, accu, memCost, &bucketPtr));

        //cout << "Insert into cache " << memCost << "/" << body_->cache->GetUsage() << endl;
        hndl = body_->cache->Insert(sliceKey, bucketPtr, memCost, &delete_cached_bucket);
        assert(hndl != nullptr);
    }

    // Return a shared_ptr to the bucket with an unusual custom deleter: it
    // releases the cache handle instead of doing anything to directly free
    // the BktT object itself. The RocksDB cache handles the reference
    // counting on the BktT object instead of the shared_ptr.
    void *val = body_->cache->Value(hndl);
    ans = shared_ptr<BktT>(static_cast<BktT*>(val),
                           [this, hndl] (BktT*) { this->body_->cache->Release(hndl); });
    return Status::OK();
}

} // namespace GLnexus

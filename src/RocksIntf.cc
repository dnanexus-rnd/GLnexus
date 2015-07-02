// Implement a KeyValue interface to a RocksDB on-disk database.
//
// Interface that we need to support.
#include "KeyValue.h"

// RocksDB header files that are needed.
#include "rocksdb/db.h"
#include "rocksdb/slice.h"
#include "rocksdb/options.h"

std::string kDBPath = "/tmp/rocksdb_simple_example";

namespace GLnexus {
namespace RocksIntf {

    // Convert a rocksdb Status into a GLnexus Status structure
    static Status convertStatus(const rocksdb::Status &s)
    {
        rocksdb::Status::Code code = s.code();
        switch (code) {
            // good case
        case kOk: return Status::OK();

            // error cases
        case kNotFound: return Status::NotFound();
        case kCorruption: return Status::Failure("corruption");
        case kNotSupported: return Status::NotImplemented();
        case kInvalidArgument: return Status::Invalid();
        case kIOError: return Status::IOError();
        case kMergeInProgress: return Status::Failure("merge in progress");
        case kIncomplete: return Status::Failure("incomplete");
        case kShutdownInProgress: return Status::Failure("shutdown in progress");
        case kTimedOut: return Status::Failure("timed out");
        case kAborted: return Status::Failure("aborted");

            // catch all for other cases that currently unlisted.
        default: return Status::Failure("other reason");
        }
    }


    class Iterator : public KeyValue::Iterator {
    private:
        rocksdb::Iterator* iter_;
        friend class Reader;

    public:
        Iterator(rocksdb::DB* db, const rocksdb::ColumnFamilyHandle& coll)
        {
            // FIXME: type safety problem here.
            // coll should be:
            //   rocksdb::ColumnFamilyHandle*
            // but is:
            //   const rocksdb::ColumnFamilyHandle&
            //
            rocksdb::ReadOptions options;  // default values
            iter_ = db->NewIterator(options, coll);
            iter_.SeekToFirst();
        }

        Status next(std::string& key, std::string& value) override {
            // It seems kind of silly to call Valid twice here. I am
            // not sure how to avoid doing this, while maintaining safety.
            //
            if (!iter_.Valid())
                return Status::NotFound();
            iter_.Next();
            if (!iter_.Valid())
                return Status::NotFound();
            key = iter_.key();
            value = iter_.value();
            return Status::OK();
        }

        Status Seek(std::string& key) {
            if (!iter_.Valid())
                return Status::NotFound();
            iter_.Seek(key);
            if (!iter_.Valid())
                return Status::NotFound();
            return Status::OK();
        }
    };


    class Reader : public KeyValue::Reader<rocksdb::ColumnFamilyHandle,Iterator> {
    private:
        rocksdb::DB* db_;

    public:
        Status get(const rocksdb::ColumnFamilyHandle* coll,
                   const std::string& key,
                   std::string& value) const override {
            const ReadOptions r_options; // what should this be set to?
            rocksdb::Status s = db_->Get(r_options, coll, key, value);
            return convertStatus(s);
        }

        Status iterator(const rocksdb::ColumnFamilyHandle& coll,
                        std::unique_ptr<Iterator>& it) const override {
            // FIXME: type safety problem here.
            // coll should be:
            //   rocksdb::ColumnFamilyHandle*
            // but is:
            //   const rocksdb::ColumnFamilyHandle&

            //
            // Make sure this column family actually exists
            rocksdb::ColumnFamilyMetaData *colMd;
            rocksdb::Status s = db_->GetColumnFamilyMetaData(coll, colMd);
            if (!s.ok())
                return convertStatus(s);

            it = std::make_unique<Iterator>(db_, coll);
            return Status::OK();
        }

        Status iterator(const rocksdb::ColumnFamilyHandle& coll,
                        const std::string& key,
                        std::unique_ptr<Iterator>& it) const override {
            //
            // Make sure this column family actually exists
            rocksdb::ColumnFamilyMetaData *colMd;
            rocksdb::Status s = db_->GetColumnFamilyMetaData(coll, colMd);
            if (!s.ok())
                return convertStatus(s);

            it = std::make_unique<Iterator>(db_, coll);
            return it->Seek(key);
        }
    };


    class WriteBatch : public KeyValue::WriteBatch<uint64_t> {
    public:
        Status put(const uint64_t& coll, const std::string& key, const std::string& value) override {
            assert(false);
            return Status::OK();
        }
    };


    class DB : public KeyValue::DB<uint64_t, RocksIntf::Reader, RocksIntf::Iterator, RocksIntf::WriteBatch> {
        // a collection is a partition of the key value space. Think, separate
        // databases. They all share the same log, so can be modified atomically.
        //
        // In RocksDB, this is called a column family.
        //
        // map a collection-name to an ID.
        // std::map<std::string,uint64_t> collections_;

        //std::vector<std::map<std::string,std::string>> data_;
    private:
        rocksdb::DB* db = null;

    public:
        DB() {
            // Default optimization choices, could use a deeper look.
            rocksdb::Options options;
            options.IncreaseParallelism();
            options.OptimizeLevelStyleCompaction();
            options.create_if_missing = true;

            // open DB
            printf("Opening database [...");
            rocksdb::Status s = DB::Open(options, kDBPath, &db);
            assert(s.ok());
            printf("]\n");
        }

        Status collection(const std::string& name, uint64_t& coll) const override {
            assert(false);
        }

        Status create_collection(const std::string& name) override {
            // FIXME: what do we do with the handle?
            ColumnFamilyOptions options; // defaults for now
            ColumnFamilyHandle *handle;
            rocksdb::Status s = rocksdb::CreateColumnFamily(options, name, handle);
            return Status::OK();
        }

        Status current(std::unique_ptr<RocksIntf::Reader>& reader) const override {
            assert(false);
        }

        Status begin_writes(std::unique_ptr<RocksIntf::WriteBatch>& writes) override {
            assert(false);
        }

        Status commit_writes(RocksIntf::WriteBatch* writes) override {
            assert(false);
        }

        void wipe() {
        }
    };
}}

// Implement a KeyValue interface to a RocksDB on-disk database.
//
#include <map>
#include "KeyValue.h"
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

        // No copying allowed
        Iterator(const Iterator&);
        void operator=(const Iterator&);

    public:
        // Note: the NewIterator call allocates heap-memory
        Iterator(rocksdb::DB* db, rocksdb::ColumnFamilyHandle* coll)
        {
            rocksdb::ReadOptions options;  // default values
            iter_ = db->NewIterator(options, coll);
            iter_.SeekToFirst();
        }

        // desctructor. We need to release the heap memory used by the
        // the iterator.
        ~Iterator() {
            delete iter_;
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

        Status seek(std::string& key) {
            if (!iter_.Valid())
                return Status::NotFound();
            iter_.Seek(key);
            if (!iter_.Valid())
                return Status::NotFound();
            return Status::OK();
        }
    };


    class Reader : public KeyValue::Reader<rocksdb::ColumnFamilyHandle*,Iterator> {
    private:
        rocksdb::DB* db_;

        // No copying allowed
        Reader(const Reader&);
        void operator=(const Reader&);
        
    public:
        Reader(rocksdb::DB *db) {
            db_ = db;
        }

        ~Reader() {}

        Status get(rocksdb::ColumnFamilyHandle* coll,
                   const std::string& key,
                   std::string& value) const override {
            const ReadOptions r_options; // what should this be set to?
            rocksdb::Status s = db_->Get(r_options, coll, key, value);
            return convertStatus(s);
        }

        Status iterator(rocksdb::ColumnFamilyHandle* coll,
                        std::unique_ptr<Iterator>& it) const override {
            // Make sure this column family actually exists
            rocksdb::ColumnFamilyMetaData *colMd;
            rocksdb::Status s = db_->GetColumnFamilyMetaData(coll, colMd);
            if (!s.ok())
                return convertStatus(s);

            it = std::make_unique<Iterator>(db_, coll);
            return Status::OK();
        }

        Status iterator(rocksdb::ColumnFamilyHandle* coll,
                        const std::string& key,
                        std::unique_ptr<Iterator>& it) const override {
            //
            // Make sure this column family actually exists
            rocksdb::ColumnFamilyMetaData *colMd;
            rocksdb::Status s = db_->GetColumnFamilyMetaData(coll, colMd);
            if (!s.ok())
                return convertStatus(s);

            it = std::make_unique<Iterator>(db_, coll);
            return it->seek(key);
        }
    };


    class WriteBatch : public KeyValue::WriteBatch<rocksdb::ColumnFamilyHandle*> {
    private:
        rocksdb::WriteBatch *wb_;

        // No copying allowed
        WriteBatch(const WriteBatch&);
        void operator=(const WriteBatch&);

    public:
        WriteBatch() {
            wb_ = new rocskdb::WriteBatch();
        }

        ~WriteBatch() {
            delete wb_;
        }

        Status put(rocksdb::ColumnFamilyHandle* coll,
                   const std::string& key,
                   const std::string& value) override {
            wb_->Put(coll, key, value);
        }
    };


    class DB : public KeyValue::DB<rocksdb::ColumnFamilyHandle*, Reader, Iterator, WriteBatch> {
    private:
        rocksdb::DB* db_ = null;
        std::map<const std::string& name, rocksdb::ColumnFamilyHandle*> coll2handle_;

        // No copying allowed
        DB(const DB&);
        void operator=(const DB&);

    public:
        DB() {
            // Default optimization choices, could use a deeper look.
            rocksdb::Options options;
            options.IncreaseParallelism();
            options.OptimizeLevelStyleCompaction();
            options.create_if_missing = true;

            // open DB
            rocksdb::Status s = DB::Open(options, kDBPath, &db_);
            assert(s.ok());

            coll2handle_ = new std::map<const std::string&, rocksdb::ColumnFamilyHandle*> ();
        }

        ~DB() {
            delete db;
            delete coll2handle_;
        }

        Status collection(const std::string& name, 
                          rocskdb::ColumnFamilyHandle*& coll) const override {
            auto p = coll2handle_.find(name);
            if (p != coll2handle_.end()) {
                coll = p->second;
                return Status::OK();
            }
            return Status::NotFound("column family does not exist", name);
        }

        Status create_collection(const std::string& name) override {
            if (coll2handle_.find(name) != coll2handle_.end()) {
                return Status::Exists("column family already exists", name);
            }

            // create new column family in rocksdb
            rocksdb::ColumnFamilyOptions options; // defaults for now
            rocksdb::ColumnFamilyHandle *handle;
            rocksdb::Status s = rocksdb::CreateColumnFamily(options, name, handle);
            if (!s.ok())
                return convertStatus(s);

            // success, add a mapping from the column family name to the handle.
            coll2handle_.put(name, handle);
            return Status::OK();
        }

        Status current(std::unique_ptr<Reader>& reader) const override {
            reader = std::make_unique<Reader>(db_);
            return Status::OK();
        }

        Status begin_writes(std::unique_ptr<WriteBatch>& writes) override {
            writes = std::make_unique<WriteBatch>();
            return Status::OK();
        }

        Status commit_writes(WriteBatch* updates) override {
            rocksdb::WriteOptions options;
            options.sync = true;
            rocksdb::Status s = db_->Write(options, updates);
            return convertStatus(s);
        }

        void wipe() {
            // Not sure what the semantics are supposed to be here. 
            rocksdb::Options options;
            db_->DestroyDB(kDBPath, options);
        }
    };
}}

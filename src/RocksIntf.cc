// Implement a KeyValue interface to a RocksDB on-disk database.
//
#include <cstddef>
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream> // for ostringstream
#include <string>
#include "KeyValue.h"
#include "RocksIntf.h"
#include "rocksdb/db.h"
#include "rocksdb/slice.h"
#include "rocksdb/options.h"
#include "rocksdb/write_batch.h"

std::string kDBPathBase = "/tmp/rocksdb_dir";
static int NUM_LIMIT = 1024 * 1024 * 1024;

namespace GLnexus {
namespace RocksIntf {

    // generate a random number in the range [0 .. n-1]
    static int genRandNumber(int n)
    {
        static bool firstTime = true;

        // initialization
        if (firstTime) {
            firstTime = false;
            srand (time(NULL));
        }

        int i = rand() % n;
        return i;
    }


    // Convert a rocksdb Status into a GLnexus Status structure
    static Status convertStatus(const rocksdb::Status &s)
    {
        rocksdb::Status::Code code = s.code();
        switch (code) {
            // good case
        case rocksdb::Status::kOk: return Status::OK();

            // error cases
        case rocksdb::Status::kNotFound:
            return Status::NotFound();
        case rocksdb::Status::kCorruption:
            return Status::Failure("corruption");
        case rocksdb::Status::kNotSupported:
            return Status::NotImplemented();
        case rocksdb::Status::kInvalidArgument:
            return Status::Invalid();
        case rocksdb::Status::kIOError:
            return Status::IOError();
        case rocksdb::Status::kMergeInProgress:
            return Status::Failure("merge in progress");
        case rocksdb::Status::kIncomplete:
            return Status::Failure("incomplete");
        case rocksdb::Status::kShutdownInProgress:
            return Status::Failure("shutdown in progress");
        case rocksdb::Status::kTimedOut:
            return Status::Failure("timed out");
        case rocksdb::Status::kAborted:
            return Status::Failure("aborted");

            // catch all for unlisted cases, all errors
        default: return Status::Failure("other reason");
        }
    }

    /*    void destroy()
    {
        rocksdb::Options options;
        rocksdb::Status s = rocksdb::DestroyDB(kDBPath, options);
        assert(s.ok());
        }*/


    class Iterator : public KeyValue::Iterator {
    private:
        rocksdb::Iterator* iter_;
        bool atStart_ = false;
        
        // No copying allowed
        Iterator(const Iterator&);
        void operator=(const Iterator&);

    public:
        // Note: the NewIterator call allocates heap-memory
        Iterator(rocksdb::DB* db, KeyValue::CollectionHandle _coll)
        {
            auto coll = reinterpret_cast<rocksdb::ColumnFamilyHandle*>(_coll);
            rocksdb::ReadOptions options;  // default values
            iter_ = db->NewIterator(options, coll);
            iter_->SeekToFirst();
            atStart_ = true;
        }

        // desctructor. We need to release the heap memory used by the
        // the iterator.
        ~Iterator() {
            delete iter_;
        }

        Status next(std::string& key, std::string& value) override {
            if (atStart_) {
                // We are at the first element
                atStart_ = false;
            }
            else {
                // We are beyong the first element 
                iter_->Next();
            }
            if (!iter_->Valid())
                return Status::NotFound();
            key = iter_->key().ToString();
            value = iter_->value().ToString();
            return Status::OK();
        }

        Status seek(const std::string& key) {
            if (!iter_->Valid())
                return Status::NotFound();
            iter_->Seek(key);
            if (!iter_->Valid())
                return Status::NotFound();
            return Status::OK();
        }
    };


    class Reader : public KeyValue::Reader {
    private:
        rocksdb::DB* db_ = NULL;

        // No copying allowed
        Reader(const Reader&);
        void operator=(const Reader&);

    public:
        Reader() {}

        Reader(rocksdb::DB *db) {
            db_ = db;
        }

        ~Reader() {}

        Status get(KeyValue::CollectionHandle _coll,
                   const std::string& key,
                   std::string& value) const override {
            auto coll = reinterpret_cast<rocksdb::ColumnFamilyHandle*>(_coll);
            const rocksdb::ReadOptions r_options; // what should this be set to?
            std::string v_tmp; // convert from pointer to reference, can we
            rocksdb::Status s = db_->Get(r_options, coll, key, &v_tmp);
            value = v_tmp;
            return convertStatus(s);
        }

        Status iterator(KeyValue::CollectionHandle _coll,
                        std::unique_ptr<KeyValue::Iterator>& it) const override {
            auto coll = reinterpret_cast<rocksdb::ColumnFamilyHandle*>(_coll);
            it = std::make_unique<Iterator>(db_, coll);
            return Status::OK();
        }

        Status iterator(KeyValue::CollectionHandle _coll,
                        const std::string& key,
                        std::unique_ptr<KeyValue::Iterator>& it) const override {
            auto coll = reinterpret_cast<rocksdb::ColumnFamilyHandle*>(_coll);
            auto ans = std::make_unique<Iterator>(db_, coll);
            Status s;
            S(ans->seek(key));
            it = std::move(ans);
            return Status::OK();
        }
    };

    class WriteBatch : public KeyValue::WriteBatch {
    private:
        rocksdb::WriteBatch *wb_;
        rocksdb::DB* db_;
        friend class DB;

        // No copying allowed
        WriteBatch(const WriteBatch&);
        void operator=(const WriteBatch&);

    public:
        WriteBatch(rocksdb::DB* db) : db_(db) {
            wb_ = new rocksdb::WriteBatch();
        }

        ~WriteBatch() {
            delete wb_;
            wb_ = NULL; // extra sanitation
            db_ = NULL;
        }

        Status put(KeyValue::CollectionHandle _coll,
                   const std::string& key,
                   const std::string& value) override {
            auto coll = reinterpret_cast<rocksdb::ColumnFamilyHandle*>(_coll);
            wb_->Put(coll, key, value);
            return Status::OK();
        }

        Status commit() override {
            rocksdb::WriteOptions options;
            options.sync = true;
            rocksdb::Status s = db_->Write(options, wb_);
            return convertStatus(s);
        }
    };


    class DB : public KeyValue::DB {
    private:
        rocksdb::DB* db_;
        std::map<const std::string, rocksdb::ColumnFamilyHandle*> coll2handle_;

        // No copying allowed
        DB(const DB&);
        void operator=(const DB&);

    public:
        DB(const std::vector<std::string>& collections) {
            // Default optimization choices, could use a deeper look.
            rocksdb::Options options;
            options.IncreaseParallelism();
            options.OptimizeLevelStyleCompaction();
            options.create_if_missing = true;

            // open DB
            std::ostringstream out;
            int rndNum = genRandNumber(NUM_LIMIT);
            out << kDBPathBase << "_" << std::setfill('0') << std::setw(10) << rndNum;
            std::string dbFullPath = out.str();
            rocksdb::Status s = rocksdb::DB::Open(options, dbFullPath, &db_);
            assert(s.ok());

            std::cout << "Construct RocksDB at "<< dbFullPath << "\n";
            std::cout.flush();

            // create initial collections
            for (const auto &colName : collections) {
                rocksdb::ColumnFamilyOptions options; // defaults for now
                rocksdb::ColumnFamilyHandle *handle;
                rocksdb::Status s = db_->CreateColumnFamily(options, colName, &handle);
                assert(s.ok());

                // success, add a mapping from the column family name to the handle.
                coll2handle_[colName] = handle;
            }
        }

        ~DB() override {
            printf("Delete RocksDB\n"); fflush(stdout);
            delete db_;
        }

        Status collection(const std::string& name,
                          KeyValue::CollectionHandle& coll) const override {
            auto p = coll2handle_.find(name);
            if (p != coll2handle_.end()) {
                coll = reinterpret_cast<KeyValue::CollectionHandle>(p->second);
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
            rocksdb::Status s = db_->CreateColumnFamily(options, name, &handle);
            if (!s.ok())
                return convertStatus(s);

            // success, add a mapping from the column family name to the handle.
            coll2handle_[name] = handle;
            return Status::OK();
        }

        Status current(std::unique_ptr<KeyValue::Reader>& reader) const override {
            reader = std::make_unique<RocksIntf::Reader>(db_);
            return Status::OK();
        }

        Status begin_writes(std::unique_ptr<KeyValue::WriteBatch>& writes) override {
            writes = std::make_unique<RocksIntf::WriteBatch>(db_);
            return Status::OK();
        }
    };

    Status make(const std::vector<std::string>& collections,
                std::unique_ptr<KeyValue::DB> &db,
                const std::string dbPath)
    {
        auto db_uq = std::make_unique<RocksIntf::DB>(collections);
        db = std::move(db_uq);
        if (db.get() != NULL)
            return Status::OK();
            return Status::Failure();
        return Status::OK();        
    }
}}

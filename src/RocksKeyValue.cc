// Implement a KeyValue interface to a RocksDB on-disk database.
//
#include <cstddef>
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream> // for ostringstream
#include <string>
#include "KeyValue.h"
#include "RocksKeyValue.h"
#include "rocksdb/db.h"
#include "rocksdb/slice.h"
#include "rocksdb/options.h"
#include "rocksdb/write_batch.h"

namespace GLnexus {
namespace RocksKeyValue {

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
        return Status::Invalid("RocksDB kInvalidArgument", s.ToString());
    case rocksdb::Status::kIOError:
        return Status::IOError("RocksDB kIOError", s.ToString());
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
    default: return Status::Failure("other reason", s.ToString());
    }
}

rocksdb::ColumnFamilyOptions GLnexusColumnFamilyOptions() {
    rocksdb::ColumnFamilyOptions opts;
    opts.OptimizeUniversalStyleCompaction();
    return opts;
}

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

    // Destructor. We need to release the heap memory used by the
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
        std::string v_tmp;
        rocksdb::Status s = db_->Get(r_options, coll, key, &v_tmp);
        value = std::move(v_tmp);
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

    DB(rocksdb::DB *db, std::map<const std::string, rocksdb::ColumnFamilyHandle*>& coll2handle)
        : db_(db), coll2handle_(std::move(coll2handle))
        {}

public:
    static Status Initialize(const std::string& dbPath,
                       std::unique_ptr<KeyValue::DB> &db) {
        rocksdb::Options options;
        options.IncreaseParallelism();
        options.OptimizeUniversalStyleCompaction();
        options.create_if_missing = true;
        options.error_if_exists = true;

        rocksdb::DB *rawdb = nullptr;
        rocksdb::Status s = rocksdb::DB::Open(options, dbPath, &rawdb);
        if (!s.ok()) {
            return convertStatus(s);
        }
        assert(rawdb != nullptr);

        std::map<const std::string, rocksdb::ColumnFamilyHandle*> coll2handle;
        db.reset(new DB(rawdb, coll2handle));
        if (!db) {
            delete rawdb;
            return Status::Failure();
        }
        return Status::OK();
    }

    static Status Open(const std::string& dbPath,
                       std::unique_ptr<KeyValue::DB> &db) {
        rocksdb::Options options;
        options.IncreaseParallelism();
        options.OptimizeLevelStyleCompaction();
        options.create_if_missing = false;

        // detect the database's column families
        std::vector<std::string> column_family_names;
        rocksdb::Status s = rocksdb::DB::ListColumnFamilies(options, dbPath, &column_family_names);
        if (!s.ok()) {
            return convertStatus(s);
        }
        std::vector<rocksdb::ColumnFamilyDescriptor> column_families;
        for (const auto& nm : column_family_names) {
            rocksdb::ColumnFamilyDescriptor cfd;
            cfd.name = nm;
            cfd.options = GLnexusColumnFamilyOptions();
            column_families.push_back(std::move(cfd));
        }

        // open the database (all column families)
        rocksdb::DB *rawdb = nullptr;
        std::vector<rocksdb::ColumnFamilyHandle*> column_family_handles;
        s = rocksdb::DB::Open(options, dbPath, column_families,
                              &column_family_handles, &rawdb);
        if (!s.ok()) {
            return convertStatus(s);
        }
        assert(rawdb != nullptr);

        // create the database object with coll2handle_ pre-filled
        std::map<const std::string, rocksdb::ColumnFamilyHandle*> coll2handle;
        for (size_t i = 0; i < column_families.size(); i++) {
            coll2handle[column_family_names[i]] = column_family_handles[i];
        }
        db.reset(new DB(rawdb, coll2handle));
        if (!db) {
            for (auto h : column_family_handles) {
                delete h;
            }
            delete rawdb;
            return Status::Failure();
        }
        return Status::OK();
    }

    ~DB() override {
        // free RocksDB column family handles
        for (const auto& p : coll2handle_) {
            delete p.second;
        }
        // delete database
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
        rocksdb::ColumnFamilyHandle *handle;
        rocksdb::Status s = db_->CreateColumnFamily(GLnexusColumnFamilyOptions(), name, &handle);
        if (!s.ok()) {
            return convertStatus(s);
        }
        assert(handle != nullptr);

        // success, add a mapping from the column family name to the handle.
        coll2handle_[name] = handle;
        return Status::OK();
    }

    Status current(std::unique_ptr<KeyValue::Reader>& reader) const override {
        reader = std::make_unique<RocksKeyValue::Reader>(db_);
        return Status::OK();
    }

    Status begin_writes(std::unique_ptr<KeyValue::WriteBatch>& writes) override {
        writes = std::make_unique<RocksKeyValue::WriteBatch>(db_);
        return Status::OK();
    }

    Status get(KeyValue::CollectionHandle _coll,
               const std::string& key,
               std::string& value) const override {
        auto coll = reinterpret_cast<rocksdb::ColumnFamilyHandle*>(_coll);
        const rocksdb::ReadOptions r_options; // what should this be set to?
        std::string v_tmp;
        rocksdb::Status s = db_->Get(r_options, coll, key, &v_tmp);
        value = std::move(v_tmp);
        return convertStatus(s);
    }

    Status put(KeyValue::CollectionHandle _coll,
               const std::string& key,
               const std::string& value) override {
        auto coll = reinterpret_cast<rocksdb::ColumnFamilyHandle*>(_coll);
        const rocksdb::WriteOptions options;
        rocksdb::Status s = db_->Put(options, coll, key, value);
        return convertStatus(s);
    }
};

Status Initialize(const std::string& dbPath, std::unique_ptr<KeyValue::DB>& db)
{
    return DB::Initialize(dbPath, db);
}

Status Open(const std::string& dbPath, std::unique_ptr<KeyValue::DB>& db)
{
    return DB::Open(dbPath, db);
}

Status destroy(const std::string dbPath)
{
    rocksdb::Options options;
    return convertStatus(rocksdb::DestroyDB(dbPath, options));
}



}}

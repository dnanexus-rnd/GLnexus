// Implement a KeyValue interface to a RocksDB on-disk database.
//
#include <cstddef>
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream> // for ostringstream
#include <string>
#include <thread>
#include <algorithm>
#include <unistd.h>
#include "KeyValue.h"
#include "RocksKeyValue.h"
#include "rocksdb/db.h"
#include "rocksdb/slice.h"
#include "rocksdb/options.h"
#include "rocksdb/write_batch.h"
#include "rocksdb/table.h"
#include "rocksdb/memtablerep.h"
#include "rocksdb/cache.h"
#include "rocksdb/slice_transform.h"

namespace GLnexus {
namespace RocksKeyValue {

static size_t totalRAM() {
    // http://nadeausoftware.com/articles/2012/09/c_c_tip_how_get_physical_memory_size_system
    static size_t memoized = 0;
    if (!memoized) {
        memoized = (size_t)sysconf( _SC_PHYS_PAGES ) * (size_t)sysconf( _SC_PAGESIZE );
        if (!memoized) {
            memoized = 4<<30;
        }
    }
    return memoized;
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
        return Status::NotImplemented("RocksDB kNotSupported", s.ToString());
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

// Reference for RocksDB tuning: https://github.com/facebook/rocksdb/wiki/RocksDB-Tuning-Guide
void ApplyColumnFamilyOptions(OpenMode mode, size_t prefix_length,
                              rocksdb::ColumnFamilyOptions& opts) {
    // universal compaction, 1GiB memtable budget
    opts.OptimizeUniversalStyleCompaction(1<<30);
    opts.num_levels = 4;
    opts.target_file_size_base = 10 * size_t(1<<30);
    opts.level0_file_num_compaction_trigger = 4;

    opts.compaction_options_universal.compression_size_percent = -1;
    opts.compaction_options_universal.max_size_amplification_percent = 300;
    opts.compaction_options_universal.size_ratio = 10;
    opts.compaction_options_universal.min_merge_width = 2;
    opts.compaction_options_universal.max_merge_width = 6;

    // 256 KiB blocks, with a large sharded cache
    rocksdb::BlockBasedTableOptions bbto;
    bbto.format_version = 2;
    bbto.block_size = 256 * 1024;
    bbto.block_cache = rocksdb::NewLRUCache(totalRAM() / 2, 8);

    // compress all files with Zstandard
    opts.compression_per_level.clear();
    opts.compression = rocksdb::kZSTD;
    opts.compression_opts.level = 2;
    opts.compression_opts.max_dict_bytes = bbto.block_size * 4;
    // TODO when released: https://github.com/facebook/rocksdb/commit/24ad4306001e0e9ec13ec29f11e50b4a7b070bf4
    // opts.compression_opts.zstd_max_train_bytes = opts.compression_opts.max_dict_bytes;

    if (prefix_length) {
        // prefix-based hash indexing for this column family
        opts.prefix_extractor.reset(rocksdb::NewFixedPrefixTransform(prefix_length));
        opts.memtable_factory.reset(rocksdb::NewHashSkipListRepFactory());
        bbto.index_type = rocksdb::BlockBasedTableOptions::kHashSearch;
    }

    opts.table_factory.reset(rocksdb::NewBlockBasedTableFactory(bbto));

    if (mode == OpenMode::BULK_LOAD) {
        // Use RocksDB's vector memtable implementation instead of the default
        // skiplist. Insertion to the vector memtable is much faster than the
        // skiplist (even though the latter supports concurrent writes), at the
        // cost of making interleaved reads prohibitively expensive, which is
        // a good tradeoff for bulk loading.
        opts.memtable_factory = std::make_shared<rocksdb::VectorRepFactory>();

        // Increase memtable size (and shrink block cache to compensate)
        opts.write_buffer_size = totalRAM() / 6;
        opts.max_write_buffer_number = 4;
        opts.min_write_buffer_number_to_merge = 1;
        bbto.block_cache = rocksdb::NewLRUCache(totalRAM() / 10, 6);

        // Never slowdown ingest since we'll wait for compaction to converge
        // at the end of the bulk load operation
        opts.level0_slowdown_writes_trigger = (1<<30);
        opts.level0_stop_writes_trigger = (1<<30);
        opts.soft_pending_compaction_bytes_limit = 0;
        opts.hard_pending_compaction_bytes_limit = 0;

        // Size amplification isn't really a thing during bulk loading because
        // nothing is getting deleted. The heuristic can also lead to merges
        // above max_merge_width, so disable it.
        opts.compaction_options_universal.max_size_amplification_percent = (1<<30);
    }
}

void ApplyDBOptions(OpenMode mode, rocksdb::Options& opts) {
    ApplyColumnFamilyOptions(mode, 0, static_cast<rocksdb::ColumnFamilyOptions&>(opts));

    if (mode == OpenMode::READ_ONLY) {
        // override db_log_dir when we're read-only, so that the RocksDB dir
        // stays untouched. https://github.com/facebook/rocksdb/issues/478
        opts.db_log_dir = "/tmp";
    }

    opts.max_open_files = -1;

    // configure parallelism
    opts.max_background_jobs = std::max(std::thread::hardware_concurrency(), 2U);
    opts.max_subcompactions = std::min(opts.max_background_jobs, 10);
    opts.env->SetBackgroundThreads(opts.max_background_jobs, rocksdb::Env::LOW);
    opts.env->SetBackgroundThreads(opts.max_background_jobs, rocksdb::Env::HIGH);

    opts.access_hint_on_compaction_start = rocksdb::Options::AccessHint::SEQUENTIAL;
    opts.compaction_readahead_size = 16 << 20;

    // legacy issue -- feature not supported by the memtable implemetations we select
    opts.allow_concurrent_memtable_write = false;

    if (mode == OpenMode::BULK_LOAD) {
        // don't throttle ingest, and free up disk space more frequently than the
        // 6h default interval
        opts.delayed_write_rate = 10ULL * (1<<30);
        opts.delete_obsolete_files_period_micros = 15ULL * 60 * 1000000;
    }
}

class Iterator : public KeyValue::Iterator {
private:
    std::unique_ptr<rocksdb::Iterator> iter_;
    rocksdb::Slice key_, value_;

    // No copying allowed
    Iterator(const Iterator&) = delete;
    void operator=(const Iterator&) = delete;

public:

    Iterator(std::unique_ptr<rocksdb::Iterator>&& iter) : iter_(move(iter)) {
        if (iter_->Valid()) {
            key_ = iter_->key();
            value_ = iter_->value();
        }
    }

    bool valid() const override {
        return iter_->Valid();
    }

    std::pair<const char*, size_t> key() const override {
        return std::make_pair(key_.data(), key_.size());
    }
    std::pair<const char*, size_t> value() const override {
        return std::make_pair(value_.data(), value_.size());
    }

    Status next() override {
        if (!iter_->status().ok()) {
            return convertStatus(iter_->status());
        }
        iter_->Next();
        if (!iter_->status().ok()) {
            return convertStatus(iter_->status());
        }
        if (iter_->Valid()) {
            key_ = iter_->key();
            value_ = iter_->value();
        }
        return Status::OK();
    }
};


class Reader : public KeyValue::Reader {
private:
    rocksdb::DB* db_ = nullptr;

    // No copying allowed
    Reader(const Reader&) = delete;
    void operator=(const Reader&) = delete;

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
                    const std::string& key,
                    std::unique_ptr<KeyValue::Iterator>& it) const override {
        auto coll = reinterpret_cast<rocksdb::ColumnFamilyHandle*>(_coll);
        rocksdb::ReadOptions options;  // default values
        std::unique_ptr<rocksdb::Iterator> rit(db_->NewIterator(options, coll));
        if (!rit) {
            return Status::Failure("rocksdb::DB::NewIterator()");
        }
        if (key.empty()) {
            rit->SeekToFirst();
        } else {
            rit->Seek(key);
        }
        if (!rit->status().ok()) {
            return convertStatus(rit->status());
        }
        it = std::make_unique<Iterator>(move(rit));
        return Status::OK();
    }
};

class WriteBatch : public KeyValue::WriteBatch {
private:
    rocksdb::WriteBatch *wb_;
    rocksdb::DB* db_;
    const rocksdb::WriteOptions& batch_write_options_;
    friend class DB;

    // No copying allowed
    WriteBatch(const WriteBatch&);
    void operator=(const WriteBatch&);

public:
    WriteBatch(rocksdb::DB* db, const rocksdb::WriteOptions& batch_write_options)
        : db_(db), batch_write_options_(batch_write_options) {
        wb_ = new rocksdb::WriteBatch();
    }

    ~WriteBatch() {
        delete wb_;
        wb_ = nullptr; // extra sanitation
        db_ = nullptr;
    }

    Status put(KeyValue::CollectionHandle _coll,
               const std::string& key,
               const std::string& value) override {
        auto coll = reinterpret_cast<rocksdb::ColumnFamilyHandle*>(_coll);
        wb_->Put(coll, key, value);
        return Status::OK();
    }

    Status commit() override {
        rocksdb::Status s = db_->Write(batch_write_options_, wb_);
        return convertStatus(s);
    }
};

class DB : public KeyValue::DB {
private:
    rocksdb::DB* db_;
    std::map<const std::string, rocksdb::ColumnFamilyHandle*> coll2handle_;
    OpenMode mode_;
    prefix_spec prefix_spec_;
    rocksdb::WriteOptions write_options_, batch_write_options_;

    // No copying allowed
    DB(const DB&);
    void operator=(const DB&);

    DB(rocksdb::DB *db, std::map<const std::string, rocksdb::ColumnFamilyHandle*>& coll2handle,
       OpenMode mode, prefix_spec* pfx)
        : db_(db), coll2handle_(std::move(coll2handle)),
          mode_(mode) {
            if (pfx) {
                prefix_spec_ = *pfx;
            }
            // prepare write options
            if (mode_ == OpenMode::BULK_LOAD) {
                write_options_.disableWAL = true;
                batch_write_options_.disableWAL = true;
            } else {
                batch_write_options_.sync = true;
            }
        }

public:
    static Status Initialize(const std::string& dbPath,
                             std::unique_ptr<KeyValue::DB> &db,
                             prefix_spec* pfx) {
        rocksdb::Options options;
        ApplyDBOptions(OpenMode::NORMAL, options);
        options.create_if_missing = true;
        options.error_if_exists = true;

        rocksdb::DB *rawdb = nullptr;
        rocksdb::Status s = rocksdb::DB::Open(options, dbPath, &rawdb);
        if (!s.ok()) {
            return convertStatus(s);
        }
        assert(rawdb != nullptr);

        std::map<const std::string, rocksdb::ColumnFamilyHandle*> coll2handle;
        db.reset(new DB(rawdb, coll2handle, OpenMode::NORMAL, pfx));
        if (!db) {
            delete rawdb;
            return Status::Failure();
        }
        return Status::OK();
    }

    static Status Open(const std::string& dbPath,
                       std::unique_ptr<KeyValue::DB> &db,
                       prefix_spec *pfx,
                       OpenMode mode) {
        // prepare options
        rocksdb::Options options;
        ApplyDBOptions(mode, options);
        options.create_if_missing = false;

        // detect the database's column families
        std::vector<std::string> column_family_names;
        rocksdb::Status s = rocksdb::DB::ListColumnFamilies(options, dbPath, &column_family_names);
        if (!s.ok()) {
            return convertStatus(s);
        }
        std::vector<rocksdb::ColumnFamilyDescriptor> column_families;
        for (const auto& nm : column_family_names) {
            rocksdb::ColumnFamilyOptions colopts;
            size_t effective_pfx = 0;
            if (pfx && nm == pfx->first) {
                effective_pfx = pfx->second;
            }
            ApplyColumnFamilyOptions(mode, effective_pfx, colopts);
            rocksdb::ColumnFamilyDescriptor cfd;
            cfd.name = nm;
            cfd.options = colopts;
            column_families.push_back(std::move(cfd));
        }

        // open the database (all column families)
        rocksdb::DB *rawdb = nullptr;
        std::vector<rocksdb::ColumnFamilyHandle*> column_family_handles;

        if (mode == OpenMode::READ_ONLY) {
            s = rocksdb::DB::OpenForReadOnly(options, dbPath, column_families,
                                             &column_family_handles, &rawdb);
        } else {
            s = rocksdb::DB::Open(options, dbPath, column_families,
                                  &column_family_handles, &rawdb);
        }
        if (!s.ok()) {
            return convertStatus(s);
        }
        assert(rawdb != nullptr);

        // create the database object with coll2handle_ pre-filled
        std::map<const std::string, rocksdb::ColumnFamilyHandle*> coll2handle;
        for (size_t i = 0; i < column_families.size(); i++) {
            coll2handle[column_family_names[i]] = column_family_handles[i];
        }
        db.reset(new DB(rawdb, coll2handle, mode, pfx));
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
        flush();
        if (mode_ == OpenMode::BULK_LOAD) {
            // Wait for compactions to converge. Specifically, wait until
            // there's no more than one background compaction running.
            // Argument: once that's the case, the number of sorted runs can't
            // be far above level0_file_num_compaction_trigger, or else a
            // second background compaction would start. And, while it'd
            // be nice to let the 'final' compaction finish, we don't want to
            // sit around waiting for an often-lengthy, single-threaded
            // operation.
            uint64_t num_running_compactions = 0;
            do {
                if (num_running_compactions) {
                    sleep(10);
                }
                num_running_compactions = 0;
                db_->GetIntProperty(rocksdb::DB::Properties::kNumRunningCompactions, &num_running_compactions);
            } while(num_running_compactions>1);
        }
        // Free column handles
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
        assert(name.size());
        if (coll2handle_.find(name) != coll2handle_.end()) {
            return Status::Exists("column family already exists", name);
        }

        // create new column family in rocksdb
        rocksdb::ColumnFamilyOptions colopts;
        size_t pfx = 0;
        if (name == prefix_spec_.first) {
            pfx = prefix_spec_.second;
        }
        ApplyColumnFamilyOptions(mode_, pfx, colopts);
        rocksdb::ColumnFamilyHandle *handle;
        rocksdb::Status s = db_->CreateColumnFamily(colopts, name, &handle);
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
        writes = std::make_unique<RocksKeyValue::WriteBatch>(db_, batch_write_options_);
        return Status::OK();
    }

    Status get(KeyValue::CollectionHandle _coll,
               const std::string& key,
               std::string& value) const override {
        auto coll = reinterpret_cast<rocksdb::ColumnFamilyHandle*>(_coll);
        static const rocksdb::ReadOptions r_options; // what should this be set to?
        std::string v_tmp;
        rocksdb::Status s = db_->Get(r_options, coll, key, &v_tmp);
        value = std::move(v_tmp);
        return convertStatus(s);
    }

    Status put(KeyValue::CollectionHandle _coll,
               const std::string& key,
               const std::string& value) override {
        auto coll = reinterpret_cast<rocksdb::ColumnFamilyHandle*>(_coll);
        rocksdb::Status s = db_->Put(write_options_, coll, key, value);
        return convertStatus(s);
    }

    Status flush() override {
        if (mode_ != OpenMode::READ_ONLY) {
            Status s;
            S(convertStatus(db_->SyncWAL()));
            for (const auto& p : coll2handle_) {
                S(convertStatus(db_->Flush(rocksdb::FlushOptions(), p.second)));
            }
        }
        return Status::OK();
    }
};

Status Initialize(const std::string& dbPath, std::unique_ptr<KeyValue::DB>& db, prefix_spec* pfx)
{
    return DB::Initialize(dbPath, db, pfx);
}

Status Open(const std::string& dbPath, std::unique_ptr<KeyValue::DB>& db,
            prefix_spec* pfx, OpenMode mode)
{
    return DB::Open(dbPath, db, pfx, mode);
}

Status destroy(const std::string dbPath)
{
    rocksdb::Options options;
    Status s = convertStatus(rocksdb::DestroyDB(dbPath, options));
    ignore_retval(system(("rm -rf " + dbPath).c_str()));
    return s;
}

}}

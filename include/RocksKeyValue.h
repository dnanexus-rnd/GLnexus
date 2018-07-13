#ifndef GLNEXUS_ROCKS_INTF_H
#define GLNEXUS_ROCKS_INTF_H

// Implement a KeyValue interface to a RocksDB on-disk database.
//
#include "KeyValue.h"
namespace GLnexus {
namespace RocksKeyValue {

/// Specification for collection key prefix. RocksDB can maintain a special
/// hash index to speed up lookups in collections in which the key space can
/// be partitioned into a number of buckets according to a key prefix of fixed
/// size, with the tradeoff that in-order iterators are restricted to one key
/// prefix. Activate this mode by providing a [prefix_spec] containing the
/// collection (only one is supported right now) and prefix size. Importantly,
/// the database must be initialized and always used with the same
/// [prefix_spec].
using prefix_spec = std::pair<std::string, size_t>;

/// Database open mode
enum class OpenMode {
    /// Online, read-write operations
    NORMAL,

    /// Read-only. Write operations will fail. Concurrent reads will be faster
    /// because less locking is needed.
    READ_ONLY,

    /// Offline bulk loading. Write operations will be faster but read
    /// performance may be reduced drastically. Transaction durability
    /// features are disabled; the database is liable to be corrupted if the
    /// DB object's destructor is not executed (e.g. crash or power failure)
    BULK_LOAD
};

struct config {
    prefix_spec *pfx = nullptr;
    OpenMode mode = OpenMode::NORMAL;
    size_t mem_budget = 0;
    size_t thread_budget = 0;
};

/// Initialize a new database. The parent directory must exist. Fails if the
/// path already exists.
Status Initialize(const std::string& dbpath, const config& cfg, std::unique_ptr<KeyValue::DB>& db);

/// Open an existing database.
Status Open(const std::string& dbPath, const config& cfg, std::unique_ptr<KeyValue::DB>& db);

// Delete an existing database.
Status destroy(const std::string dbPath);
}}

#endif

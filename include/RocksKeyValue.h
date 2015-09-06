#ifndef GLNEXUS_ROCKS_INTF_H
#define GLNEXUS_ROCKS_INTF_H

// Implement a KeyValue interface to a RocksDB on-disk database.
//
#include "KeyValue.h"

namespace GLnexus {
namespace RocksKeyValue {

/// Initialize a new database. The parent directory must exist. Fails if the
/// path already exists.
Status Initialize(const std::string& dbpath, std::unique_ptr<KeyValue::DB>& db);

/// Open an existing database.
Status Open(const std::string& dbPath, std::unique_ptr<KeyValue::DB>& db);

// Delete an existing database.
Status destroy(const std::string dbPath); 
}}

#endif

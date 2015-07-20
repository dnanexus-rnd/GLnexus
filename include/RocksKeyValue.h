#ifndef GLNEXUS_ROCKS_INTF_H
#define GLNEXUS_ROCKS_INTF_H

// Implement a KeyValue interface to a RocksDB on-disk database.
//
#include "KeyValue.h"

namespace GLnexus {
namespace RocksKeyValue {

Status Open(const std::vector<std::string>& collections,
            const std::string dbPath,
            std::unique_ptr<KeyValue::DB> &db);

void destroy(const std::string dbPath); 
}}

#endif

#ifndef GLNEXUS_ROCKS_INTF_H
#define GLNEXUS_ROCKS_INTF_H

// Implement a KeyValue interface to a RocksDB on-disk database.
//
#include "KeyValue.h"

namespace GLnexus {
    namespace RocksIntf {

        Status make(const std::vector<std::string>& collections,
                    std::unique_ptr<KeyValue::DB> &db,
                    const std::string dbPath);
    }
}

#endif

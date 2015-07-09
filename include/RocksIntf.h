#ifndef GLNEXUS_ROCKS_INTF_H
#define GLNEXUS_ROCKS_INTF_H

// Implement a KeyValue interface to a RocksDB on-disk database.
//
#include "KeyValue.h"

namespace GLnexus {
    KeyValue::DB* makeRocksIntf(const std::vector<std::string>& collections);
}

#endif

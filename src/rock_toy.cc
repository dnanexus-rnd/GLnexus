// Copyright (c) 2013, Facebook, Inc.  All rights reserved.
// This source code is licensed under the BSD-style license found in the
// LICENSE file in the root directory of this source tree. An additional grant
// of patent rights can be found in the PATENTS file in the same directory.
#include <cstdio>
#include <string>
#include <stdio.h>

#include "rocksdb/db.h"
#include "rocksdb/slice.h"
#include "rocksdb/options.h"

using namespace rocksdb;

std::string kDBPath = "/tmp/rocksdb_simple_example";

int main()
{
  DB* db;
  Options options;
  // Optimize RocksDB. This is the easiest way to get RocksDB to perform well
  options.IncreaseParallelism();
  options.OptimizeLevelStyleCompaction();
  // create the DB if it's not already present
  options.create_if_missing = true;

  // open DB
  printf("Starting database ...\n");
  Status s = DB::Open(options, kDBPath, &db);
  assert(s.ok());
  printf("done\n");

  // Put key-value
  printf("put\n");
  s = db->Put(WriteOptions(), "key1", "value");
  assert(s.ok());

  // get value
  printf("get\n");
  std::string value;
  s = db->Get(ReadOptions(), "key1", &value);
  assert(s.ok());
  assert(value == "value");

  // atomically apply a set of updates
  printf("atomic update\n");
  {
    WriteBatch batch;
    batch.Delete("key1");
    batch.Put("key2", value);
    s = db->Write(WriteOptions(), &batch);
  }

  printf("get key1\n");
  s = db->Get(ReadOptions(), "key1", &value);
  assert(s.IsNotFound());

  printf("get key2\n");
  db->Get(ReadOptions(), "key2", &value);
  assert(value == "value");

  delete db;

  return 0;
}

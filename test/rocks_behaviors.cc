// Basic behaviors of RocksDB

#include <cstdio>
#include <string>
#include <stdio.h>
#include <iostream>

#include "rocksdb/db.h"
#include "rocksdb/slice.h"
#include "rocksdb/options.h"
#include "catch.hpp"

using namespace rocksdb;

static std::string kDBPath = "/tmp/rocksdb_simple_example";
static int NUM_KEYS = 20;
static std::string COLL_NAME = "rocky2";

TEST_CASE("Basic operations on RocksDB") {
  // Remove the files under the DB path
  REQUIRE(system(("rm -rf " + kDBPath).c_str()) == 0);

  DB* db;
  Options options;

  // create the DB if it's not already present
  options.create_if_missing = true;

  // open DB
  Status s = DB::Open(options, kDBPath, &db);
  REQUIRE(s.ok());

  // create column 
  ColumnFamilyHandle* coll;
  s = db->CreateColumnFamily(options, COLL_NAME, &coll);
  REQUIRE(s.ok());
      
  // put operations
  for (int i=0; i < NUM_KEYS; i++) {
      std::string k = "key" + std::to_string(i);
      std::string v = "value" + std::to_string(i);
      const rocksdb::WriteOptions w_options;
      s = db->Put(w_options, coll, k, v);
      REQUIRE(s.ok());
  }

  // get operations
  for (int i=0; i < NUM_KEYS; i++) {
      std::string k = "key" + std::to_string(i);
      std::string v;
      const rocksdb::ReadOptions r_options;
      s = db->Get(r_options, coll, k, &v);
      REQUIRE(s.ok());
      REQUIRE(v == "value" + std::to_string(i));
  }

  // check that iteration works
  {
      rocksdb::ReadOptions options;  // default values
      rocksdb::Iterator* iter = db->NewIterator(options, coll);
      iter->SeekToFirst();

      std::string key = iter->key().ToString();
      std::string value = iter->value().ToString();
      //std::cout << "key=" << key << "   value=" << value << std::endl;
      for (int i=1; i < NUM_KEYS; i++) {
          REQUIRE(iter->Valid());

          iter->Next();
          std::string key = iter->key().ToString();
          std::string value = iter->value().ToString();
          //std::cout << "key=" << key << "   value=" << value << std::endl;
      }

      delete iter;
  }

  // update batch
  {
      WriteBatch batch;
      batch.Delete("key1");
      batch.Put("key2", "xyz");
      s = db->Write(WriteOptions(), &batch);
      REQUIRE(s.ok());
  }

  // Check that the batch was successful
  {
      std::string value;
      s = db->Get(ReadOptions(), "key1", &value);
      REQUIRE(s.IsNotFound());

      db->Get(ReadOptions(), "key2", &value);
      REQUIRE(value == "xyz");
  }
  
  delete coll;
  delete db;

  // destroy the database, so we can reuse the path
  DestroyDB(kDBPath, options);
  //REQUIRE(s.ok());

  // destroy is not sufficient, we also need to remove the files 
  // under the path
  REQUIRE(system(("rm -rf " + kDBPath).c_str()) == 0);
}

#include <cstddef>
#include <map>
#include <iostream>
#include <map>
#include "KeyValue.h"
#include "BCFKeyValueData.h"

#include "rocksdb/db.h"
#include "rocksdb/slice.h"
#include "rocksdb/options.h"
#include "catch.hpp"

using namespace std;
using namespace GLnexus;

// Implement a KeyValue interface to a RocksDB on-disk database.
//
std::string kDBPath = "/tmp/rocksdb_simple_example";

namespace RocksIntf {

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


    class Iterator : public KeyValue::Iterator {
    private:
        rocksdb::Iterator* iter_;

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
        }

        // desctructor. We need to release the heap memory used by the
        // the iterator.
        ~Iterator() {
            delete iter_;
        }

        Status next(std::string& key, std::string& value) override {
            // It seems kind of silly to call Valid twice here. I am
            // not sure how to avoid doing this, while maintaining safety.
            //
            if (!iter_->Valid())
                return Status::NotFound();
            iter_->Next();
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
            std::string* v_tmp; // convert from pointer to reference, can we
            rocksdb::Status s = db_->Get(r_options, coll, key, v_tmp);
            value = *v_tmp;
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

    class DB;
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
        rocksdb::DB* db_ = NULL;
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
            rocksdb::Status s = rocksdb::DB::Open(options, kDBPath, &db_);
            assert(s.ok());

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

        ~DB() {
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

        void wipe() {
            // Not sure what the semantics are supposed to be here.
            rocksdb::Options options;
            rocksdb::DestroyDB(kDBPath, options);
        }
    };
}

using T = BCFKeyValueData;

TEST_CASE("RocksDB construction on improperly initialized database") {
    vector<string> collections = {"header","bcf"};
    RocksIntf::DB db(collections);
    unique_ptr<T> data;
    REQUIRE(T::Open(&db, data) == StatusCode::INVALID);
}

TEST_CASE("RocksDB initialization") {
    RocksIntf::DB db({});
    auto contigs = {make_pair<string,uint64_t>("21", 1000000), make_pair<string,uint64_t>("22", 1000001)};
    REQUIRE(T::InitializeDB(&db, contigs).ok());
    unique_ptr<T> data;
    REQUIRE(T::Open(&db, data).ok());

    SECTION("contigs") {
        vector<pair<string,size_t>> contigs;
        Status s = data->contigs(contigs);
        REQUIRE(s.ok());
        REQUIRE(contigs.size() == 2);
        REQUIRE(contigs[0].first == "21");
        REQUIRE(contigs[0].second == 1000000);
        REQUIRE(contigs[1].first == "22");
        REQUIRE(contigs[1].second == 1000001);
    }

    SECTION("sampleset_samples") {
        KeyValue::CollectionHandle coll;
        string null(1, '\0');
        REQUIRE(db.collection("sampleset", coll).ok());
        REQUIRE(db.put(coll, "trio1", "").ok());
        REQUIRE(db.put(coll, "trio1" + null + "fa", "").ok());
        REQUIRE(db.put(coll, "trio1" + null + "mo", "").ok());
        REQUIRE(db.put(coll, "trio1" + null + "ch", "").ok());
        REQUIRE(db.put(coll, "trio2", "").ok());
        REQUIRE(db.put(coll, "trio2" + null + "fa2", "").ok());
        REQUIRE(db.put(coll, "trio2" + null + "mo2", "").ok());
        REQUIRE(db.put(coll, "trio2" + null + "ch2", "").ok());

        shared_ptr<const set<string>> samples;
        REQUIRE(data->sampleset_samples("trio1", samples).ok());
        REQUIRE(samples->size() == 3);
        REQUIRE(samples->find("fa") != samples->end());
        REQUIRE(samples->find("mo") != samples->end());
        REQUIRE(samples->find("ch") != samples->end());
        REQUIRE(samples->find("fa2") == samples->end());

        REQUIRE(data->sampleset_samples("trio2", samples).ok());
        REQUIRE(samples->size() == 3);
        REQUIRE(samples->find("fa2") != samples->end());
        REQUIRE(samples->find("mo2") != samples->end());
        REQUIRE(samples->find("ch2") != samples->end());
        REQUIRE(samples->find("fa") == samples->end());

        REQUIRE(data->sampleset_samples("bogus", samples) == StatusCode::NOT_FOUND);
    }

    SECTION("sample_dataset") {
        KeyValue::CollectionHandle coll;
        REQUIRE(db.collection("sample_dataset", coll).ok());
        REQUIRE(db.put(coll, "fa", "trio1").ok());
        REQUIRE(db.put(coll, "mo", "trio1").ok());
        REQUIRE(db.put(coll, "ch", "trio1").ok());
        REQUIRE(db.put(coll, "fa2", "trio2").ok());
        REQUIRE(db.put(coll, "mo2", "trio2").ok());
        REQUIRE(db.put(coll, "ch2", "trio2").ok());

        string dataset;
        REQUIRE(data->sample_dataset("fa", dataset).ok());
        REQUIRE(dataset == "trio1");
        REQUIRE(data->sample_dataset("ch", dataset).ok());
        REQUIRE(dataset == "trio1");
        REQUIRE(data->sample_dataset("mo2", dataset).ok());
        REQUIRE(dataset == "trio2");
        REQUIRE(data->sample_dataset("bogus", dataset) == StatusCode::NOT_FOUND);
    }
}

TEST_CASE("RocksDB::import_gvcf") {
    RocksIntf::DB db({});
    auto contigs = {make_pair<string,uint64_t>("21", 1000000)};
    REQUIRE(T::InitializeDB(&db, contigs).ok());
    unique_ptr<T> data;
    REQUIRE(T::Open(&db, data).ok());
    unique_ptr<DataCache> cache;
    REQUIRE(DataCache::Start(data.get(), cache).ok());

    SECTION("NA12878D_HiSeqX.21.10009462-10009469.gvcf") {
        Status s = data->import_gvcf(cache.get(), "NA12878D", "test/data/NA12878D_HiSeqX.21.10009462-10009469.gvcf");
        REQUIRE(s.ok());

        string dataset;
        REQUIRE(data->sample_dataset("NA12878", dataset).ok());
        REQUIRE(dataset == "NA12878D");
    }

    SECTION("incompatible contigs") {
        db.wipe();
        contigs = { make_pair<string,uint64_t>("21", 1000000), make_pair<string,uint64_t>("22", 1000000) };
        Status s = T::InitializeDB(&db, contigs);
        REQUIRE(s.ok());

        REQUIRE(DataCache::Start(data.get(), cache).ok());
        s = data->import_gvcf(cache.get(), "NA12878D", "test/data/NA12878D_HiSeqX.21.10009462-10009469.gvcf");
        REQUIRE(s == StatusCode::INVALID);
    }
}

TEST_CASE("RocksDB BCF retrieval") {
    RocksIntf::DB db({});
    auto contigs = {make_pair<string,uint64_t>("21", 1000000)};
    REQUIRE(T::InitializeDB(&db, contigs).ok());
    unique_ptr<T> data;
    REQUIRE(T::Open(&db, data).ok());
    unique_ptr<DataCache> cache;
    REQUIRE(DataCache::Start(data.get(), cache).ok());

    Status s = data->import_gvcf(cache.get(), "NA12878D", "test/data/NA12878D_HiSeqX.21.10009462-10009469.gvcf");
    REQUIRE(s.ok());

    SECTION("dataset_bcf_header") {
        shared_ptr<const bcf_hdr_t> hdr;
        s = data->dataset_bcf_header("NA12878D", hdr);
        REQUIRE(s.ok());

        vector<string> samples;
        unsigned n = bcf_hdr_nsamples(hdr.get());
        for (unsigned i = 0; i < n; i++) {
            samples.push_back(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, i)));
        }
        REQUIRE(samples.size() == 1);
        REQUIRE(samples[0] == "NA12878");
    }

    SECTION("dataset_bcf") {
        // get all records
        shared_ptr<const bcf_hdr_t> hdr;
        vector<shared_ptr<bcf1_t>> records;
        s = data->dataset_bcf("NA12878D", range(0, 0, 1000000000), hdr, records);
        REQUIRE(s.ok());

        REQUIRE(records.size() == 5);

        REQUIRE(records[0]->pos == 10009461);
        REQUIRE(records[1]->rlen == 2);
        REQUIRE(records[0]->n_allele == 2);
        REQUIRE(string(records[0]->d.allele[0]) == "T");
        REQUIRE(string(records[0]->d.allele[1]) == "<NON_REF>");
        REQUIRE(bcf_get_info(hdr.get(), records[0].get(), "END")->v1.i == 10009463); // nb END stays 1-based!

        REQUIRE(records[1]->pos == 10009463);
        REQUIRE(records[1]->rlen == 2);
        REQUIRE(records[1]->n_allele == 3);
        REQUIRE(string(records[1]->d.allele[0]) == "TA");
        REQUIRE(string(records[1]->d.allele[1]) == "T");
        REQUIRE(string(records[1]->d.allele[2]) == "<NON_REF>");

        REQUIRE(records[2]->rlen == 1);

        REQUIRE(records[4]->pos == 10009468);
        REQUIRE(records[4]->n_allele == 2);
        REQUIRE(string(records[4]->d.allele[0]) == "A");
        REQUIRE(string(records[4]->d.allele[1]) == "<NON_REF>");
        REQUIRE(bcf_get_info(hdr.get(), records[4].get(), "END")->v1.i == 10009471); // nb END stays 1-based!

        // subset of records
        s = data->dataset_bcf("NA12878D", range(0, 10009463, 10009466), hdr, records);
        REQUIRE(s.ok());

        REQUIRE(records.size() == 2);

        REQUIRE(records[0]->pos == 10009463);
        REQUIRE(records[0]->n_allele == 3);
        REQUIRE(string(records[0]->d.allele[0]) == "TA");
        REQUIRE(string(records[0]->d.allele[1]) == "T");
        REQUIRE(string(records[0]->d.allele[2]) == "<NON_REF>");

        REQUIRE(records[1]->pos == 10009465);
        REQUIRE(records[1]->n_allele == 2);
        REQUIRE(string(records[1]->d.allele[0]) == "A");
        REQUIRE(string(records[1]->d.allele[1]) == "<NON_REF>");

        // empty results
        s = data->dataset_bcf("NA12878D", range(0, 0, 1000), hdr, records);
        REQUIRE(s.ok());
        REQUIRE(records.size() == 0);

        s = data->dataset_bcf("NA12878D", range(1, 10009463, 10009466), hdr, records);
        REQUIRE(s.ok());
        REQUIRE(records.size() == 0);

        // bogus dataset
        s = data->dataset_bcf("bogus", range(1, 10009463, 10009466), hdr, records);
        REQUIRE(s == StatusCode::NOT_FOUND);
    }
}

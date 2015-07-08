#include <iostream>
#include <map>
#include "BCFKeyValueData.h"
#include "catch.hpp"
using namespace std;
using namespace GLnexus;

// Trivial in-memory KeyValue implementation used in unit tests.
namespace KeyValueMem {
    using namespace KeyValue;

    class Iterator : public KeyValue::Iterator {
        std::map<std::string,std::string> data_;
        std::map<std::string,std::string>::const_iterator it_;
        friend class Reader;

    public:
        Status next(std::string& key, std::string& value) override {
            if (it_ == data_.end()) return Status::NotFound();
            key = it_->first;
            value = it_->second;
            it_++;
            return Status::OK();
        }
    };

    class Reader : public KeyValue::Reader {
        std::vector<std::map<std::string,std::string>> data_;
        friend class DB;

    public:
        Status get(CollectionHandle _coll, const std::string& key, std::string& value) const override {
            auto coll = reinterpret_cast<uint64_t>(_coll);
            assert(coll < data_.size());
            const auto& m = data_[coll];
            auto p = m.find(key);
            if (p == m.end()) return Status::NotFound("key", key);
            value = p->second;
            return Status::OK();
        }

        Status iterator(CollectionHandle _coll, std::unique_ptr<KeyValue::Iterator>& it) const override {
            auto coll = reinterpret_cast<uint64_t>(_coll);
            assert(coll < data_.size());
            auto it2 = std::make_unique<Iterator>();
            it2->data_ = data_[coll];
            it2->it_ = it2->data_.begin();
            it.reset(it2.release());
            return Status::OK();
        }

        Status iterator(CollectionHandle _coll, const std::string& key, std::unique_ptr<KeyValue::Iterator>& it) const override {
            auto coll = reinterpret_cast<uint64_t>(_coll);
            assert(coll < data_.size());
            auto it2 = std::make_unique<Iterator>();
            it2->data_ = data_[coll];
            it2->it_ = it2->data_.lower_bound(key);
            it.reset(it2.release());
            return Status::OK();
        }
    };

    class DB;

    class WriteBatch : public KeyValue::WriteBatch {
        std::vector<std::map<std::string,std::string>> data_;
        DB* db_;
        friend class DB;

    public:
        Status put(CollectionHandle _coll, const std::string& key, const std::string& value) override {
            auto coll = reinterpret_cast<uint64_t>(_coll);
            assert(coll < data_.size());
            data_[coll][key] = value;
            return Status::OK();
        }
        Status commit() override;
    };

    class DB : public KeyValue::DB {
        std::map<std::string,uint64_t> collections_;
        std::vector<std::map<std::string,std::string>> data_;
        friend class WriteBatch;

    public:
        DB(const std::vector<std::string>& collections) {
            for (uint64_t i = 0; i < collections.size(); i++) {
                assert(collections_.find(collections[i]) == collections_.end());
                collections_[collections[i]] = i;
            }
            data_ = std::vector<std::map<std::string,std::string>>(collections_.size());
        }

        Status collection(const std::string& name, CollectionHandle& coll) const override {
            auto p = collections_.find(name);
            if (p != collections_.end()) {
                coll = reinterpret_cast<CollectionHandle>(p->second);
                return Status::OK();
            }
            return Status::NotFound("KeyValueMem::collection", name);
        }

        Status create_collection(const std::string& name) override {
            if (collections_.find(name) != collections_.end()) {
                return Status::Exists("key-value collection already exists", name);
            }
            collections_[name] = data_.size();
            data_.emplace_back();
            return Status::OK();
        }

        Status current(std::unique_ptr<KeyValue::Reader>& reader) const override {
            auto p = std::make_unique<KeyValueMem::Reader>();
            p->data_ = data_;
            reader = std::move(p);
            return Status::OK();
        }

        Status begin_writes(std::unique_ptr<KeyValue::WriteBatch>& writes) override {
            auto p = std::make_unique<KeyValueMem::WriteBatch>();
            p->db_ = this;
            p->data_ = std::vector<std::map<std::string,std::string>>(data_.size());
            writes = std::move(p);
            return Status::OK();
        }

        void wipe() {
            collections_.clear();
            data_.clear();
        }
    };

    Status WriteBatch::commit() {
        assert(data_.size() <= db_->data_.size());
        for (size_t i = 0; i < data_.size(); i++) {
            for (const auto& p : data_[i]) {
                db_->data_[i][p.first] = p.second;
            }
        }
        return Status::OK();
    }
}

using T = BCFKeyValueData;

TEST_CASE("BCFKeyValueData construction on improperly initialized database") {
    vector<string> collections = {"header","bcf"};
    KeyValueMem::DB db(collections);
    unique_ptr<T> data;
    REQUIRE(T::Open(&db, data) == StatusCode::INVALID);
}

TEST_CASE("BCFKeyValueData initialization") {
    KeyValueMem::DB db({});
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

TEST_CASE("BCFKeyValueData::import_gvcf") {
    KeyValueMem::DB db({});
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

TEST_CASE("BCFKeyValueData BCF retrieval") {
    KeyValueMem::DB db({});
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

#include <iostream>
#include <map>
#include "BCFKeyValueData.h"
#include "BCFSerialize.h"
#include "compare_iter.h"
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
        bool valid() const override {
            return it_ != data_.end();
        }

        pair<const char*,size_t> key() const override {
            return make_pair(it_->first.c_str(), it_->first.size());
        }

        pair<const char*,size_t> value() const override {
            return make_pair(it_->second.c_str(), it_->second.size());
        }

        Status next() override {
            if (it_ != data_.end()) {
                it_++;
            }
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

        Status iterator(CollectionHandle _coll, const std::string& key, std::unique_ptr<KeyValue::Iterator>& it) const override {
            auto coll = reinterpret_cast<uint64_t>(_coll);
            assert(coll < data_.size());
            auto it2 = std::make_unique<Iterator>();
            it2->data_ = data_[coll];
            if (key.empty()) {
                it2->it_ = it2->data_.begin();
            } else {
                it2->it_ = it2->data_.lower_bound(key);
            }
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

    SECTION("metadata") {
        vector<pair<string,size_t>> contigs;
        Status s = data->contigs(contigs);
        REQUIRE(s.ok());
        REQUIRE(contigs.size() == 2);
        REQUIRE(contigs[0].first == "21");
        REQUIRE(contigs[0].second == 1000000);
        REQUIRE(contigs[1].first == "22");
        REQUIRE(contigs[1].second == 1000001);

        KeyValue::CollectionHandle coll;
        REQUIRE(db.collection("sampleset", coll).ok());
        string version;
        REQUIRE(db.get(coll, "*", version).ok());
        REQUIRE(version == "0");
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
    vector<pair<string,uint64_t>> contigs = {make_pair<string,uint64_t>("21", 48129895)};
    REQUIRE(T::InitializeDB(&db, contigs).ok());
    unique_ptr<T> data;
    REQUIRE(T::Open(&db, data).ok());
    unique_ptr<MetadataCache> cache;
    REQUIRE(MetadataCache::Start(*data, cache).ok());
    set<string> samples_imported;
    size_t ct;

    SECTION("empty all_samples_sampleset") {
        REQUIRE(cache->sample_count(ct).ok());
        REQUIRE(ct == 0);

        shared_ptr<const set<string>> all;
        Status s = cache->sampleset_samples("*", all);
        REQUIRE(s == StatusCode::NOT_FOUND);

        string sampleset;
        s = cache->all_samples_sampleset(sampleset);
        REQUIRE(s.ok());

        s = cache->sampleset_samples(sampleset, all);
        REQUIRE(s.ok());
        REQUIRE(all->size() == 0);

        KeyValue::CollectionHandle coll;
        REQUIRE(db.collection("sampleset", coll).ok());
        string version;
        REQUIRE(db.get(coll, "*", version).ok());
        REQUIRE(version == "0");
    }

    SECTION("NA12878D_HiSeqX.21.10009462-10009469.gvcf") {
        Status s = data->import_gvcf(*cache, "NA12878D", "test/data/NA12878D_HiSeqX.21.10009462-10009469.gvcf", samples_imported);
        REQUIRE(s.ok());

        REQUIRE(cache->sample_count(ct).ok());
        REQUIRE(ct == 1);

        // check that internal version number of all-samples sampleset was
        // incremented
        KeyValue::CollectionHandle coll;
        REQUIRE(db.collection("sampleset", coll).ok());
        string version;
        REQUIRE(db.get(coll, "*", version).ok());
        REQUIRE(version == "1");

        // check the * sample set
        string dataset;
        REQUIRE(data->sample_dataset("NA12878", dataset).ok());
        REQUIRE(dataset == "NA12878D");

        shared_ptr<const set<string>> all;
        s = cache->sampleset_samples("*", all);
        REQUIRE(s == StatusCode::NOT_FOUND);

        string sampleset;
        s = cache->all_samples_sampleset(sampleset);
        REQUIRE(s.ok());

        s = cache->sampleset_samples(sampleset, all);
        REQUIRE(s.ok());
        REQUIRE(all->size() == 1);
        REQUIRE(*(all->begin()) == "NA12878");
    }

    SECTION("all_samples_sampleset") {
        Status s = data->import_gvcf(*cache, "1", "test/data/sampleset_range1.gvcf", samples_imported);
        REQUIRE(s.ok());

        REQUIRE(cache->sample_count(ct).ok());
        REQUIRE(ct == 1);

        KeyValue::CollectionHandle coll;
        REQUIRE(db.collection("sampleset", coll).ok());
        string version;
        REQUIRE(db.get(coll, "*", version).ok());
        REQUIRE(version == "1");

        string sampleset, sampleset2;
        s = cache->all_samples_sampleset(sampleset);
        REQUIRE(s.ok());

        shared_ptr<const set<string>> all;
        s = cache->sampleset_samples(sampleset, all);
        REQUIRE(s.ok());
        REQUIRE(all->size() == 1);

        // the sample set should now be memoized so long as no new samples are added
        s = cache->all_samples_sampleset(sampleset2);
        REQUIRE(s.ok());
        REQUIRE(sampleset == sampleset2);


        s = data->import_gvcf(*cache, "2", "test/data/sampleset_range2.gvcf", samples_imported);
        REQUIRE(s.ok());

        REQUIRE(cache->sample_count(ct).ok());
        REQUIRE(ct == 2);

        REQUIRE(db.get(coll, "*", version).ok());
        REQUIRE(version == "2");

        // now we should get a new sample set
        s = cache->all_samples_sampleset(sampleset);
        REQUIRE(s.ok());
        REQUIRE(sampleset != sampleset2);

        s = cache->sampleset_samples(sampleset, all);
        REQUIRE(s.ok());
        REQUIRE(all->size() == 2);

        s = cache->all_samples_sampleset(sampleset2);
        REQUIRE(s.ok());
        REQUIRE(sampleset == sampleset2);
    }

    SECTION("new_sampleset") {
        Status s = data->import_gvcf(*cache, "1", "test/data/sampleset_range1.gvcf", samples_imported);
        REQUIRE(s.ok());

        s = data->import_gvcf(*cache, "2", "test/data/sampleset_range2.gvcf", samples_imported);
        REQUIRE(s.ok());

        s = data->new_sampleset(*cache, "x", set<string>{"HX0001"});
        REQUIRE(s.ok());

        shared_ptr<const set<string>> samples;
        s = cache->sampleset_samples("x", samples);
        REQUIRE(s.ok());
        REQUIRE(*samples == set<string>{"HX0001"});

        s = data->new_sampleset(*cache, "y", set<string>{"HX0001","HX0002"});
        REQUIRE(s.ok());

        s = cache->sampleset_samples("y", samples);
        REQUIRE(s.ok());
        REQUIRE(*samples == set<string>({"HX0001","HX0002"}));

        // empty samples
        s = data->new_sampleset(*cache, "z", set<string>());
        REQUIRE(s == StatusCode::INVALID);

        // nonexistent sample
        s = data->new_sampleset(*cache, "z", set<string>{"HX0001","hX0002"});
        REQUIRE(s == StatusCode::NOT_FOUND);

        // duplicate sample set
        s = data->new_sampleset(*cache, "x", set<string>{"HX0001","HX0002"});
        REQUIRE(s == StatusCode::EXISTS);

        // TODO: test invalid sample set names
    }

    SECTION("incompatible contigs") {
        db.wipe();
        contigs = { make_pair<string,uint64_t>("21", 1000000), make_pair<string,uint64_t>("22", 1000000) };
        Status s = T::InitializeDB(&db, contigs);
        REQUIRE(s.ok());

        REQUIRE(MetadataCache::Start(*data, cache).ok());
        s = data->import_gvcf(*cache, "NA12878D", "test/data/NA12878D_HiSeqX.21.10009462-10009469.gvcf", samples_imported);
        REQUIRE(s == StatusCode::INVALID);

        // * sample set version number should NOT change
        KeyValue::CollectionHandle coll;
        REQUIRE(db.collection("sampleset", coll).ok());
        string version;
        REQUIRE(db.get(coll, "*", version).ok());
        REQUIRE(version == "0");
    }

    SECTION("detect bogus END field") {
        db.wipe();
        contigs = { make_pair<string,uint64_t>("21", 1000000) };
        Status s = T::InitializeDB(&db, contigs);
        REQUIRE(s.ok());

        REQUIRE(MetadataCache::Start(*data, cache).ok());
        s = data->import_gvcf(*cache, "NA12878D", "test/data/bogus_END.gvcf", samples_imported);
        REQUIRE(s == StatusCode::INVALID);

        // * sample set version number should NOT change
        KeyValue::CollectionHandle coll;
        REQUIRE(db.collection("sampleset", coll).ok());
        string version;
        REQUIRE(db.get(coll, "*", version).ok());
        REQUIRE(version == "0");
    }

    SECTION("reject duplicate data set") {
        Status s = data->import_gvcf(*cache, "x", "test/data/NA12878D_HiSeqX.21.10009462-10009469.gvcf", samples_imported);
        REQUIRE(s.ok());

        s = data->import_gvcf(*cache, "x", "test/data/NA12878D_HiSeqX.21.10009462-10009469.gvcf", samples_imported);
        REQUIRE(s == StatusCode::EXISTS);

        KeyValue::CollectionHandle coll;
        REQUIRE(db.collection("sampleset", coll).ok());
        string version;
        REQUIRE(db.get(coll, "*", version).ok());
        REQUIRE(version == "1");
    }

    SECTION("reject duplicate sample") {
        Status s = data->import_gvcf(*cache, "y", "test/data/NA12878D_HiSeqX.21.10009462-10009469.gvcf", samples_imported);
        REQUIRE(s.ok());

        s = data->import_gvcf(*cache, "z", "test/data/NA12878D_HiSeqX.21.10009462-10009469.gvcf", samples_imported);
        REQUIRE(s == StatusCode::EXISTS);

        KeyValue::CollectionHandle coll;
        REQUIRE(db.collection("sampleset", coll).ok());
        string version;
        REQUIRE(db.get(coll, "*", version).ok());
        REQUIRE(version == "1");
    }
}

TEST_CASE("BCFKeyValueData BCF retrieval") {
    KeyValueMem::DB db({});
    auto contigs = {make_pair<string,uint64_t>("21", 48129895)};
    REQUIRE(T::InitializeDB(&db, contigs).ok());
    unique_ptr<T> data;
    REQUIRE(T::Open(&db, data).ok());
    unique_ptr<MetadataCache> cache;
    REQUIRE(MetadataCache::Start(*data, cache).ok());
    set<string> samples_imported;

    Status s = data->import_gvcf(*cache, "NA12878D", "test/data/NA12878D_HiSeqX.21.10009462-10009469.gvcf", samples_imported);
    cout << s.str() << endl;
    REQUIRE(s.ok());

    SECTION("dataset_header") {
        shared_ptr<const bcf_hdr_t> hdr;
        s = data->dataset_header("NA12878D", hdr);
        REQUIRE(s.ok());

        vector<string> samples;
        unsigned n = bcf_hdr_nsamples(hdr.get());
        for (unsigned i = 0; i < n; i++) {
            samples.push_back(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, i)));
        }
        REQUIRE(samples.size() == 1);
        REQUIRE(samples[0] == "NA12878");
    }

    SECTION("dataset_range") {
        // get all records
        shared_ptr<const bcf_hdr_t> hdr;
        s = data->dataset_header("NA12878D", hdr);
        REQUIRE(s.ok());
        vector<shared_ptr<bcf1_t>> records;
        s = data->dataset_range("NA12878D", hdr.get(), range(0, 0, 1000000000), records);
        cout << s.str() << endl;
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
        s = data->dataset_range("NA12878D", hdr.get(), range(0, 10009463, 10009466), records);
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
        s = data->dataset_range("NA12878D", hdr.get(), range(0, 0, 1000), records);
        REQUIRE((records.size() == 0));

        s = data->dataset_range("NA12878D", hdr.get(), range(1, 10009463, 10009466), records);
        //REQUIRE(s == StatusCode::NOT_FOUND);
        REQUIRE((records.size() == 0));

        // bogus dataset
        s = data->dataset_range("bogus", hdr.get(), range(1, 10009463, 10009466), records);
        //REQUIRE(s == StatusCode::NOT_FOUND);
        REQUIRE((records.size() == 0));
    }
}

// Test subtle cases of record search. The terminology is, search for range X,
// where the records in the database are {A1, A2, ..} for dataset A.
//
TEST_CASE("BCFKeyValueData range overlap with a single dataset") {
    std::vector<int> intervals = {9, 11, 13, 10000};

    for (int ilen : intervals) {
        //cout << "interval_len=" << ilen << endl;
        KeyValueMem::DB db({});
        auto contigs = {make_pair<string,uint64_t>("21", 48129895)};

        // Buckets of size 9 break the ranges [1005 -- 1010] and [3004 -- 3006]
        // in two.
        REQUIRE(T::InitializeDB(&db, contigs, ilen).ok());
        unique_ptr<T> data;
        REQUIRE(T::Open(&db, data).ok());
        unique_ptr<MetadataCache> cache;
        REQUIRE(MetadataCache::Start(*data, cache).ok());
        set<string> samples_imported;

        Status s = data->import_gvcf(*cache, "synth_A", "test/data/synthetic_A.21.gvcf", samples_imported);
        if (s.bad()) {
            cout << s.str() << endl;
        }
        REQUIRE(s.ok());
        shared_ptr<const bcf_hdr_t> hdr;
        s = data->dataset_header("synth_A", hdr);
        REQUIRE(s.ok());

        vector<shared_ptr<bcf1_t>> records;
        /* Case 1
           |<- X ->)
           |<-A1->)       |<-A2->)
           Expected result: empty set
        */
        s = data->dataset_range("synth_A", hdr.get(), range(0, 1005, 1010), records);
        REQUIRE(s.ok());
        REQUIRE(records.size() == 0);

        /* Case 2
           |<- X ->)
           |<-A1->)  |<-A3->)
           |<-A2->)  |<-A4->)
           Expected result: {A1, A2, A3}
        */
        s = data->dataset_range("synth_A", hdr.get(), range(0, 2003, 2006), records);
        REQUIRE(s.ok());
//        for (auto r : records) {
//            cout << "r= " << r->rid << "," << r->pos << "," << r->rlen << "," << r->shared.s  << endl;
//        }
        REQUIRE(records.size() == 3);

        /*
          Case 3
          |<- X  ->)
          |<-        A1   ->)
          |<-A2 - >)  |<-A4 ->|
          |<-A3->)
          Expected result: {A1, A2, A3}
        */
        s = data->dataset_range("synth_A", hdr.get(), range(0, 3004, 3006), records);
        REQUIRE(s.ok());
        REQUIRE(records.size() == 3);

        std::shared_ptr<GLnexus::StatsRangeQuery> statsRq = data->getRangeStats();
        cout << statsRq->str() << endl;
    }
}

// --------------------------------------------------------------------
// This is a design for a test that will be useful with query primitives
// that work on multiple datasets.
//
// Test subtle cases of record search. The terminology is, search for range X,
// where the records in the database are {A1, A2, ..} for dataset A, and {B1, B2, ..}
// for dataset B.
//
// Case 1
//                        |<- X ->|
//      |<-A1->| |<-A2->|          |<-B1->| |<-B2->|
// Expected result: empty set
//
// Case 2
//           |<- X ->|
//       |<-A1->|  |<-A2->|
//              |<-B1->|  |<-B2->|
// Expected result: {A1, A2, B1}
//
// Case 3
//          |<-    X     ->|
//          |<-A1       ->|  |<- A2 ->|
// |<-B1->|   |< -B2   ->|
// Expected result: {A1, B2}
// --------------------------------------------------------------------

TEST_CASE("BCFData::sampleset_range") {
    // This tests the base implemention of sampleset_range in BCFData, which
    // returns iterators over 100kbp slices.

    KeyValueMem::DB db({});
    auto contigs = {make_pair<string,uint64_t>("21", 48129895)};
    REQUIRE(T::InitializeDB(&db, contigs).ok());
    unique_ptr<T> data;
    REQUIRE(T::Open(&db, data).ok());
    unique_ptr<MetadataCache> cache;
    REQUIRE(MetadataCache::Start(*data, cache).ok());
    set<string> samples_imported;

    Status s = data->import_gvcf(*cache, "1", "test/data/sampleset_range1.gvcf", samples_imported);
    REQUIRE(s.ok());
    s = data->import_gvcf(*cache, "2", "test/data/sampleset_range2.gvcf", samples_imported);
    REQUIRE(s.ok());

    size_t ct;
    REQUIRE(cache->sample_count(ct).ok());
    REQUIRE(ct == 2);

    // check * version number
    KeyValue::CollectionHandle coll;
    REQUIRE(db.collection("sampleset", coll).ok());
    string version;
    REQUIRE(db.get(coll, "*", version).ok());
    REQUIRE(version == "2");

    string sampleset;
    s = cache->all_samples_sampleset(sampleset);
    REQUIRE(s.ok());

    shared_ptr<const set<string>> samples, datasets;
    vector<unique_ptr<RangeBCFIterator>> iterators;

    // We'll slice various ranges and ensure we get the correct number of
    // iterators and that each iterator returns the appropriate data for each
    // dataset

    range rng(0, 100000, 200000);
    s = data->sampleset_range_base(*cache, sampleset, rng,
                              samples, datasets, iterators);
    REQUIRE(s.ok());
    REQUIRE(iterators.size() == 1);

    string dataset;
    shared_ptr<const bcf_hdr_t> hdr;
    vector<shared_ptr<bcf1_t>> records, all_records;
    s = iterators[0]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    #define check() s = data->dataset_range(dataset, hdr.get(), range(0,0,1000000), all_records); \
                    REQUIRE(s.ok()); \
                    REQUIRE(records.size() <= count_if(all_records.begin(), all_records.end(), [&](shared_ptr<bcf1_t>& r){return rng.overlaps(r.get());})); \
                    REQUIRE(all_of(records.begin(), records.end(), [&](shared_ptr<bcf1_t>& r){return rng.overlaps(r.get());}))
    REQUIRE(dataset == "1");
    REQUIRE(records.size() == 3);
    check();
    s = iterators[0]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(dataset == "2");
    REQUIRE(records.size() == 0);
    s = iterators[0]->next(dataset, hdr, records);
    REQUIRE(s == StatusCode::NOT_FOUND);

    rng = range(0, 100000, 200001);
    s = data->sampleset_range_base(*cache, sampleset, rng,
                              samples, datasets, iterators);
    REQUIRE(s.ok());
    REQUIRE(iterators.size() == 2);
    s = iterators[0]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(dataset == "1");
    REQUIRE(records.size() == 3);
    check();
    s = iterators[0]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(dataset == "2");
    REQUIRE(records.size() == 0);
    s = iterators[0]->next(dataset, hdr, records);
    REQUIRE(s == StatusCode::NOT_FOUND);
    s = iterators[1]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(dataset == "1");
    REQUIRE(records.size() == 1);
    check();
    s = iterators[1]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(dataset == "2");
    REQUIRE(records.size() == 0);
    s = iterators[1]->next(dataset, hdr, records);
    REQUIRE(s == StatusCode::NOT_FOUND);

    rng = range(0, 100000, 200100);
    s = data->sampleset_range_base(*cache, sampleset, rng,
                              samples, datasets, iterators);
    REQUIRE(s.ok());
    REQUIRE(iterators.size() == 2);
    s = iterators[0]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(dataset == "1");
    REQUIRE(records.size() == 3);
    check();
    s = iterators[0]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(dataset == "2");
    REQUIRE(records.size() == 0);
    s = iterators[0]->next(dataset, hdr, records);
    REQUIRE(s == StatusCode::NOT_FOUND);
    s = iterators[1]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(dataset == "1");
    REQUIRE(records.size() == 2);
    check();
    s = iterators[1]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(dataset == "2");
    REQUIRE(records.size() == 0);
    s = iterators[1]->next(dataset, hdr, records);
    REQUIRE(s == StatusCode::NOT_FOUND);

    rng = range(0, 100000, 300500);
    s = data->sampleset_range_base(*cache, sampleset, rng,
                              samples, datasets, iterators);
    REQUIRE(s.ok());
    REQUIRE(iterators.size() == 3);
    s = iterators[0]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(dataset == "1");
    REQUIRE(records.size() == 3);
    check();
    s = iterators[0]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(dataset == "2");
    REQUIRE(records.size() == 0);
    check();
    s = iterators[0]->next(dataset, hdr, records);
    REQUIRE(s == StatusCode::NOT_FOUND);
    s = iterators[1]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(dataset == "1");
    REQUIRE(records.size() == 5);
    check();
    s = iterators[1]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(dataset == "2");
    REQUIRE(records.size() == 3);
    check();
    s = iterators[1]->next(dataset, hdr, records);
    REQUIRE(s == StatusCode::NOT_FOUND);
    s = iterators[2]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(dataset == "1");
    REQUIRE(records.size() == 2);
    check();
    s = iterators[2]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(dataset == "2");
    REQUIRE(records.size() == 2);
    check();
    s = iterators[2]->next(dataset, hdr, records);
    REQUIRE(s == StatusCode::NOT_FOUND);
}

TEST_CASE("BCFKeyValueData::sampleset_range") {

    // This tests the optimized bucket-based range slicing in BCFKeyValueData

    KeyValueMem::DB db({});
    auto contigs = {make_pair<string,uint64_t>("21", 48129895)};
    REQUIRE(T::InitializeDB(&db, contigs, 25000).ok());
    unique_ptr<T> data;
    REQUIRE(T::Open(&db, data).ok());
    unique_ptr<MetadataCache> cache;
    REQUIRE(MetadataCache::Start(*data, cache).ok());
    set<string> samples_imported;

    Status s = data->import_gvcf(*cache, "1", "test/data/sampleset_range1.gvcf", samples_imported);
    REQUIRE(s.ok());
    s = data->import_gvcf(*cache, "2", "test/data/sampleset_range2.gvcf", samples_imported);
    REQUIRE(s.ok());
    s = data->import_gvcf(*cache, "3", "test/data/sampleset_range3.gvcf", samples_imported);
    REQUIRE(s.ok());

    size_t ct;
    REQUIRE(cache->sample_count(ct).ok());
    REQUIRE(ct == 3);

    string sampleset;
    s = cache->all_samples_sampleset(sampleset);
    REQUIRE(s.ok());

    // We'll slice various ranges and ensure we get the correct number of
    // iterators and that each iterator returns the appropriate data for each
    // dataset

    shared_ptr<const set<string>> samples, datasets;
    vector<unique_ptr<RangeBCFIterator>> iterators;

    range rng(0, 190000, 200000);
    s = data->sampleset_range(*cache, sampleset, rng,
                              samples, datasets, iterators);
    REQUIRE(s.ok());
    // TODO: probably should be 2; the bucket range calculator is a little too
    // inclusive when the range ends exactly on a boundary. Not a bug because
    // records not overlapping the query range are always filtered out at the
    // end.
    REQUIRE(iterators.size() == 3);

    string dataset;
    shared_ptr<const bcf_hdr_t> hdr;
    vector<shared_ptr<bcf1_t>> records, all_records;
    s = iterators[0]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(dataset == "1");
     #define check() s = data->dataset_range(dataset, hdr.get(), range(0,0,1000000), all_records); \
                    REQUIRE(s.ok()); \
                    REQUIRE(records.size() <= count_if(all_records.begin(), all_records.end(), [&](shared_ptr<bcf1_t>& r){return rng.overlaps(r.get());})); \
                    REQUIRE(all_of(records.begin(), records.end(), [&](shared_ptr<bcf1_t>& r){return rng.overlaps(r.get());}))
    REQUIRE(records.size() == 0);
    check();
    REQUIRE(iterators[0]->next(dataset, hdr, records).ok());
    REQUIRE(dataset == "2");
    REQUIRE(records.empty());
    REQUIRE(iterators[0]->next(dataset, hdr, records).ok());
    REQUIRE(dataset == "3");
    REQUIRE(records.empty());
    REQUIRE(iterators[0]->next(dataset, hdr, records) == StatusCode::NOT_FOUND);

    s = iterators[1]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(dataset == "1");
    REQUIRE(records.size() == 3);
    check();
    REQUIRE(iterators[1]->next(dataset, hdr, records).ok());
    REQUIRE(dataset == "2");
    REQUIRE(records.empty());
    REQUIRE(iterators[1]->next(dataset, hdr, records).ok());
    REQUIRE(dataset == "3");
    REQUIRE(records.size() == 3);
    check();
    REQUIRE(iterators[1]->next(dataset, hdr, records) == StatusCode::NOT_FOUND);

    s = iterators[2]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(records.empty());
    REQUIRE(iterators[2]->next(dataset, hdr, records).ok());
    REQUIRE(records.empty());
    REQUIRE(iterators[2]->next(dataset, hdr, records).ok());
    REQUIRE(records.empty());
    REQUIRE(iterators[2]->next(dataset, hdr, records) == StatusCode::NOT_FOUND);

    rng = range(0, 190000, 200050);
    s = data->sampleset_range(*cache, sampleset, rng,
                              samples, datasets, iterators);
    REQUIRE(s.ok());
    REQUIRE(iterators.size() == 3);

    s = iterators[0]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(records.empty());
    REQUIRE(iterators[0]->next(dataset, hdr, records).ok());
    REQUIRE(records.empty());
    REQUIRE(iterators[0]->next(dataset, hdr, records).ok());
    REQUIRE(records.empty());
    REQUIRE(iterators[0]->next(dataset, hdr, records) == StatusCode::NOT_FOUND);

    s = iterators[1]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(records.size() == 3);
    check();
    REQUIRE(iterators[1]->next(dataset, hdr, records).ok());
    REQUIRE(records.empty());
    REQUIRE(iterators[1]->next(dataset, hdr, records).ok());
    REQUIRE(records.size() == 3);
    check();
    REQUIRE(iterators[1]->next(dataset, hdr, records) == StatusCode::NOT_FOUND);

    s = iterators[2]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(records.size() == 1);
    check();
    REQUIRE(iterators[2]->next(dataset, hdr, records).ok());
    REQUIRE(records.empty());
    REQUIRE(iterators[2]->next(dataset, hdr, records).ok());
    REQUIRE(records.size() == 1);
    check();
    REQUIRE(iterators[2]->next(dataset, hdr, records) == StatusCode::NOT_FOUND);

    rng = range(0, 290000, 300050);
    s = data->sampleset_range(*cache, sampleset, rng,
                              samples, datasets, iterators);
    REQUIRE(s.ok());
    REQUIRE(iterators.size() == 3);

    s = iterators[0]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(records.empty());
    REQUIRE(iterators[0]->next(dataset, hdr, records).ok());
    REQUIRE(records.empty());
    REQUIRE(iterators[0]->next(dataset, hdr, records).ok());
    REQUIRE(records.empty());
    REQUIRE(iterators[0]->next(dataset, hdr, records) == StatusCode::NOT_FOUND);

    s = iterators[1]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(records.size() == 3);
    check();
    REQUIRE(iterators[1]->next(dataset, hdr, records).ok());
    REQUIRE(records.size() == 3);
    check();
    REQUIRE(iterators[1]->next(dataset, hdr, records).ok());
    REQUIRE(records.size() == 3);
    check();
    REQUIRE(iterators[1]->next(dataset, hdr, records) == StatusCode::NOT_FOUND);

    s = iterators[2]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(records.size() == 1);
    check();
    REQUIRE(iterators[2]->next(dataset, hdr, records).ok());
    REQUIRE(records.size() == 1);
    check();
    REQUIRE(iterators[2]->next(dataset, hdr, records).ok());
    REQUIRE(records.size() == 1);
    check();
    REQUIRE(iterators[2]->next(dataset, hdr, records) == StatusCode::NOT_FOUND);

    // This query exercises a code path where the KeyValue iterator advances
    // to the end of the database
    rng = range(0, 399999, 400000);
    s = data->sampleset_range(*cache, sampleset, rng,
                              samples, datasets, iterators);
    REQUIRE(s.ok());
    REQUIRE(iterators.size() == 3);

    s = iterators[0]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(records.empty());
    REQUIRE(iterators[0]->next(dataset, hdr, records).ok());
    REQUIRE(records.empty());
    REQUIRE(iterators[0]->next(dataset, hdr, records).ok());
    REQUIRE(records.empty());
    REQUIRE(iterators[0]->next(dataset, hdr, records) == StatusCode::NOT_FOUND);

    s = iterators[1]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(records.size() == 1);
    check();
    REQUIRE(iterators[1]->next(dataset, hdr, records).ok());
    REQUIRE(records.empty());
    REQUIRE(iterators[1]->next(dataset, hdr, records).ok());
    REQUIRE(records.empty());
    REQUIRE(iterators[1]->next(dataset, hdr, records) == StatusCode::NOT_FOUND);

    s = iterators[2]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(records.empty());
    REQUIRE(iterators[2]->next(dataset, hdr, records).ok());
    REQUIRE(records.empty());
    REQUIRE(iterators[2]->next(dataset, hdr, records).ok());
    REQUIRE(records.empty());
    REQUIRE(iterators[2]->next(dataset, hdr, records) == StatusCode::NOT_FOUND);

    // now test with a subset of the samples
    REQUIRE(data->new_sampleset(*cache, "two", set<string>{"HX0002", "HX0003"}).ok());

    rng = range(0, 199899, 199900);
    s = data->sampleset_range(*cache, "two", rng,
                              samples, datasets, iterators);
    REQUIRE(s.ok());
    REQUIRE(*samples == set<string>({"HX0002", "HX0003"}));
    REQUIRE(*datasets == set<string>({"2", "3"}));
    REQUIRE(iterators.size() == 2);

    s = iterators[0]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(dataset == "2");
    REQUIRE(records.empty());
    REQUIRE(iterators[0]->next(dataset, hdr, records).ok());
    REQUIRE(dataset == "3");
    REQUIRE(records.empty());
    REQUIRE(iterators[0]->next(dataset, hdr, records) == StatusCode::NOT_FOUND);

    s = iterators[1]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(dataset == "2");
    REQUIRE(records.empty());
    REQUIRE(iterators[1]->next(dataset, hdr, records).ok());
    REQUIRE(dataset == "3");
    REQUIRE(records.size() == 2);
    check();
    REQUIRE(iterators[1]->next(dataset, hdr, records) == StatusCode::NOT_FOUND);

    // and finally just one sample (tests the code path using sampleset_range_base)
    REQUIRE(data->new_sampleset(*cache, "one", set<string>{"HX0002"}).ok());

    rng = range(0, 299899, 299900);
    s = data->sampleset_range(*cache, "one", rng,
                              samples, datasets, iterators);
    REQUIRE(s.ok());
    REQUIRE(*samples == set<string>({"HX0002"}));
    REQUIRE(*datasets == set<string>({"2"}));
    // we only get one iterator, not iterators for two buckets, because
    // sampleset_range_base is used
    REQUIRE(iterators.size() == 1);

    s = iterators[0]->next(dataset, hdr, records);
    REQUIRE(s.ok());
    REQUIRE(dataset == "2");
    REQUIRE(records.size() == 2);
    check();
    REQUIRE(iterators[0]->next(dataset, hdr, records) == StatusCode::NOT_FOUND);
}

TEST_CASE("BCFKeyValueData compare iterator implementations") {
    // This tests the optimized bucket-based range slicing in BCFKeyValueData
    int nRegions = 13;
    int nIter = 10;
    int lenChrom = 48129895;

    KeyValueMem::DB db({});
    auto contigs = {make_pair<string,uint64_t>("21", lenChrom)};
    REQUIRE(T::InitializeDB(&db, contigs, 1011).ok());
    unique_ptr<T> data;
    REQUIRE(T::Open(&db, data).ok());
    unique_ptr<MetadataCache> cache;
    REQUIRE(MetadataCache::Start(*data, cache).ok());
    set<string> samples_imported;

    Status s = data->import_gvcf(*cache, "1", "test/data/sampleset_range1.gvcf", samples_imported);
    cout << s.str() << endl;
    REQUIRE(s.ok());
    s = data->import_gvcf(*cache, "2", "test/data/sampleset_range2.gvcf", samples_imported);
    REQUIRE(s.ok());
    s = data->import_gvcf(*cache, "3", "test/data/sampleset_range3.gvcf", samples_imported);
    REQUIRE(s.ok());


    string sampleset;
    s = cache->all_samples_sampleset(sampleset);
    REQUIRE(s.ok());

    // split the range [0 ... 1000000] into [nRegions] areas.
    //
    for (int i = 0; i < nIter; i++) {
        int beg = genRandDouble(nRegions) * lenChrom;
        int end = genRandDouble(nRegions) * lenChrom;
        if (end < beg)
            std::swap(beg, end);
        range rng(0, beg, end);
        int rc = compare_query(*data, *cache, sampleset, rng);
        switch (rc) {
        case 1: break; // comparison succeeded
        case 0: REQUIRE(false); // ERROR
        case -1: break;  // Query used too much memory, continue
        }
    }

    cout << "Compared " << nIter << " range queries between the two iterators" << endl;
}

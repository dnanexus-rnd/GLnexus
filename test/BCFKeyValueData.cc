#include <iostream>
#include <map>
#include "BCFKeyValueData.h"
#include "catch.hpp"
using namespace std;
using namespace GLnexus;

using T = BCFKeyValueData<KeyValue::Mem::DB>;

TEST_CASE("BCFKeyValueData construction on improperly initialized database") {
    vector<string> collections = {"header","bcf"};
    KeyValue::Mem::DB db(collections);
    unique_ptr<T> data;
    REQUIRE(T::Open(&db, data) == StatusCode::INVALID);
}

TEST_CASE("BCFKeyValueData initialization") {
    KeyValue::Mem::DB db({});
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
        typename KeyValue::Mem::DB::collection_handle_type coll;
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
        typename KeyValue::Mem::DB::collection_handle_type coll;
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
    KeyValue::Mem::DB db({});
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
    KeyValue::Mem::DB db({});
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

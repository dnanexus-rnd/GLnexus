#include <cstddef>
#include <map>
#include <iostream>
#include <map>
#include <thread>
#include <unistd.h>

#include "KeyValue.h"
#include "BCFKeyValueData.h"
#include "RocksKeyValue.h"

#include "rocksdb/db.h"
#include "rocksdb/slice.h"
#include "rocksdb/options.h"
#include "catch.hpp"

using namespace std;
using namespace GLnexus;

using T = BCFKeyValueData;

std::string kDBPathBase = "/tmp/rocksdb_dir";
static int NUM_LIMIT = 1024 * 1024 * 1024;

// generate a random number in the range [0 .. n-1]
static int gen_rand_number(int n)
{
    static bool firstTime = true;

    // initialization
    if (firstTime) {
        firstTime = false;
        srand (time(NULL));
    }

    int i = rand() % n;
    return i;
}

static const std::string createRandomDBFileName()
{
    std::ostringstream out;
    int rndNum = gen_rand_number(NUM_LIMIT);
    out << kDBPathBase << "_" << std::setfill('0') << std::setw(10) << rndNum;
    return out.str();
}

TEST_CASE("RocksDB open/initialize states") {
    // attempt to open a non-existent database
    string dbPath = createRandomDBFileName();
    std::unique_ptr<KeyValue::DB> db;
    Status s = RocksKeyValue::Open(dbPath, db);
    REQUIRE(s.bad());

    // create the database
    s = RocksKeyValue::Initialize(dbPath, db);
    REQUIRE(s.ok());
    REQUIRE((bool)db);
    unique_ptr<T> data;
    REQUIRE(T::Open(db.get(), data) == StatusCode::INVALID);

    // write something into it
    REQUIRE(db->create_collection("test").ok());
    KeyValue::CollectionHandle coll;
    REQUIRE(db->collection("test",coll).ok());
    REQUIRE(db->put(coll, "foo", "bar").ok());
    std::string v;
    REQUIRE(db->get(coll, "foo", v).ok());
    REQUIRE(v == "bar");

    data.reset();
    db.reset();

    // attempt to initialize an already-existing database
    s = RocksKeyValue::Initialize(dbPath, db);
    REQUIRE(s.bad());
    db.reset();

    // open read-only
    s = RocksKeyValue::Open(dbPath, db, nullptr, RocksKeyValue::OpenMode::READ_ONLY);
    REQUIRE(s.ok());
    REQUIRE(db->collection("test",coll).ok());
    v.clear();
    REQUIRE(db->get(coll, "foo", v).ok());
    REQUIRE(v == "bar");
    REQUIRE(db->put(coll, "foo", "bar").bad());
    REQUIRE(db->put(coll, "bar", "baz").bad());
    db.reset();

    // open in "bulk load" mode
    s = RocksKeyValue::Open(dbPath, db, nullptr, RocksKeyValue::OpenMode::BULK_LOAD);
    REQUIRE(s.ok());
    REQUIRE(db->collection("test",coll).ok());
    v.clear();
    REQUIRE(db->get(coll, "foo", v).ok());
    REQUIRE(v == "bar");
    REQUIRE(db->put(coll, "foo", "bar").ok());
    REQUIRE(db->put(coll, "bar", "baz").ok());
    db.reset();

    RocksKeyValue::destroy(dbPath);
}

TEST_CASE("RocksDB initialization") {
    std::string dbPath = createRandomDBFileName();
    std::unique_ptr<KeyValue::DB> db;
    Status s = RocksKeyValue::Initialize(dbPath, db);
    REQUIRE(s.ok());

    auto contigs = {make_pair<string,uint64_t>("21", 1000000), make_pair<string,uint64_t>("22", 1000001)};
    REQUIRE(T::InitializeDB(db.get(), contigs).ok());
    db.reset();

    s = RocksKeyValue::Open(dbPath, db);
    REQUIRE(s.ok());

    unique_ptr<T> data;
    REQUIRE(T::Open(db.get(), data).ok());

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
        REQUIRE(db->collection("sampleset", coll).ok());
        REQUIRE(db->put(coll, "trio1", "").ok());
        REQUIRE(db->put(coll, "trio1" + null + "fa", "").ok());
        REQUIRE(db->put(coll, "trio1" + null + "mo", "").ok());
        REQUIRE(db->put(coll, "trio1" + null + "ch", "").ok());
        REQUIRE(db->put(coll, "trio2", "").ok());
        REQUIRE(db->put(coll, "trio2" + null + "fa2", "").ok());
        REQUIRE(db->put(coll, "trio2" + null + "mo2", "").ok());
        REQUIRE(db->put(coll, "trio2" + null + "ch2", "").ok());

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
        REQUIRE(db->collection("sample_dataset", coll).ok());
        REQUIRE(db->put(coll, "fa", "trio1").ok());
        REQUIRE(db->put(coll, "mo", "trio1").ok());
        REQUIRE(db->put(coll, "ch", "trio1").ok());
        REQUIRE(db->put(coll, "fa2", "trio2").ok());
        REQUIRE(db->put(coll, "mo2", "trio2").ok());
        REQUIRE(db->put(coll, "ch2", "trio2").ok());

        string dataset;
        REQUIRE(data->sample_dataset("fa", dataset).ok());
        REQUIRE(dataset == "trio1");
        REQUIRE(data->sample_dataset("ch", dataset).ok());
        REQUIRE(dataset == "trio1");
        REQUIRE(data->sample_dataset("mo2", dataset).ok());
        REQUIRE(dataset == "trio2");
        REQUIRE(data->sample_dataset("bogus", dataset) == StatusCode::NOT_FOUND);
    }

    RocksKeyValue::destroy(dbPath);
}

TEST_CASE("RocksDB::import_gvcf") {
    std::unique_ptr<KeyValue::DB> db;
    std::string dbPath = createRandomDBFileName();
    Status s = RocksKeyValue::Initialize(dbPath, db);
    REQUIRE(s.ok());

    auto contigs = {make_pair<string,uint64_t>("21", 48129895)};
    REQUIRE(T::InitializeDB(db.get(), contigs).ok());
    unique_ptr<T> data;
    REQUIRE(T::Open(db.get(), data).ok());
    unique_ptr<MetadataCache> cache;
    REQUIRE(MetadataCache::Start(*data, cache).ok());

    SECTION("NA12878D_HiSeqX.21.10009462-10009469.gvcf") {
        set<string> samples_imported;
        Status s = data->import_gvcf(*cache, "NA12878D", "test/data/NA12878D_HiSeqX.21.10009462-10009469.gvcf", samples_imported);
        REQUIRE(s.ok());

        string dataset;
        REQUIRE(data->sample_dataset("NA12878", dataset).ok());
        REQUIRE(dataset == "NA12878D");
    }

    RocksKeyValue::destroy(dbPath);
}


TEST_CASE("RocksDB::import_gvcf incompatible") {
    std::unique_ptr<KeyValue::DB> db;
    std::string dbPath = createRandomDBFileName();
    Status s = RocksKeyValue::Initialize(dbPath, db);
    REQUIRE(s.ok());

    auto contigs = { make_pair<string,uint64_t>("21", 1000000),
                     make_pair<string,uint64_t>("22", 1000000) };
    REQUIRE(T::InitializeDB(db.get(), contigs).ok());
    unique_ptr<T> data;
    REQUIRE(T::Open(db.get(), data).ok());
    unique_ptr<MetadataCache> cache;
    REQUIRE(MetadataCache::Start(*data, cache).ok());

    SECTION("incompatible contigs") {
        set<string> samples_imported;
        s = data->import_gvcf(*cache, "NA12878D", "test/data/NA12878D_HiSeqX.21.10009462-10009469.gvcf", samples_imported);
        REQUIRE(s == StatusCode::INVALID);
    }

    RocksKeyValue::destroy(dbPath);
}

TEST_CASE("RocksDB BCF retrieval") {
    std::unique_ptr<KeyValue::DB> db;
    std::string dbPath = createRandomDBFileName();
    Status s = RocksKeyValue::Initialize(dbPath, db);
    REQUIRE(s.ok());

    auto contigs = {make_pair<string,uint64_t>("21", 48129895)};
    REQUIRE(T::InitializeDB(db.get(), contigs).ok());
    unique_ptr<T> data;
    REQUIRE(T::Open(db.get(), data).ok());
    unique_ptr<MetadataCache> cache;
    REQUIRE(MetadataCache::Start(*data, cache).ok());
    set<string> samples_imported;

    s = data->import_gvcf(*cache, "NA12878D", "test/data/NA12878D_HiSeqX.21.10009462-10009469.gvcf", samples_imported);
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
        s = data->dataset_range("NA12878D", hdr.get(), range(0, 0, 1000000000), 0, records);
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
        s = data->dataset_range("NA12878D", hdr.get(), range(0, 10009463, 10009466), 0, records);
        REQUIRE(s.ok());
        REQUIRE(records.size() == 2);
        std::shared_ptr<StatsRangeQuery> srq = data->getRangeStats();
        //cout << srq->str() << endl;
        REQUIRE(srq->nBCFRecordsRead == 9);
        REQUIRE(srq->nBCFRecordsInRange == 7);

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
        s = data->dataset_range("NA12878D", hdr.get(), range(0, 0, 1000), 0, records);
        REQUIRE(records.size() == 0);

        s = data->dataset_range("NA12878D", hdr.get(), range(1, 10009463, 10009466), 0, records);
        //REQUIRE(s == StatusCode::NOT_FOUND);
        REQUIRE(records.size() == 0);

        // bogus dataset
        s = data->dataset_range("bogus", hdr.get(), range(1, 10009463, 10009466), 0, records);
        //REQUIRE(s == StatusCode::NOT_FOUND);
        REQUIRE(records.size() == 0);

    }

    RocksKeyValue::destroy(dbPath);
}

TEST_CASE("RocksKeyValue prefix mode") {
    RocksKeyValue::prefix_spec prefix_spec("bcf", BCFKeyValueDataPrefixLength());
    std::unique_ptr<KeyValue::DB> db;
    std::string dbPath = createRandomDBFileName();
    Status s = RocksKeyValue::Initialize(dbPath, db, &prefix_spec);
    REQUIRE(s.ok());

    auto contigs = {make_pair<string,uint64_t>("21", 48129895)};
    REQUIRE(T::InitializeDB(db.get(), contigs).ok());
    unique_ptr<T> data;
    REQUIRE(T::Open(db.get(), data).ok());
    unique_ptr<MetadataCache> cache;
    REQUIRE(MetadataCache::Start(*data, cache).ok());
    set<string> samples_imported;

    s = data->import_gvcf(*cache, "NA12878D", "test/data/NA12878D_HiSeqX.21.10009462-10009469.gvcf", samples_imported);
    REQUIRE(s.ok());

    SECTION("dataset_range") {
        // get all records
        shared_ptr<const bcf_hdr_t> hdr;
        s = data->dataset_header("NA12878D", hdr);
        REQUIRE(s.ok());
        vector<shared_ptr<bcf1_t>> records;

        // subset of records
        s = data->dataset_range("NA12878D", hdr.get(), range(0, 10009463, 10009466), 0, records);
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
        s = data->dataset_range("NA12878D", hdr.get(), range(0, 0, 1000), 0, records);
        REQUIRE(records.size() == 0);

        s = data->dataset_range("NA12878D", hdr.get(), range(1, 10009463, 10009466), 0, records);
        REQUIRE(records.size() == 0);

        // bogus dataset
        s = data->dataset_range("bogus", hdr.get(), range(1, 10009463, 10009466), 0, records);
        REQUIRE(records.size() == 0);

    }

    RocksKeyValue::destroy(dbPath);
}

static std::string getBasename(const std::string &_fname) {
    std::string filename = _fname;

    // remove directory path
    const size_t last_slash_idx = filename.find_last_of("\\/");
    if (std::string::npos != last_slash_idx) {
        filename.erase(0, last_slash_idx + 1);
    }

// Remove extension if present.
    const size_t period_idx = filename.rfind('.');
    if (std::string::npos != period_idx) {
        filename.erase(period_idx);
    }
    return filename;
}

static void importGVCF(T *data,
                       MetadataCache *cache,
                       const std::string &filename,
                       bool flag) {
    set<string> samples_imported;
    std::string dataset = getBasename(filename);
    Status s = data->import_gvcf(*cache, dataset, filename, samples_imported);
    if (flag)
        assert(s.ok());
    else
        assert(s.bad());
    cout << "imported " << filename << endl;;
}

// Query a named dataset. This should work once the data is in the DB;
// it will never be removed.
static void queryDataset(T *data, const std::string &dataset) {
    shared_ptr<const bcf_hdr_t> hdr;
    Status s = data->dataset_header(dataset, hdr);

    vector<shared_ptr<bcf1_t>> records;
    s = data->dataset_range(dataset, hdr.get(), range(0, 1005, 1010), 0,
                            records);
    assert(s.ok());
    assert(records.size() == 0);

    s = data->dataset_range(dataset, hdr.get(), range(0, 2003, 2006), 0,
                            records);
    assert(s.ok());
    assert(records.size() == 3);
}

// Query the '*' sampleset.
static void querySampleset(T *data, int num_iter) {
    Status s;

    for (int i=0; i < num_iter; i++) {
//        if ((i % 10) == 0)
//            cout << "querySampleset " << std::to_string(i) << endl;

        // map global sampleset to samples
        string sampleset;
        s = data->all_samples_sampleset(sampleset);
        assert(s.ok());
        shared_ptr<const set<string> > samples;
        s = data->sampleset_samples(sampleset, samples);
        assert(s.ok());

        // map samples to datasets
        vector<string> datasets;
        for (auto sample : *samples) {
            string ds;
            s = data->sample_dataset(sample, ds);
            assert(s.ok());
            datasets.push_back(ds);
        }

        // query the datasets
        for (auto dataset : datasets)
            queryDataset(data, dataset);

        usleep(3 * 1000); // sleep for a few milliseconds
    }
}

// print the samples
static void printDBSamples(T *data) {
    string sampleset;
    Status s = data->all_samples_sampleset(sampleset);
    assert(s.ok());

    std::shared_ptr<const std::set<std::string> > samples;
    s = data->sampleset_samples(sampleset, samples);
    assert(s.ok());
    //const char* const delim = ", ";
    //std::ostringstream imploded;
    //std::copy(samples->begin(), samples->end(),
    //          std::ostream_iterator<std::string>(imploded, delim));
    //cout << "samples = [" << imploded.str() << "]" << endl;
    assert(samples->size() == 4);
}


// Test concurrent upload, and then concurrent queries.
// Do not test upload failures, and query during upload.
TEST_CASE("Multi-threading") {
    std::unique_ptr<KeyValue::DB> db;
    std::string dbPath = createRandomDBFileName();
    Status s = RocksKeyValue::Initialize(dbPath, db);
    REQUIRE(s.ok());

    auto contigs = {make_pair<string,uint64_t>("21", 48129895)};
    REQUIRE(T::InitializeDB(db.get(), contigs).ok());
    unique_ptr<T> data;
    REQUIRE(T::Open(db.get(), data).ok());
    unique_ptr<MetadataCache> cache;
    REQUIRE(MetadataCache::Start(*data, cache).ok());

    std::vector<shared_ptr<std::thread> > threads;
    std::vector<std::string> files = { "test/data/mt/synthetic_A.21.gvcf",
                                       "test/data/mt/synthetic_B.21.gvcf",
                                       "test/data/mt/synthetic_C.21.gvcf",
                                       "test/data/mt/synthetic_D.21.gvcf" };

    // Spawn a thread for each file to import. All these calls should succeed.
    for (auto filename : files) {
        shared_ptr<thread> thr = make_shared<thread>(importGVCF,
                                                     data.get(), cache.get(), filename,
                                                     true);
        threads.push_back(thr);
    }

    // wait for all import threads to complete
    for (auto thr : threads)
        thr->join();
    threads.clear();
    //std::cout << "All import threads completed\n";

    // sanity check the resulting database
    printDBSamples(data.get());

    // Spawn a thread for each file to import. All these calls should fail, because
    // the data is already in the database.
    //cout << "Trying to add the same data ..." << endl;
    for (auto filename : files) {
        shared_ptr<thread> thr = make_shared<thread>(importGVCF,
                                                     data.get(), cache.get(), filename,
                                                     false);
        threads.push_back(thr);
    }

    // wait for all import threads to complete
    for (auto thr : threads)
        thr->join();
    //cout << "Success" << endl;
    threads.clear();

    // sanity checks
    printDBSamples(data.get());

    //cout << "Running queries serially" << endl;
    {
        // Now run some queries
        for (auto filename : files) {
            queryDataset(data.get(), getBasename(filename));
        }
    }

    //cout << "Running queries in parallel" << endl;
    {
        for (auto filename : files) {
            shared_ptr<thread> thr = make_shared<thread>(queryDataset,
                                                         data.get(),
                                                         getBasename(filename));
            threads.push_back(thr);
        }

        for (auto thr : threads)
            thr->join();
        threads.clear();
    }
    //cout << "Success" << endl;

    RocksKeyValue::destroy(dbPath);
}

//
// N threads that import gVCFs
// M threads that perform queries
//
TEST_CASE("Concurrent import/query") {
    // setup
    std::unique_ptr<KeyValue::DB> db;
    std::string dbPath = createRandomDBFileName();
    Status s = RocksKeyValue::Initialize(dbPath, db);
    REQUIRE(s.ok());

    auto contigs = {make_pair<string,uint64_t>("21", 48129895)};
    REQUIRE(T::InitializeDB(db.get(), contigs).ok());
    unique_ptr<T> data;
    REQUIRE(T::Open(db.get(), data).ok());
    unique_ptr<MetadataCache> cache;
    REQUIRE(MetadataCache::Start(*data, cache).ok());

    std::vector<shared_ptr<std::thread> > threads;
    std::vector<std::string> files = { "test/data/mt/synthetic_A.21.gvcf",
                                       "test/data/mt/synthetic_B.21.gvcf",
                                       "test/data/mt/synthetic_C.21.gvcf",
                                       "test/data/mt/synthetic_D.21.gvcf" };


    // Spawn a thread for each file to import. All these calls should succeed.
    for (auto filename : files) {
        shared_ptr<thread> thr = make_shared<thread>(importGVCF,
                                                     data.get(), cache.get(), filename,
                                                     true);
        threads.push_back(thr);
    }

    // Add a thread that does queries
    shared_ptr<thread> qThr = make_shared<thread>(querySampleset, data.get(), 40);

    // wait for all import threads to complete
    for (auto thr : threads) {
        thr->join();
    }
    qThr->join();

    // cleanup
    RocksKeyValue::destroy(dbPath);
}

#include "BCFKeyValueData.h"
#include "BCFSerialize.h"
#include "yaml-cpp/yaml.h"
#include "vcf.h"
#include "hfile.h"
#include <sstream>
#include <iomanip>
#include <iostream>
#include <assert.h>
#include <math.h>
#include <thread>
#include <mutex>
#include <sys/time.h>
#include "fcmm.hpp"
#include "khash.h"
using namespace std;

namespace GLnexus {


// Memory efficient representation of a bucket range. This could
// be turned into a standard C++ iterator, although, that might be
// a bit of an overkill.
class BucketExtent {
private:
    int rid_ = 0;
    int step_ = 0;
    int bgn_ = 0;
    int end_ = 0;
    int current_ = 0;

    // disable copy and assignment constructors
    BucketExtent(const BucketExtent&);
    BucketExtent& operator=(const BucketExtent&);

public:
    BucketExtent(const range &query, int interval_len) {
        rid_ = query.rid;
        step_ = interval_len;
        bgn_ = std::max((query.beg / interval_len) -1, 0);
        bgn_ *= interval_len;
        end_ = ((query.end / interval_len)) * interval_len;
        current_ = 0;
    }

    range begin() {
        range r = range(rid_, bgn_, bgn_ + step_);
        current_ = bgn_;
        return r;
    }

    range next() {
        current_ += step_;
        range r = range(rid_, current_, current_ + step_);
        return r;
    }

    range end() {
        range r = range(rid_, end_, end_ + step_);
        return r;
    }
};

// Map a range into a set of buckets. The records in the range
// can be found by scanning all the buckets. Note that a BCF record
// is placed in a bucket based on its start position.
// It could start in one bucket, and extend into an adjacent bucket(s).
//
// This class separates out the logic for answering the following questions:
//   1) Which buckets should I scan for this query range?
//   2) Which bucket does a bcf1_t with this range go into?
class BCFBucketRange {
private:
    // disable copy and assignment constructors
    BCFBucketRange(const BCFBucketRange&);
    BCFBucketRange& operator=(const BCFBucketRange&);

public:
    static const size_t PREFIX_LENGTH = 23;
    int interval_len;

    // constructor
    BCFBucketRange(int interval_len) : interval_len(interval_len) {};

    // Given the range of a bucket, produce the key prefix for the bucket.
    // Important: the range must be exactly that of the bucket.
    // BCFBucketRange::bucket below translates an arbitrary range into a
    // bucket's range.
    std::string bucket_prefix(const range& rng) {
        stringstream ss;
        // We add leading zeros to ensure that string lexicographic ordering will sort
        // keys in ascending order.
        ss << setw(3) << setfill('0') << rng.rid
           << setw(10) << setfill('0') << rng.beg
           << setw(10) << setfill('0') << rng.end;
        std::string ans = ss.str();
        assert(ans.size() == PREFIX_LENGTH);
        return ans;
    }

    // Produce the complete key for a bucket (given the prefix) in a dataset
    std::string bucket_key(const std::string& prefix, const std::string& dataset) {
        assert(prefix.size() == 23);
        return prefix+dataset;
    }

    // Same as bucket_key(bucket_prefix(rng), dataset)
    // Important: the range must be exactly that of the bucket.
    // BCFBucketRange::bucket below translates an arbitrary range into a
    // bucket's range.
    std::string bucket_key(const range& rng, const std::string& dataset) {
        return bucket_key(bucket_prefix(rng), dataset);
    }

    // Decompose the key into bucket prefix and dataset
    Status parse_key(const string& key, string& bucket, string& dataset) {
        if (key.size() < PREFIX_LENGTH) {
            return Status::Invalid("BCFBucketRange::parse_key: key too small", key);
        }
        bucket = key.substr(0, PREFIX_LENGTH);
        dataset = key.substr(PREFIX_LENGTH);
        return Status::OK();
    }

    // Given arbitrary [query] range, return a structure describing one or
    // more buckets to search through in order to find all records overlapping
    // query. This may be multiple buckets, even for small query ranges, to
    // account for the possibility of records spanning multiple buckets.
    std::shared_ptr<BucketExtent> scan(const range& query) {
        return make_shared<BucketExtent>(query, interval_len);
    }

    // Which bucket should a BCF record be placed in?
    range bucket(bcf1_t *rec) {
        int bgn = (rec->pos / interval_len) * interval_len;
        return range(rec->rid, bgn, bgn + interval_len);
    }
};

size_t BCFKeyValueDataPrefixLength() { return BCFBucketRange::PREFIX_LENGTH; }


// Metadata that is being added to the DB
struct ActiveMetadata {
    set<string> datasets;
    set<string> samples;

    void erase(const string dataset_, set<string> samples_) {
        datasets.erase(dataset_);
        for (const auto& sample : samples_)
            samples.erase(sample);
    }

    void add(const string dataset_, set<string> samples_) {
        datasets.insert(dataset_);
        for (const auto& sample : samples_)
            samples.insert(sample);

    }
};

// std::hash<string> using the string hash function from htslib
class KStringHash {
public:
    std::size_t operator()(string const& s) const  {
        return (size_t) kh_str_hash_func(s.c_str());
    }
};
using BCFHeaderCache = fcmm::Fcmm<string,shared_ptr<const bcf_hdr_t>,hash<string>,KStringHash>;
// this is not a hard limit but the FCMM performance degrades if it's too low
const size_t BCF_HEADER_CACHE_SIZE = 65536;

// pImpl idiom
struct BCFKeyValueData_body {
    KeyValue::DB* db;
    unique_ptr<BCFHeaderCache> header_cache;
    std::unique_ptr<BCFBucketRange> rangeHelper;
    std::mutex mutex;
    ActiveMetadata amd;
    std::mutex statsMutex;
    StatsRangeQuery statsRq; // statistics for range queries
};

auto collections = { "config", "sampleset", "sample_dataset", "header", "bcf" };

BCFKeyValueData::BCFKeyValueData() = default;
BCFKeyValueData::~BCFKeyValueData() = default;

Status BCFKeyValueData::InitializeDB(KeyValue::DB* db,
                                     const vector<pair<string,size_t>>& contigs,
                                     int interval_len) {
    Status s;

    // create collections
    for (const auto& coll : collections) {
        S(db->create_collection(coll));
    }

    KeyValue::CollectionHandle config;
    S(db->collection("config", config));

    // store contigs
    {
        set<string> prev_contigs;
        YAML::Emitter yaml;
        yaml << YAML::BeginSeq;
        for (const auto& p : contigs) {
            if (prev_contigs.find(p.first) != prev_contigs.end()) {
                return Status::Invalid("duplicate reference contig", p.first);
            }
            yaml << YAML::BeginMap;
            yaml << YAML::Key << p.first;
            yaml << YAML::Value << p.second;
            yaml << YAML::EndMap;
        }
        yaml << YAML::EndSeq;
        S(db->put(config, "contigs", yaml.c_str()));
    }

    // store parameters
    {
        YAML::Emitter yaml;
        yaml << YAML::BeginSeq;
        yaml << YAML::BeginMap;
        yaml << YAML::Key << "interval_len";
        yaml << YAML::Value << interval_len;
        yaml << YAML::EndMap;
        yaml << YAML::EndSeq;
        S(db->put(config, "param", yaml.c_str()));
    }

    // create * sample set, with version number 0
    KeyValue::CollectionHandle sampleset;
    S(db->collection("sampleset", sampleset));
    return db->put(sampleset, "*", "0");
}

Status BCFKeyValueData::Open(KeyValue::DB* db, unique_ptr<BCFKeyValueData>& ans) {
    assert(db != nullptr);

    // check database has been initialized
    KeyValue::CollectionHandle coll;
    for (const auto& collnm : collections) {
        if (db->collection(collnm, coll).bad()) {
            return Status::Invalid("database hasn't been properly initialized");
        }
    }

    ans.reset(new BCFKeyValueData());
    ans->body_.reset(new BCFKeyValueData_body);
    ans->body_->db = db;

    // Read the parameters from the DB
    const char *unexpected = "BCFKeyValueData::Open unexpected YAML";
    Status s;
    vector<pair<string,size_t> > param;
    S(ans->body_->db->collection("config", coll));
    try {
        string param_yaml;
        S(ans->body_->db->get(coll, "param", param_yaml));
        YAML::Node n = YAML::Load(param_yaml);
        if (!n.IsSequence()) {
            return Status::Invalid(unexpected, param_yaml);
        }
        for (const auto& item : n) {
            if (!item.IsMap() || item.size() != 1) {
                return Status::Invalid(unexpected, param_yaml);
            }
            auto m = item.as<map<string,size_t>>();
            assert (m.size() == 1);
            param.push_back(*(m.begin()));
        }
    } catch(YAML::Exception& exn) {
        return Status::Invalid("BCFKeyValueData::Open YAML parse error in param", exn.msg);
    }
    if (param.empty()) {
        return Status::Invalid("database has empty parameter metadata");
    }

    // Sift through parameters, sanity check
    int interval_len = -1;
    for (auto item : param) {
        if (item.first == "interval_len") {
            interval_len = item.second;
            if (interval_len <= 0)
                return Status::Invalid("bad interval length ", std::to_string(interval_len));
        }
    }

    ans->body_->rangeHelper = make_unique<BCFBucketRange>(interval_len);
    ans->body_->header_cache = make_unique<BCFHeaderCache>(BCF_HEADER_CACHE_SIZE);
    return Status::OK();
}

Status BCFKeyValueData::contigs(vector<pair<string,size_t> >& ans) const {
    Status s;
    KeyValue::CollectionHandle coll;
    S(body_->db->collection("config", coll));

    // the contigs entry in config contains a YAML list of contigname-size pairs:
    // - 21: 1000000
    // - 22: 1234567
    // ...

    const char *unexpected = "BCFKeyValueData::contigs unexpected YAML";
    ans.clear();
    try {
        string contigs_yaml;
        S(body_->db->get(coll, "contigs", contigs_yaml));
        YAML::Node n = YAML::Load(contigs_yaml);
        if (!n.IsSequence()) {
            return Status::Invalid(unexpected, contigs_yaml);
        }
        for (const auto& item : n) {
            if (!item.IsMap() || item.size() != 1) {
                return Status::Invalid(unexpected, contigs_yaml);
            }
            auto m = item.as<map<string,size_t>>();
            assert (m.size() == 1);
            ans.push_back(*(m.begin()));
        }
    } catch(YAML::Exception& exn) {
        return Status::Invalid("BCFKeyValueData::contigs YAML parse error", exn.msg);
    }
    if (ans.empty()) {
        return Status::Invalid("database has empty contigs metadata");
    }

    return Status::OK();
}

Status BCFKeyValueData::sampleset_samples(const string& sampleset,
                                          shared_ptr<const set<string> >& ans) const {
    if (sampleset == "*") {
        // * is a special reserved sample set representing all available
        // samples in the database. It must be hidden from callers because
        // it's mutable, while sample sets are supposed to be immutable.
        return Status::NotFound();
    }

    // samplesets collection key scheme:
    // sampleset_id
    // sampleset_id\0sample_1
    // sampleset_id\0sample_2
    // ...
    // sampleset_id\0sample_n
    // next_sampleset
    // next_sampleset\0sample_1
    // ...
    // the corresponding values are empty.

    Status s;
    KeyValue::CollectionHandle coll;
    S(body_->db->collection("sampleset",coll));

    unique_ptr<KeyValue::Iterator> it;
    S(body_->db->iterator(coll, sampleset, it));

    if (!it->valid() || it->key_str() != sampleset) {
        return Status::NotFound("sample set not found", sampleset);
    }

    auto samples = make_shared<set<string>>();
    for (s = it->next(); s.ok() && it->valid(); s = it->next()) {
        auto key = it->key_str();
        size_t nullpos = key.find('\0');
        if (nullpos == string::npos || key.substr(0, nullpos) != sampleset) {
            break;
        }
        samples->insert(key.substr(nullpos+1));
    }
    if (s.bad()) return s;

    ans = samples;

    return Status::OK();
}

Status BCFKeyValueData::sample_dataset(const string& sample, string& ans) const {
    Status s;
    KeyValue::CollectionHandle coll;
    S(body_->db->collection("sample_dataset",coll));
    return body_->db->get(coll, sample, ans);
}

Status BCFKeyValueData::all_samples_sampleset(string& ans) {
    Status s;

    // Get the current * sample set version number.
    KeyValue::CollectionHandle coll;
    S(body_->db->collection("sampleset",coll));
    unique_ptr<KeyValue::Iterator> it;
    S(body_->db->iterator(coll, "*", it));
    if (!it->valid() || it->key_str() != "*") return Status::NotFound("BCFKeyValueData::all_samples_sampleset: improperly initialized database");
    uint64_t version = strtoull(it->value().first, nullptr, 10);
    ans = "*@" + to_string(version); // this is the desired sample set

    // Does the desired sample set exist already? If so, we are done.
    string ignore;
    s = body_->db->get(coll, ans, ignore);
    if (s.ok() || s != StatusCode::NOT_FOUND) {
        return s;
    }

    // Otherwise, continue to read all the samples in the * sample set, and
    // prepare a write batch creating the desired sample set.
    unique_ptr<KeyValue::WriteBatch> wb;
    S(body_->db->begin_writes(wb));
    S(wb->put(coll, ans, string()));
    for (s = it->next(); s.ok() && it->valid(); s = it->next()) {
        auto key = it->key_str();
        size_t nullpos = key.find('\0');
        if (nullpos == string::npos || key.substr(0, nullpos) != "*") {
            break;
        }
        string sample = key.substr(nullpos+1);
        S(wb->put(coll, ans + string(1,'\0') + sample, string()));
    }
    if (s.bad()) return s;
    it.reset();

    // Commit the new sample set iff no other thread has done so in the
    // meantime.
    {
        lock_guard<mutex> lock(body_->mutex);
        s = body_->db->get(coll, ans, ignore);
        if (s != StatusCode::NOT_FOUND) {
            return s;
        }
        return wb->commit();
    }
}

shared_ptr<StatsRangeQuery> BCFKeyValueData::getRangeStats() {
    // return a copy of the current statistics
    std::lock_guard<std::mutex> lock(body_->statsMutex);
    auto statsCopy = make_shared<StatsRangeQuery>(body_->statsRq);
    return statsCopy;
}

Status BCFKeyValueData::dataset_header(const string& dataset,
                                       shared_ptr<const bcf_hdr_t>& hdr) const {
    auto cached = body_->header_cache->end();
    if ((cached = body_->header_cache->find(dataset)) != body_->header_cache->end()) {
        // Return memoized header
        hdr = cached->second;
        assert(hdr);
        return Status::OK();
    }

    // Retrieve the header
    Status s;
    KeyValue::CollectionHandle coll;
    S(body_->db->collection("header",coll));
    string data;
    S(body_->db->get(coll, dataset, data));

    // Parse the header
    shared_ptr<bcf_hdr_t> ans;
    int consumed;
    S(BCFScanner::read_header(data.c_str(), data.size(), consumed, ans));
    hdr = ans;

    // Memoize it
    body_->header_cache->insert(make_pair(dataset, hdr));;
    return Status::OK();
}

// Add a <key,value> pair to the database.
// The key is a concatenation of the dataset name and the chromosome and genomic range.
Status write_bucket(BCFBucketRange& rangeHelper, KeyValue::DB* db,
                    BCFWriter *writer, const string& dataset,
                    const range& rng) {
    // extract the data
    string data;
    Status s;
    S(writer->contents(data));

    // Generate the key
    string key = rangeHelper.bucket_key(rng, dataset);

    // write to the database
    KeyValue::CollectionHandle coll_bcf;
    S(db->collection("bcf", coll_bcf));
    S(db->put(coll_bcf, key, data));
    return Status::OK();
}


// Parse the records and extract those overlapping the query range
static Status scan_bucket(
    const string &dataset,
    const pair<const char*, size_t>& data,
    const bcf_hdr_t* hdr,
    const range& query,
    StatsRangeQuery &srq,
    vector<shared_ptr<bcf1_t> >& records)
{
    Status s;
    unique_ptr<BCFScanner> scanner;
    S(BCFScanner::Open(data.first, data.second, scanner));

    // statistics counter for BCF records
    do {
        srq.nBCFRecordsRead++;
        bool flag;
        S(scanner->overlaps(query, flag));

        if (flag) {
            shared_ptr<bcf1_t> vt;
            S(scanner->read(vt));
            if (bcf_unpack(vt.get(), BCF_UN_ALL) != 0) {
                return Status::IOError("BCFKeyValueData::dataset_bcf bcf_unpack",
                                       dataset + "@" + query.str());
            }
            records.push_back(vt);
        }
        s = scanner->next();
    } while (s.ok());

    if (s != StatusCode::NOT_FOUND) {
        return s;
    }
    return Status::OK();
}

// Search all the buckets that may hold records within the query range. The tricky
// corner case is the bucket immediately before the beginning of the query range.
// It may hold records that start outside the range, but end inside it. We want
// to return those as well.
//
// Return value: list of records that overlap with the query
// range. This list may be empty, if no overlapping records were
// found. If an invalid rid is requested, the return value will be an
// empty list.
Status BCFKeyValueData::dataset_range(const string& dataset,
                                      const bcf_hdr_t* hdr,
                                      const range& query,
                                      vector<shared_ptr<bcf1_t> >& records) {
    Status s;
    records.clear();

    // basic sanity checks
    if (query.rid < 0 || query.beg < 0 || query.end < 0)
        return Status::Invalid("BCFKeyValueData::dataset_bcf: invalid query range", query.str());

    // Retrieve the pertinent DB entries
    KeyValue::CollectionHandle coll;
    S(body_->db->collection("bcf",coll));

    // iterate through the buckets in range
    shared_ptr<BucketExtent> bkExt = body_->rangeHelper->scan(query);

    StatsRangeQuery accu;
    for (range r = bkExt->begin(); r <= bkExt->end(); r = bkExt->next()) {
        //cout << "scanning bucket " << r.str() << endl;
        string key = body_->rangeHelper->bucket_key(r, dataset);
        string data;
        s = body_->db->get(coll, key, data);
        if (s.ok()) {
            S(scan_bucket(dataset, make_pair(data.c_str(), data.size()), hdr, query,
                          accu, records));
        } else if (s != StatusCode::NOT_FOUND) {
            return s;
        }
    }
    accu.nBCFRecordsInRange += records.size();

    // update database statistics
    {
        std::lock_guard<mutex> lock(body_->statsMutex);
        body_->statsRq += accu;
    }

    return Status::OK();
}

// BCFKeyValueData::sampleset_range optimized implementation: if the sample
// set covers >=25% of the samples in the database, produces RangeBCFIterators
// that use underlying KeyValue::Iterators instead of repeated point lookups
// (as in the base implementation). One iterator per underlying storage bucket
// is produced.

class BCFBucketIterator : public RangeBCFIterator {
    BCFData& data_;
    BCFKeyValueData_body& body_;

    bool first_ = true;

    range range_;
    shared_ptr<const set<string>> datasets_;
    set<string>::const_iterator dataset_;

    string bucket_prefix_;
    shared_ptr<KeyValue::Reader> reader_;
    unique_ptr<KeyValue::Iterator> it_;

    StatsRangeQuery stats_;

    Status next_impl(string& dataset, shared_ptr<const bcf_hdr_t>& hdr,
                      vector<shared_ptr<bcf1_t>>& records) {
        // precondition: dataset_ != datasets_.end()

        // pull the desired data set ID (and increment the iterator for the next call)
        dataset = *dataset_++;

        // get the data set header
        Status s;
        S(data_.dataset_header(dataset, hdr));

        if (first_) {
            // first call to next(): begin the iteration at the first dataset
            assert(!it_);
            KeyValue::CollectionHandle coll;
            S(body_.db->collection("bcf",coll));
            S(body_.db->iterator(coll,
                                 body_.rangeHelper->bucket_key(bucket_prefix_, dataset),
                                 it_));
            assert(it_);
            first_ = false;
        }

        records.clear();
        if (!it_ || !it_->valid()) {
            // we've already advanced the KeyValue iterator past the end of
            // the bucket, i.e. the bucket contains no further records for any
            // dataset, so we're now just returning empty results for each
            // dataset.
            return Status::OK();
        }

        // advance the KeyValue iterator to the desired dataset
        string key_dataset;
        for (; s.ok() && it_->valid(); s = it_->next()) {
            string key_prefix;
            S(body_.rangeHelper->parse_key(it_->key_str(), key_prefix, key_dataset));

            if (key_prefix != bucket_prefix_) {
                // we've now advanced past the end of the bucket, so there are
                // no records for this data set (or subsequent data sets)
                assert(key_prefix > bucket_prefix_);
                it_.reset();
                return Status::OK();
            }

            if (key_dataset >= dataset) {
                break;
            }
        }
        if (s.bad()) return s;
        if (!it_->valid()) {
            // wow, we've reached the end of the whole bcf collection
            it_.reset();
            return Status::OK();
        }
        if (key_dataset != dataset) {
            // the database contains no bucket corresponding to this dataset,
            // but might have them for subsequent data sets
            assert(key_dataset > dataset);
            return Status::OK();
        }

        // extract the records overlapping range_
        s = scan_bucket(dataset, it_->value(), hdr.get(), range_, stats_, records);
        if (s.ok()) {
            stats_.nBCFRecordsInRange += records.size();
        }
        return s;
    }

public:
    BCFBucketIterator(BCFData& data, BCFKeyValueData_body& body, range range,
                      const std::string& bucket_prefix, shared_ptr<const set<string>>& datasets,
                      const shared_ptr<KeyValue::Reader>& reader)
        : data_(data), body_(body), range_(range), datasets_(datasets),
          dataset_(datasets->begin()), bucket_prefix_(bucket_prefix), reader_(reader) {}

    virtual ~BCFBucketIterator() {
        lock_guard<mutex> lock(body_.statsMutex);
        body_.statsRq += stats_;
    }

    Status next(string& dataset, shared_ptr<const bcf_hdr_t>& hdr,
                vector<shared_ptr<bcf1_t>>& records) override {
        if (dataset_ == datasets_->end()) {
            // we've finished returning all desired results
            it_.reset();
            return Status::NotFound();
        }

        Status s = next_impl(dataset, hdr, records);
        if (s == StatusCode::NOT_FOUND) {
            // censor NotFound errors so that the caller doesn't misinterpret
            // them as normal EOF.
            return Status::Failure("BCFBucketIterator::next()", s.str());
        }
        return s;
    }
};

Status BCFKeyValueData::sampleset_range(const MetadataCache& metadata, const string& sampleset,
                                        const range& pos,
                                        shared_ptr<const set<string>>& samples,
                                        shared_ptr<const set<string>>& datasets,
                                        vector<unique_ptr<RangeBCFIterator>>& iterators) {
    Status s;

    // resolve samples and datasets
    S(metadata.sampleset_datasets(sampleset, samples, datasets));

    // TODO: dispatch to sampleset_range_base if the sample set is "too
    // small". need a way to know N of the whole database

    // get a KeyValue::Reader so that all iterators read from the same
    // snapshot (this isn't strictly necessary since datasets are immutable,
    // but seems nice to have)
    unique_ptr<KeyValue::Reader> ureader;
    S(body_->db->current(ureader));
    shared_ptr<KeyValue::Reader> reader(move(ureader));

    // create one iterator per bucket
    shared_ptr<BucketExtent> bkExt = body_->rangeHelper->scan(pos);
    iterators.clear();
    for (range r = bkExt->begin(); r <= bkExt->end(); r = bkExt->next()) {
        // Calculate the key prefix for this bucket. The BCFBucketIterator
        // object will use KeyValue::iterator() to position itself to scan all
        // keys with this prefix, stopping upon reaching a key with a
        // different prefix.
        string bucket = body_->rangeHelper->bucket_prefix(r);

        iterators.push_back(make_unique<BCFBucketIterator>
                            (*this, *body_, pos, bucket, datasets, reader));
    }

    return Status::OK();
}

// Provide a way to call the non-optimized base implementation of
// sampleset_range. Mostly for unit testing.
Status BCFKeyValueData::sampleset_range_base(const MetadataCache& metadata, const string& sampleset,
                                             const range& pos,
                                             shared_ptr<const set<string>>& samples,
                                             shared_ptr<const set<string>>& datasets,
                                             vector<unique_ptr<RangeBCFIterator>>& iterators) {
    return BCFData::sampleset_range(metadata, sampleset, pos, samples, datasets, iterators);
}

// test whether a gVCF file is compatible for deposition into the database
static bool gvcf_compatible(const MetadataCache& metadata, const bcf_hdr_t *hdr) {
    Status s;
    const auto& contigs = metadata.contigs();

    // verify contigs match exactly. even the order matters

    int ncontigs = 0;
    const char **contignames = bcf_hdr_seqnames(hdr, &ncontigs);

    if (((uint)ncontigs) != contigs.size()) return false;

    for (int i = 0; i < ncontigs; i++) {
        if (string(contignames[i]) != contigs[i].first) return false;
        // TODO: check contig lengths too
        // REQUIRE(hdr->id[BCF_DT_CTG][0].val != nullptr);
        // REQUIRE(hdr->id[BCF_DT_CTG][0].val->info[0] == 1000000);
    }

    return true;
}

// Sanity-check an individual bcf1_t record before ingestion.
static Status validate_bcf(BCFBucketRange& rangeHelper, const bcf_hdr_t *hdr,
                           bcf1_t *bcf,
                           int prev_rid, int prev_pos) {
    // Check that bcf->rlen is calculated correctly based on POS,END if
    // available or POS,strlen(REF) otherwise
    bcf_info_t *info = bcf_get_info(hdr, bcf, "END");
    if (info) {
        if (info->type != BCF_BT_INT8 && info->type != BCF_BT_INT16 && info->type != BCF_BT_INT32) {
            return Status::Invalid("gVCF record's END field has unexpected type", range(bcf).str());
        }
        if (info->len != 1) {
            return Status::Invalid("gVCF record has multiple END fields", range(bcf).str());
        }
        if (info->v1.i < bcf->pos) {
            return Status::Invalid("gVCF record has END < POS", range(bcf).str());
        }
        if (info->v1.i - bcf->pos != bcf->rlen) {
            return Status::Invalid("gVCF record END-POS doesn't match rlen", range(bcf).str());
        }
    } else {
        if (bcf->rlen != (int) strlen(bcf->d.allele[0])) {
            return Status::Invalid("gVCF rlen doesn't match strlen(REF) (and no END field)", range(bcf).str());
        }
    }

    if (bcf->rlen > rangeHelper.interval_len) {
        return Status::Invalid("gVCF has record that is longer than ",
                               std::to_string(rangeHelper.interval_len));
    }

    // verify record ordering is monotonically increasing within a contig
    if (prev_rid == bcf->rid &&
        prev_pos >= bcf->pos) {
        return Status::Invalid("gVCF record ordering is wrong ",
                               std::to_string(prev_pos) + " >= " + range(bcf).str());
    }

    return Status::OK();
}


// Verify that a VCF file is well formed.
// AND, fill in the [samples_out]
static Status vcf_validate_basic_facts(MetadataCache& metadata,
                                       const string& dataset,
                                       const string& filename,
                                       bcf_hdr_t *hdr,
                                       vcfFile *vcf,
                                       set<string>& samples_out)
{
    if (!hdr) return Status::IOError("reading gVCF header", filename);
    if (!gvcf_compatible(metadata, hdr)) {
        return Status::Invalid("Incompatible gVCF. The reference contigs must match the database configuration exactly.", filename);
    }

    vector<string> samples;
    unsigned n = bcf_hdr_nsamples(hdr);
    if (n == 0) {
        return Status::Invalid("gVCF contains no samples", dataset + " (" + filename + ")");
    }
    for (unsigned i = 0; i < n; i++) {
        string sample(bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i));
        if (sample.size() == 0) {
            return Status::Invalid("gVCF contains empty sample name", dataset + " (" + filename + ")");
        }
        samples.push_back(move(sample));
    }
    samples_out.clear();
    samples_out.insert(samples.begin(), samples.end());
    if (samples.size() != samples_out.size()) {
        return Status::Invalid("gVCF sample names are not unique", dataset + " (" + filename + ")");
    }

    return Status::OK();
}

// Make sure that we the database doesn't already include these datasets and samples.
static Status verify_dataset_and_samples(BCFKeyValueData_body *body_,
                                         MetadataCache& metadata,
                                         const string& dataset,
                                         const string& filename,
                                         const set<string>& samples) {
    Status s;

    KeyValue::CollectionHandle coll;
    S(body_->db->collection("header",coll));
    string data;
    s = body_->db->get(coll, dataset, data);
    if (s.ok()) {
        return Status::Exists("data set already exists",
                              dataset + " (" + filename + ")");
    } else if (s != StatusCode::NOT_FOUND) {
        return s;
    }

    for (const auto& sample : samples) {
        string ignored_dataset;
        s = metadata.sample_dataset(sample, ignored_dataset);
        if (s.ok()) {
            return Status::Exists("sample already exists", sample + " " + dataset + " (" + filename + ")");
        } else if (s != StatusCode::NOT_FOUND) {
            return s;
        }
    }
    return Status::OK();
}

static Status bulk_insert_gvcf_key_values(BCFBucketRange& rangeHelper,
                                          KeyValue::DB* db,
                                          const string& dataset,
                                          const string& filename,
                                          const bcf_hdr_t *hdr,
                                          vcfFile *vcf) {
    Status s;
    unique_ptr<BCFWriter> writer;
    unique_ptr<bcf1_t, void(*)(bcf1_t*)> vt(bcf_init(), &bcf_destroy);
    int prev_pos = -1;
    int prev_rid = -1;

    range bucket(-1, 0, rangeHelper.interval_len);
    int c;
    for(c = bcf_read(vcf, hdr, vt.get());
        c == 0;
        c = bcf_read(vcf, hdr, vt.get())) {
        // Make sure a record is not longer than [BCFBucketRange::interval_len].
        // If this is not true, then we need to compensate in the search
        // routine.
        S(validate_bcf(rangeHelper, hdr, vt.get(), prev_rid, prev_pos));

        // should we start a new bucket?
        if (vt->rid != bucket.rid || vt->pos >= bucket.end) {
            // write old bucket to DB
            if (writer != nullptr && writer->get_num_entries() > 0)
                S(write_bucket(rangeHelper, db, writer.get(), dataset, bucket));

            // start a new in-memory chunk
            writer = nullptr;
            S(BCFWriter::Open(writer));

            // Move to a new genomic range. Round down the start address to
            // a natural multiple of the interval length.
            bucket = rangeHelper.bucket(vt.get());
        }
        S(writer->write(vt.get()));
        prev_rid = vt->rid;
        prev_pos = vt->pos;
    }
    if (c != -1) return Status::IOError("reading from gVCF file", filename);

    // close last bucket
    if (writer != nullptr && writer->get_num_entries() > 0) {
        S(write_bucket(rangeHelper, db, writer.get(), dataset, bucket));
    }
    writer = nullptr;  // clear the memory state
    return Status::OK();
}


// Temporary notes on DB schema, to be moved over to wiki.
//
//  dataset -> header
//       for each dataset, there is a header stored
//  sample -> dataset
//       mapping from sample to dataset, each dataset can store multiple samples.
//
static Status import_gvcf_inner(BCFKeyValueData_body *body_,
                                MetadataCache& metadata,
                                const string& dataset,
                                const string& filename,
                                set<string>& samples_out) {
    Status s;
    unique_ptr<vcfFile, void(*)(vcfFile*)> vcf(bcf_open(filename.c_str(), "r"),
                                               [](vcfFile* f) { bcf_close(f); });
    if (!vcf) return Status::IOError("opening gVCF file", filename);
    unique_ptr<bcf_hdr_t, void(*)(bcf_hdr_t*)> hdr(bcf_hdr_read(vcf.get()), &bcf_hdr_destroy);

    S(vcf_validate_basic_facts(metadata, dataset, filename, hdr.get(), vcf.get(),
                               samples_out));

    // Atomically verify metadata and prepare
    {
        std::lock_guard<std::mutex> lock(body_->mutex);

        // verify uniqueness of data set and samples
        S(verify_dataset_and_samples(body_, metadata, dataset, filename, samples_out));

        // Make sure dataset and samples are not being added by another thread/user
        if (body_->amd.datasets.count(dataset) > 0)
            return Status::Exists("data set is currently being added",
                                  dataset + " (" + filename + ")");

        for (const auto& sample : samples_out)
            if (body_->amd.samples.count(sample) > 0)
                return Status::Exists("sample is currently being added",
                                      sample + " (" + filename + ")");

        // Add to active MD
        body_->amd.add(dataset, samples_out);
    }

    // bulk insert, non atomic
    //
    // Note: we are not dealing at all with mid-flight failures
    S(bulk_insert_gvcf_key_values(*body_->rangeHelper, body_->db,
                                  dataset, filename,
                                  hdr.get(), vcf.get()));

    // Update metadata atomically, now it will point to all the data
    Status retval = Status::Invalid();
    {
        std::lock_guard<std::mutex> lock(body_->mutex);

        // Serialize header into a string
        string hdr_data = BCFWriter::write_header(hdr.get());

        // Get collection handles and current * sample set version number
        KeyValue::CollectionHandle coll_header, coll_sample_dataset, coll_sampleset;
        S(body_->db->collection("header", coll_header));
        S(body_->db->collection("sample_dataset", coll_sample_dataset));
        S(body_->db->collection("sampleset", coll_sampleset));
        string version_str;
        S(body_->db->get(coll_sampleset, "*", version_str));
        uint64_t version = strtoull(version_str.c_str(), nullptr, 10);

        // Store header and metadata (with updated version number)
        unique_ptr<KeyValue::WriteBatch> wb;
        S(body_->db->begin_writes(wb));
        S(wb->put(coll_header, dataset, hdr_data));
        for (const auto& sample : samples_out) {
            // place an entry for this sample in the special "*" sample set
            S(wb->put(coll_sample_dataset, sample, dataset));
            string key = "*" + string(1,'\0') + sample;
            assert(key.size() == sample.size()+2);
            S(wb->put(coll_sampleset, key, string()));
        }
        // update the * sample set version number
        S(wb->put(coll_sampleset, "*", to_string(version+1)));

        // Remove from active metadata
        body_->amd.erase(dataset, samples_out);

        retval = wb->commit();
    }

    // good path
    return retval;
}

Status BCFKeyValueData::import_gvcf(MetadataCache& metadata,
                                    const string& dataset,
                                    const string& filename,
                                    set<string>& samples_out) {
    samples_out.clear(); // hygiene
    Status s = import_gvcf_inner(body_.get(),
                                 metadata,
                                 dataset,
                                 filename,
                                 samples_out);

    if (!s.ok()) {
        // We had a failure, remove from the active metadata
        std::lock_guard<std::mutex> lock(mutex);
        body_->amd.erase(dataset, samples_out);
    }

    return s;
}


} // namespace GLnexus

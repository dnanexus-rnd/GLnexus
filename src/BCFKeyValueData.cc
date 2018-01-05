#include "BCFKeyValueData.h"
#include "BCFSerialize.h"
#include "diploid.h"
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
#include <regex>
#include <endian.h>
using namespace std;

const uint64_t MAX_NUM_CONTIGS_PER_GVCF = 16777216; // 3 bytes wide
const uint64_t MAX_CONTIG_LEN = 1099511627776;      // 5 bytes wide
const uint64_t MAX_RECORD_LEN = 100000;

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
        bgn_ = (query.beg / interval_len) * interval_len;
        end_ = std::max(bgn_, ((query.end-1) / interval_len) * interval_len);
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
    static const size_t PREFIX_LENGTH = 8;
    int interval_len;

    // constructor
    BCFBucketRange(int interval_len) : interval_len(interval_len) {};

    // Given the range of a bucket, produce the key prefix for the bucket.
    // Important: the range must be exactly that of the bucket.
    // BCFBucketRange::bucket below translates an arbitrary range into a
    // bucket's range.
    std::string bucket_prefix(const range& rng) {
        uint64_t rid_be = htobe64(rng.rid);
        uint64_t beg_be = htobe64(rng.beg);
        assert(be64toh(rid_be) < MAX_NUM_CONTIGS_PER_GVCF);
        assert(be64toh(beg_be) < MAX_CONTIG_LEN);
        char buf[8]; static_assert(PREFIX_LENGTH == 8, "assumption");
        memcpy(buf, ((char*)&rid_be)+5, 3);
        memcpy(buf+3, ((char*)&beg_be)+3, 5);
        return string(buf, 8);
    }

    // Produce the complete key for a bucket (given the prefix) in a dataset
    std::string bucket_key(const std::string& prefix, const std::string& dataset) {
        assert(prefix.size() == PREFIX_LENGTH);
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

    // Which bucket does this BCF record start in?
    range bucket(bcf1_t *rec) {
        int bgn = (rec->pos / interval_len) * interval_len;
        return range(rec->rid, bgn, bgn + interval_len);
    }
    // The bucket after [rng], assuming [rng] is a bucket.
    range inc_bucket(range &rng) {
        assert((rng.end - rng.beg) == interval_len);
        return range(rng.rid,
                     rng.beg + interval_len,
                     rng.end + interval_len);
    }

    // Create a ficticious bucket marking the end of a chromosome.
    range bucket_at_end_of_chrom(int rid,
                                 const std::vector<std::pair<std::string,size_t> >&contigs) {
        //const string &contig_name = contigs[rid].first;
        size_t contig_len = contigs[rid].second;
        int bgn = ((contig_len / interval_len) + 2) * interval_len;
        return range(rid, bgn, bgn + interval_len);
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
    atomic<size_t> sample_count; // number of samples in the database. could be
                                 // obtained from the size of the current
                                 // all-samples sampleset, but maintained here
                                 // for convenience.
};

auto collections = { "config", "sampleset", "sample_dataset", "header", "bcf" };

BCFKeyValueData::BCFKeyValueData() = default;
BCFKeyValueData::~BCFKeyValueData() = default;

Status BCFKeyValueData::InitializeDB(KeyValue::DB* db,
                                     const vector<pair<string,size_t>>& contigs,
                                     int interval_len) {
    Status s;

    // some basic sanity checks
    if (contigs.size() > MAX_NUM_CONTIGS_PER_GVCF)
        return Status::Invalid("Too many contigs ", std::to_string(contigs.size()));
    for (const auto& p : contigs) {
        size_t contig_len = p.second;
        if (contig_len > MAX_CONTIG_LEN)
            return Status::Invalid("contig is too long ", string(p.first) + " " + std::to_string(contig_len));
    }

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
    S(db->put(sampleset, "*", "0"));


    // open the database once to trigger otherwise-lazy creation of the * sample set;
    // necessary to ensure the database can subsequently be opened read-only (albeit empty)
    unique_ptr<BCFKeyValueData> nop;
    return BCFKeyValueData::Open(db, nop);
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
        }
    }
    if (interval_len <= 0) {
        return Status::Invalid("Corrupt database; bad interval length ", std::to_string(interval_len));
    }

    ans->body_->rangeHelper = make_unique<BCFBucketRange>(interval_len);
    ans->body_->header_cache = make_unique<BCFHeaderCache>(BCF_HEADER_CACHE_SIZE);

    // initialize sample_count
    string sampleset;
    S(ans->all_samples_sampleset(sampleset));
    shared_ptr<const set<string>> all_samples;
    S(ans->sampleset_samples(sampleset, all_samples));
    ans->body_->sample_count = all_samples->size();

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
    uint64_t version = strtoull(string(it->value().first, it->value().second).c_str(), nullptr, 10);
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

Status BCFKeyValueData::new_sampleset(MetadataCache& metadata,
                                      const string& sampleset,
                                      const set<string>& samples) {
    if (!regex_match(sampleset, regex_id)) {
        return Status::Invalid("BCFKeyValueData::new_sampleset: invalid sample set name");
    }
    if (samples.empty()) {
        return Status::Invalid("BCFKeyValueData::new_sampleset: no samples provided");
    }

    // Prepare a write batch for the new sample set. While doing so, verify
    // that all the samples actually exist. We do this before we take the
    // mutex, thus assuming that samples cannot be deleted.
    Status s;
    string all_samples_id;
    S(metadata.all_samples_sampleset(all_samples_id));
    shared_ptr<const set<string>> all_samples;
    S(metadata.sampleset_samples(all_samples_id, all_samples));

    KeyValue::CollectionHandle coll;
    S(body_->db->collection("sampleset",coll));
    unique_ptr<KeyValue::WriteBatch> wb;
    S(body_->db->begin_writes(wb));
    S(wb->put(coll, sampleset, string()));

    for (const string& sample : samples) {
        if (all_samples->find(sample) == all_samples->end()) {
            return Status::NotFound("BCFKeyValueData::new_sampleset: sample does not exist", sample);
        }
        S(wb->put(coll, sampleset + string(1,'\0') + sample, string()));
    }

    // Now take the mutex
    lock_guard<mutex> lock(body_->mutex);

    // verify the sample set doesn't already exist
    shared_ptr<const set<string>> dummy;
    s = metadata.sampleset_samples(sampleset, dummy);
    if (s.ok()) {
        return Status::Exists("BCFKeyValueData::new_sampleset: sample set already exists", sampleset);
    } else if (s != StatusCode::NOT_FOUND) {
        return s;
    }

    // commit the new sample set
    return wb->commit();
}

Status BCFKeyValueData::sample_count(size_t& ans) const {
    ans = body_->sample_count;
    return Status::OK();
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

// Parse the records and extract those overlapping the query range
// Assumptions: the bucket contains records only from the same reference
// contig as the query, and these records are sorted by increasing beg.
// include_danglers: if false, exclude records whose beg is below that of
//                   the bucket's. This is used when scanning multiple
//                   adjacent buckets; records spanning adjacent buckets
//                   are duplicated in each bucket, so must be
//                   de-duplicated while scanning. In practice you set
//                   include_danglers to true on the first bucket you're
//                   scanning, and false on the rest.
static Status scan_bucket(
    const range& bucket,
    const string &dataset,
    const pair<const char*, size_t>& data,
    const bcf_hdr_t* hdr,
    const range& query,
    bcf_predicate predicate,
    const bool include_danglers,
    StatsRangeQuery &srq,
    vector<shared_ptr<bcf1_t> >& records)
{
    Status s;
    BCFScanner scanner(data.first, data.second);
    range cur_range(-1,-1,-1);
    #ifndef NDEBUG
    range last_range(-1,-1,-1);
    #endif
    while (scanner.valid()) {
        srq.nBCFRecordsRead++; // statistics counter for BCF records
        unsigned n_allele;
        S(scanner.read_range(cur_range, n_allele));
        assert(cur_range.rid == query.rid);
        assert(last_range <= cur_range);
        assert(cur_range.overlaps(bucket));

        if (cur_range.overlaps(query) &&
            (include_danglers || cur_range.beg >= bucket.beg)) {
            shared_ptr<bcf1_t> vt;
            S(scanner.read(vt));
            assert(vt->n_allele == n_allele);

            bool rec_ok = false;
            if (predicate == nullptr) {
                rec_ok = true;
            } else {
                S(predicate(hdr, vt.get(), rec_ok));
            }
            if (rec_ok) {
                if (bcf_unpack(vt.get(), BCF_UN_ALL) != 0 || vt->errcode != 0) {
                    return Status::IOError("BCFKeyValueData bcf_unpack",
                                           dataset + "@" + query.str());
                }
                records.push_back(vt);
            }
        } else if (cur_range.beg >= query.end) {
            // We can quit the scan at this point since the bucket's remaining
            // records all begin past the end of the query range.
            break;
        }
        S(scanner.next());

        #ifndef NDEBUG
        last_range = cur_range;
        #endif
    }

    return s;
}

// Search all the buckets that may hold records within the query range.
//
// Return value: list of records that overlap with the query
// range. This list may be empty, if no overlapping records were
// found. If an invalid rid is requested, the return value will be an
// empty list.
Status BCFKeyValueData::dataset_range(const string& dataset,
                                      const bcf_hdr_t* hdr,
                                      const range& query,
                                      bcf_predicate predicate,
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

    bool first = true;
    StatsRangeQuery accu;
    for (range r = bkExt->begin(); r <= bkExt->end(); r = bkExt->next()) {
        assert(r.overlaps(query));
        string key = body_->rangeHelper->bucket_key(r, dataset);
        string data;
        s = body_->db->get(coll, key, data);
        if (s.ok()) {
            S(scan_bucket(r, dataset, make_pair(data.c_str(), data.size()), hdr, query, predicate,
                          first, accu, records));
        } else if (s != StatusCode::NOT_FOUND) {
            return s;
        }
        first = false;
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
// set covers >=10% of the samples in the database, produces RangeBCFIterators
// that use underlying KeyValue::Iterators instead of repeated point lookups
// (as in the base implementation). One iterator per underlying storage bucket
// is produced.

class BCFBucketIterator : public RangeBCFIterator {
    BCFData& data_;
    BCFKeyValueData_body& body_;

    bool first_ = true;
    bcf_predicate predicate_;
    bool include_danglers_ = true;

    range bucket_, query_;
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

        // extract the records overlapping query_
        s = scan_bucket(bucket_, dataset, it_->value(), hdr.get(), query_, predicate_,
                        include_danglers_, stats_, records);
        if (s.ok()) {
            stats_.nBCFRecordsInRange += records.size();
        }
        return s;
    }

public:
    BCFBucketIterator(BCFData& data, BCFKeyValueData_body& body, const range& query,
                      const range& bucket, const std::string& bucket_prefix,
                      bcf_predicate predicate, bool include_danglers,
                      shared_ptr<const set<string>>& datasets,
                      const shared_ptr<KeyValue::Reader>& reader)
        : data_(data), body_(body), predicate_(predicate), include_danglers_(include_danglers),
          bucket_(bucket), query_(query), datasets_(datasets),
          dataset_(datasets->begin()), bucket_prefix_(bucket_prefix),
          reader_(reader) {}

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
                                        const range& pos, bcf_predicate predicate,
                                        shared_ptr<const set<string>>& samples,
                                        shared_ptr<const set<string>>& datasets,
                                        vector<unique_ptr<RangeBCFIterator>>& iterators) {
    Status s;

    // resolve samples and datasets
    S(metadata.sampleset_datasets(sampleset, samples, datasets));

    // Heuristic: if the desired sample set has fewer than 10% of the samples
    // in the database, then dispatch to the repeated-lookup strategy of
    // sampleset_range_base. This heuristic is wrong if the desired samples
    // are actually contiguous in the database, though.
    size_t total_sample_count;
    S(metadata.sample_count(total_sample_count));
    if (samples->size() == 1 || samples->size() * 10 < total_sample_count) {
        return sampleset_range_base(metadata, sampleset, pos, predicate, samples, datasets, iterators);
    }

    // get a KeyValue::Reader so that all iterators read from the same
    // snapshot (this isn't strictly necessary since datasets are immutable,
    // but seems nice to have)
    unique_ptr<KeyValue::Reader> ureader;
    S(body_->db->current(ureader));
    shared_ptr<KeyValue::Reader> reader(move(ureader));

    // create one iterator per bucket
    bool first = true;
    shared_ptr<BucketExtent> bkExt = body_->rangeHelper->scan(pos);
    iterators.clear();
    for (range r = bkExt->begin(); r <= bkExt->end(); r = bkExt->next()) {
        assert(r.overlaps(pos));
        // Calculate the key prefix for this bucket. The BCFBucketIterator
        // object will use KeyValue::iterator() to position itself to scan all
        // keys with this prefix, stopping upon reaching a key with a
        // different prefix.
        string bucket = body_->rangeHelper->bucket_prefix(r);

        iterators.push_back(make_unique<BCFBucketIterator>
                            (*this, *body_, pos, r, bucket, predicate, first, datasets, reader));
        first = false;
    }

    return Status::OK();
}

// Provide a way to call the non-optimized base implementation of
// sampleset_range. Mostly for unit testing.
Status BCFKeyValueData::sampleset_range_base(const MetadataCache& metadata, const string& sampleset,
                                             const range& pos, bcf_predicate predicate,
                                             shared_ptr<const set<string>>& samples,
                                             shared_ptr<const set<string>>& datasets,
                                             vector<unique_ptr<RangeBCFIterator>>& iterators) {
    return BCFData::sampleset_range(metadata, sampleset, pos, predicate, samples, datasets, iterators);
}

// test whether a gVCF file is compatible for deposition into the database
static bool gvcf_compatible(const MetadataCache& metadata, const bcf_hdr_t *hdr) {
    Status s;
    const auto& contigs = metadata.contigs();

    // verify contigs match exactly. even the order matters

    int ncontigs = 0;
    const char **contignames = bcf_hdr_seqnames(hdr, &ncontigs);
    bool ans = true;

    if (((uint)ncontigs) != contigs.size()) {
        ans = false;
    } else {
        for (int i = 0; i < ncontigs; i++) {
            if (string(contignames[i]) != contigs[i].first ||
                hdr->id[BCF_DT_CTG][i].val->info[0] != contigs[i].second) {
                ans = false;
                break;
            }
        }
    }

    free(contignames);
    return ans;
}

// Sanity-check an individual bcf1_t record before ingestion.
static Status validate_bcf(BCFBucketRange& rangeHelper,
                           const std::vector<std::pair<std::string,size_t> >&contigs,
                           const std::string& filename,
                           const bcf_hdr_t *hdr,
                           bcf1_t *bcf,
                           int prev_rid, int prev_pos) {
    if (bcf_unpack(bcf, BCF_UN_ALL) != 0 || bcf->errcode != 0) {
        return Status::Invalid("invalid VCF record (corrupt format, or undeclared info/format fields)",
                               filename + " " + range(bcf).str(contigs));
    }

    // Check that bcf->rlen is calculated correctly based on POS,END if
    // available or POS,strlen(REF) otherwise
    bcf_info_t *info = bcf_get_info(hdr, bcf, "END");
    if (info) {
        if (info->type != BCF_BT_INT8 && info->type != BCF_BT_INT16 && info->type != BCF_BT_INT32) {
            return Status::Invalid("gVCF record's END field has unexpected type", filename + " " + range(bcf).str(contigs));
        }
        if (info->len != 1) {
            return Status::Invalid("gVCF record has multiple END fields", filename + " " + range(bcf).str(contigs));
        }
        if (info->v1.i < bcf->pos) {
            return Status::Invalid("gVCF record has END < POS", filename + " " + range(bcf).str(contigs));
        }
        if (info->v1.i - bcf->pos != bcf->rlen) {
            return Status::Invalid("gVCF record END-POS doesn't match rlen", filename + " " + range(bcf).str(contigs));
        }
    } else {
        if (bcf->d.allele == nullptr) {
            return Status::Invalid("gVCF allele is null", filename + " " + range(bcf).str(contigs));
        }
        if (bcf->rlen != (int) strlen(bcf->d.allele[0])) {
            return Status::Invalid("gVCF rlen doesn't match strlen(REF) (and no END field)", filename + " " + range(bcf).str(contigs));
        }
    }

    // verify record ordering is non-decreasing within a contig
    if (prev_rid == bcf->rid && prev_pos > bcf->pos) {
        return Status::Invalid("gVCF records are out-of-order ",
                               filename + " " + std::to_string(prev_pos+1) + " >= " + range(bcf).str(contigs));
    }

    // verify record does not go over the length of the contig
    const string &contig_name = contigs[bcf->rid].first;
    size_t contig_len = contigs[bcf->rid].second;
    if (bcf->pos + bcf->rlen > contig_len) {
        return Status::Invalid("gVCF record is longer than contig ",
                               filename + " " + range(bcf).str(contigs) + " " + to_string(contig_len) + " " + contig_name);
    }

    // check that alleles are all distinct, and that all alleles are valid strings
    // of IUPAC nucleotides, except the last ALT allele which may be symbolic.
    set<string> alleles;
    for (int i=0; i < bcf->n_allele; i++) {
        const string allele_i(bcf->d.allele[i]);
        if (!(is_iupac_nucleotides(allele_i) ||
              (i == bcf->n_allele-1 && is_symbolic_allele(allele_i.c_str())))) {
            return Status::Invalid("allele is not a DNA sequence ",
                                filename + " " + allele_i +  " " + range(bcf).str(contigs));
        }
        alleles.insert(allele_i);
    }
    if (bcf->n_allele<1 || alleles.size() != bcf->n_allele) {
        return Status::Invalid("alleles are not distinct ", filename + " " + range(bcf).str(contigs));
    }

    // validate genotypes (all entries either missing or in [0, n_allele))
    if (bcf->n_sample != bcf_hdr_nsamples(hdr)) {
        return Status::Invalid("gVCF record doesn't have expected # of samples", filename + " " + range(bcf).str(contigs));
    }
    htsvecbox<int> gt;
    int nGT = bcf_get_genotypes(hdr, bcf, &gt.v, &gt.capacity);
    if (nGT != 2*bcf->n_sample) {
        return Status::Invalid("gVCF record doesn't have expected # of GT entries", filename + " " + range(bcf).str(contigs));
    }
    for (int i = 0; i < nGT; i++) {
        if (!bcf_gt_is_missing(gt[i]) && (bcf_gt_allele(gt[i]) < 0 || bcf_gt_allele(gt[i]) >= bcf->n_allele)) {
            return Status::Invalid("invalid GT entry in gVCF record", filename + " " + range(bcf).str(contigs));
        }
    }

    // validate genotype likelihoods (n_samples*nGT entries; PL: all entries nonnegative; GL: all entries nonpositive)
    htsvecbox<int32_t> pl;
    int nPL = bcf_get_format_int32(hdr, bcf, "PL", &pl.v, &pl.capacity);
    if (nPL >= 0) {
        if (nPL != bcf->n_sample * diploid::genotypes(bcf->n_allele) && bcf->n_allele > 1) {
            // the exception when bcf->n_allele == 1 accommodates xAtlas
            return Status::Invalid("gVCF record doesn't have expected # of PL entries", filename + " " + range(bcf).str(contigs));
        }
        for (int i = 0; i < nPL; i++) {
            if (pl[i] < 0) {
                return Status::Invalid("negative PL entry in gVCF record", filename + " " + range(bcf).str(contigs));
            }
        }
    }
/*
    htsvecbox<float> gl;
    int nGL = bcf_get_format_float(hdr, bcf, "GL", &gl.v, &gl.capacity);
    if (nGL >= 0) {
        if (nGL != bcf->n_sample * diploid::genotypes(bcf->n_allele)) {
            return Status::Invalid("gVCF record doesn't have expected # of GL entries", filename + " " + range(bcf).str(contigs));
        }
        for (int i = 0; i < nGL; i++) {
            if (gl[i] > 0.0) {
                return Status::Invalid("positive GL entry in gVCF record", filename + " " + range(bcf).str(contigs));
            }
        }
    }
*/
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
        if (!regex_match(sample, regex_id)) {
            return Status::Invalid("gVCF contains invalid sample name", dataset + " (" + filename + ")" + " " + sample);
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

// helper class for bulk_insert_gvcf_key_values: accumulate sizable batches of
// key/value pairs before insertion into the KeyValue database.
// This is to reduce database write lock contention during intense multi-
// threaded bulk loads, as each thread makes fewer larger inserts instead
// of many smaller inserts.
class BulkInsertBuffer {
    const size_t LIMIT = 16777216;
    KeyValue::DB& db_;
    std::unique_ptr<KeyValue::WriteBatch> buf_;
    size_t bufsz_ = 0;

public:
    BulkInsertBuffer(KeyValue::DB& db) : db_(db) {}
    ~BulkInsertBuffer() {
        assert(!buf_);
    }

    Status put(KeyValue::CollectionHandle coll, const std::string& key, const std::string& value) {
        Status s;
        size_t delta = key.size() + value.size() + 32;
        if (bufsz_ + delta >= LIMIT) {
            S(flush());
        }
        if (!buf_) {
            S(db_.begin_writes(buf_));
        }
        S(buf_->put(coll, key, value));
        bufsz_ += delta;
        return Status::OK();
    }

    // make sure to call when finished
    Status flush() {
        if (buf_ && bufsz_) {
            Status s;
            S(buf_->commit());
        }
        buf_.reset();
        bufsz_ = 0;
        return Status::OK();
    }
};

// Add a <key,value> pair to the database.
// The key is a concatenation of the dataset name and the chromosome and genomic range.
static Status write_bucket(BCFBucketRange& rangeHelper, BulkInsertBuffer& db, KeyValue::CollectionHandle& coll_bcf,
                    BCFWriter *writer, unsigned int danglers, const string& dataset,
                    const range& rng,
                    BCFKeyValueData::import_result& rslt) {
    // Generate the key
    string key = rangeHelper.bucket_key(rng, dataset);
    Status s;
    string data;
    //assert(db->get(coll_bcf, key, data) == StatusCode::NOT_FOUND);

    // extract the data
    S(writer->contents(data));

    // write to the database
    S(db.put(coll_bcf, key, data));
    rslt.add_bucket(writer->get_num_entries(), data.size(), danglers);
    return Status::OK();
}

// The code below ingests GVCF records into database buckets. Each bucket is
// inserted as a key/value pair in the database, with the key encoding the
// dataset identifier and genome position. Buckets encompass several kilobases
// of the genome, so most GVCF records fit entirely into a single bucket.
// However, even short records can straddle the a boundary between buckets,
// and reference confidence records and structural variations can be much
// longer than a bucket. To deal with these corner cases, the *danglers
// list* is a buffer of records that dangle from one bucket over to the next.
//
// Short records that straddle a boundary will appear as regular records in
// bucket K, and as danglers in bucket K+1. Long records appear as regular
// records in the first bucket they belong to, and as danglers in all the
// others (they belong to).

// Leave in the danglers list only records that extend beyond [bucket].
static void prune_danglers(vector<shared_ptr<bcf1_t>> &danglers,
                           const range &bucket) {
    auto it = remove_if(
        danglers.begin(),
        danglers.end(),
        [&bucket] (shared_ptr<bcf1_t> &dp) {return range(dp.get()).end <= bucket.end; });
    danglers.erase(it, danglers.end());
}

// convenience macro for some sanity-checking
#ifndef NDEBUG
#define CHECK_DANGLER_BUCKET(pbcf,bkt) { \
    range dangler_rng(pbcf); \
    assert(dangler_rng.overlaps(bkt)); \
    assert(dangler_rng.beg < (bkt).beg); \
    assert(dangler_rng.end > (bkt).beg); \
}
#else
#define CHECK_DANGLER_BUCKET(pbcf,bkt)
#endif

static Status write_danglers_to_in_mem_bucket(vector<shared_ptr<bcf1_t>> &danglers,
                                              BCFWriter *writer,
                                              const range &bucket) {
    Status s;

    if (!danglers.empty()) {
        for (const auto& dp : danglers) {
            if (range(dp.get()).overlaps(bucket)) {
                CHECK_DANGLER_BUCKET(dp.get(), bucket);
                S(writer->write(dp.get()));
            }
        }
        prune_danglers(danglers, bucket);
    }
    return Status::OK();
}

// Write any existing dangling records into the range between [current_bkt]
// and [next_bkt]
static Status write_danglers_between(BCFBucketRange& rangeHelper,
                                     BulkInsertBuffer& db,
                                     KeyValue::CollectionHandle& coll_bcf,
                                     const string& dataset,
                                     range &current_bkt,
                                     BCFKeyValueData::import_result& rslt,
                                     vector<shared_ptr<bcf1_t>> &danglers,
                                     range &next_bkt) {
    Status s;
    unique_ptr<BCFWriter> writer;

    // move to bucket K+1
    range current = rangeHelper.inc_bucket(current_bkt);

    // Subtlety: if this next record is in a bucket other than
    // K+1, but there were danglers from the old bucket K into
    // bucket K+1, then we need to write out bucket K+1
    // (containing only the danglers) before moving on. For
    // very long records, many contiguous filler buckets will
    // be needed.
    while (!danglers.empty() &&
           current < next_bkt) {
        S(BCFWriter::Open(writer));
        for (const auto& dp : danglers) {
            if (range(dp.get()).overlaps(current)) {
                CHECK_DANGLER_BUCKET(dp.get(), current);
                S(writer->write(dp.get()));
            }
        }
        if (writer->get_num_entries() > 0) {
            S(write_bucket(rangeHelper, db, coll_bcf, writer.get(), writer->get_num_entries(),
                           dataset, current, rslt));
        }
        writer.reset();
        prune_danglers(danglers, current);
        current = rangeHelper.inc_bucket(current);
    }

    return Status::OK();
}

static Status bulk_insert_gvcf_key_values(BCFBucketRange& rangeHelper,
                                          MetadataCache& metadata,
                                          KeyValue::DB* db,
                                          const string& dataset,
                                          const string& filename,
                                          const set<range>& range_filter,
                                          const bcf_hdr_t *hdr,
                                          vcfFile *vcf,
                                          BCFKeyValueData::import_result& rslt) {
    Status s;
    unique_ptr<BCFWriter> writer;
    BulkInsertBuffer buffer(*db);
    unique_ptr<bcf1_t, void(*)(bcf1_t*)> vt(bcf_init(), &bcf_destroy);
    int prev_pos = -1;
    int prev_rid = -1;
    vector<shared_ptr<bcf1_t>> danglers;
    unsigned int danglers_written_to_current_bucket = 0;
    // current bucket
    range bucket(-1, 0, rangeHelper.interval_len);
    S(BCFWriter::Open(writer));

    KeyValue::CollectionHandle coll_bcf;
    S(db->collection("bcf", coll_bcf));

    // scan the BCF records
    int c;
    for(c = bcf_read(vcf, hdr, vt.get());
        c == 0 && vt->errcode == 0;
        c = bcf_read(vcf, hdr, vt.get())) {
        range vt_rng(vt.get());
        if (!range_filter.empty()) {
            // Test range filter if applicable. It would be nice to use a
            // tabix index instead.
            if (all_of(range_filter.begin(), range_filter.end(),
                       [&vt_rng](const range& r) { return !r.overlaps(vt_rng); })) {
                continue;
            }
        }

        // A few hard-coded cases where we, reluctantly, skip ingestion
        // VRFromDeletion: accessory information from xAtlas
        // MAX_RECORD_LEN: blows up database (due to repetition across buckets)
        //                 and usually stems from gVCF caller bug anyway
        if (bcf_has_filter(hdr, vt.get(), "VRFromDeletion") == 1 || vt_rng.size() >= MAX_RECORD_LEN) {
            rslt.skipped_records++;
            continue;
        }

        // Check various aspects of the record's validity; e.g. make sure the
        // records are coordinate sorted.
        S(validate_bcf(rangeHelper, metadata.contigs(), filename, hdr, vt.get(), prev_rid, prev_pos));

        // should we start a new bucket?
        if (vt->rid != bucket.rid || vt->pos >= bucket.end) {
            // write old bucket K to DB
            if (writer != nullptr && writer->get_num_entries() > 0) {
                S(write_bucket(rangeHelper, buffer, coll_bcf, writer.get(), danglers_written_to_current_bucket,
                               dataset, bucket, rslt));
            }
            writer.reset();
            range next_bucket = rangeHelper.bucket(vt.get());
            S(write_danglers_between(rangeHelper, buffer, coll_bcf, dataset, bucket, rslt,
                                   danglers, next_bucket));
            bucket = next_bucket;

            // start a new in-memory chunk
            S(BCFWriter::Open(writer));

            // write danglers at the beginning of the new bucket
            danglers_written_to_current_bucket = danglers.size();
            write_danglers_to_in_mem_bucket(danglers, writer.get(), next_bucket);
        }
        // write the record into the bucket
        S(writer->write(vt.get()));
        // if it dangles off the end of the bucket, add it to danglers for
        // inclusion in the next bucket
        if (range(vt.get()).end > bucket.end) {
            auto dangler = shared_ptr<bcf1_t>(bcf_init(), &bcf_destroy);
            bcf_copy(dangler.get(), vt.get());
            assert(range(dangler.get()) == range(vt.get()));
            danglers.push_back(dangler);
        }
        prev_rid = vt->rid;
        prev_pos = vt->pos;
    }
    if (vt->errcode != 0) {
        ostringstream msg;
        msg << filename << " bcf1_t::errcode = " << vt->errcode;
        return Status::IOError("reading from gVCF file", msg.str());
    }
    if (c != -1) return Status::IOError("reading from gVCF file", filename);

    // write out last bucket
    if (writer != nullptr && writer->get_num_entries() > 0) {
        S(write_bucket(rangeHelper, buffer, coll_bcf, writer.get(), danglers_written_to_current_bucket,
                       dataset, bucket, rslt));
    }
    writer.reset();

    // write any last danglers
    range end_bucket = rangeHelper.bucket_at_end_of_chrom(vt->rid, metadata.contigs());
    S(write_danglers_between(rangeHelper, buffer, coll_bcf, dataset, bucket, rslt,
                             danglers, end_bucket));

    return buffer.flush();
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
                                const set<range>& range_filter,
                                BCFKeyValueData::import_result& rslt) {
    Status s;
    unique_ptr<vcfFile, void(*)(vcfFile*)> vcf(bcf_open(filename.c_str(), "r"),
                                               [](vcfFile* f) { bcf_close(f); });
    if (!vcf) return Status::IOError("opening gVCF file", filename);
    unique_ptr<bcf_hdr_t, void(*)(bcf_hdr_t*)> hdr(bcf_hdr_read(vcf.get()), &bcf_hdr_destroy);

    S(vcf_validate_basic_facts(metadata, dataset, filename, hdr.get(), vcf.get(),
                               rslt.samples));

    // Atomically verify metadata and prepare
    {
        std::lock_guard<std::mutex> lock(body_->mutex);

        // verify uniqueness of data set and samples
        S(verify_dataset_and_samples(body_, metadata, dataset, filename, rslt.samples));

        // Make sure dataset and samples are not being added by another thread/user
        if (body_->amd.datasets.count(dataset) > 0)
            return Status::Exists("data set is currently being added",
                                  dataset + " (" + filename + ")");

        for (const auto& sample : rslt.samples)
            if (body_->amd.samples.count(sample) > 0)
                return Status::Exists("sample is currently being added",
                                      sample + " (" + filename + ")");

        // Add to active MD
        body_->amd.add(dataset, rslt.samples);
    }

    // bulk insert, non atomic
    //
    // Note: we are not dealing at all with mid-flight failures
    S(bulk_insert_gvcf_key_values(*body_->rangeHelper, metadata, body_->db,
                                  dataset, filename, range_filter,
                                  hdr.get(), vcf.get(), rslt));

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
        for (const auto& sample : rslt.samples) {
            // place an entry for this sample in the special "*" sample set
            S(wb->put(coll_sample_dataset, sample, dataset));
            string key = "*" + string(1,'\0') + sample;
            assert(key.size() == sample.size()+2);
            S(wb->put(coll_sampleset, key, string()));
        }
        // update the * sample set version number
        S(wb->put(coll_sampleset, "*", to_string(version+1)));

        // Remove from active metadata
        body_->amd.erase(dataset, rslt.samples);

        retval = wb->commit();
        if (retval.ok()) {
            body_->sample_count += rslt.samples.size();
        }
    }

    // good path
    return retval;
}

Status BCFKeyValueData::import_gvcf(MetadataCache& metadata,
                                    const string& dataset,
                                    const string& filename,
                                    const set<range>& range_filter,
                                    import_result& rslt) {
    rslt = import_result(); // hygiene

    if (!regex_match(dataset, regex_id)) {
        return Status::Invalid("BCFKeyValueData::import_gvcf: invalid data set name", dataset);
    }

    Status s = import_gvcf_inner(body_.get(),
                                 metadata,
                                 dataset,
                                 filename,
                                 range_filter,
                                 rslt);

    if (!s.ok()) {
        // We had a failure, remove from the active metadata
        std::lock_guard<std::mutex> lock(body_->mutex);
        body_->amd.erase(dataset, rslt.samples);
    }

    return s;
}


} // namespace GLnexus

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
    int interval_len;

    // constructor
    BCFBucketRange(int interval_len) : interval_len(interval_len) {};

// Create a key for a bucket that includes records from [dataset],
// on chromosome [rid], whose start position is in the genomic range [beg,end).
//
// Note that the end position might be outside the [beg,end) genomic range.
// The search logic needs to compensate for this, by looking at several buckets.
// The key is generated so that data for a genomic range from all samples will be
// located contiguously on disk.
    std::string gen_key(const std::string &dataset, const range& rng) {
        stringstream ss;
        // We add leading zeros to ensure that string lexicographic ordering will sort
        // keys in ascending order.
        ss << setw(3) << setfill('0') << rng.rid
           << setw(10) << setfill('0') << rng.beg
           << setw(10) << setfill('0') << rng.end
           << dataset;
        string key = ss.str();
        return key;
    }

    // Given a [query] range, return a structure describing all the buckets
    // to search through.
    std::shared_ptr<BucketExtent> scan(const std::string &dataset, const range& query) {
        return make_shared<BucketExtent>(query, interval_len);
    }

    // Which bucket should a BCF record be placed in?
    range bucket(const std::string &dataset, bcf1_t *rec) {
        int bgn = (rec->pos / interval_len) * interval_len;
        return range(rec->rid, bgn, bgn + interval_len);
    }
};


// pImpl idiom
struct BCFKeyValueData::body {
    KeyValue::DB* db;
    std::unique_ptr<BCFBucketRange> rangeHelper;
};

auto collections { "config", "sampleset", "sample_dataset", "header", "bcf" };

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

    // create * sample set
    KeyValue::CollectionHandle sampleset;
    S(db->collection("sampleset", sampleset));
    return db->put(sampleset, "*", string());
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
    ans->body_.reset(new body);
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
    Status s;
    KeyValue::CollectionHandle coll;
    S(body_->db->collection("sampleset",coll));

    unique_ptr<KeyValue::Iterator> it;
    S(body_->db->iterator(coll, sampleset, it));

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

    string key, value;
    s = it->next(key, value);
    if (s == StatusCode::NOT_FOUND) {
        return Status::NotFound("sample set not found", sampleset);
    } else if (s.bad()) {
        return s;
    } else if (key != sampleset) {
        return Status::NotFound("sample set not found", sampleset);
    }

    auto samples = make_shared<set<string>>();
    while ((s = it->next(key, value)).ok()) {
        size_t nullpos = key.find('\0');
        if (nullpos == string::npos || key.substr(0, nullpos) != sampleset) {
            break;
        }
        samples->insert(key.substr(nullpos+1));
    }
    if (s.bad() && s != StatusCode::NOT_FOUND) {
        return s;
    }
    ans = samples;

    return Status::OK();

}

Status BCFKeyValueData::sample_dataset(const string& sample, string& ans) const {
    Status s;
    KeyValue::CollectionHandle coll;
    S(body_->db->collection("sample_dataset",coll));
    return body_->db->get(coll, sample, ans);
}

Status BCFKeyValueData::dataset_header(const string& dataset,
                                       shared_ptr<const bcf_hdr_t>& hdr) const {
    // TODO: cache headers

    // Retrieve the header
    Status s;
    KeyValue::CollectionHandle coll;
    S(body_->db->collection("header",coll));
    string data;
    S(body_->db->get(coll, dataset, data));

    // Parse the header
    shared_ptr<bcf_hdr_t> ans;
    int consumed;
    S(BCFReader::read_header(data.c_str(), data.size(), consumed, ans));
    hdr = ans;
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
    string key = rangeHelper.gen_key(dataset, rng);

    // write to the database
    KeyValue::CollectionHandle coll_bcf;
    S(db->collection("bcf", coll_bcf));
    S(db->put(coll_bcf, key, data));
    return Status::OK();
}


// Parse the records and extract those overlapping the query range
static Status scan_bucket(
    const string &dataset,
    const string &key,
    const string &data,
    const bcf_hdr_t* hdr,
    const range& query,
    vector<shared_ptr<bcf1_t> >& records)
{
    Status s;
    unique_ptr<BCFReader> reader;
    S(BCFReader::Open(data.c_str(), data.size(), reader));

    shared_ptr<bcf1_t> vt;
    while ((s = reader->read(vt)).ok()) {
        assert(vt);
        range vt_rng(vt);
        if (query.overlaps(vt_rng)) {
            if (bcf_unpack(vt.get(), BCF_UN_ALL) != 0) {
                return Status::IOError("BCFKeyValueData::dataset_bcf bcf_unpack",
                                       dataset + "@" + query.str());
            }
            records.push_back(vt);
        }
        vt.reset(); // important! otherwise reader overwrites the stored copy.
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
                                      vector<shared_ptr<bcf1_t> >& records) const {
    Status s;
    records.clear();

    // basic sanity checks
    if (query.rid < 0 || query.beg < 0 || query.end < 0)
        return Status::Invalid("BCFKeyValueData::dataset_bcf: invalid query range", query.str());

    // Retrieve the pertinent DB entries
    KeyValue::CollectionHandle coll;
    S(body_->db->collection("bcf",coll));

    // iterate through the buckets in range
    shared_ptr<BucketExtent> bkExt = body_->rangeHelper->scan(dataset, query);

    for (range r = bkExt->begin(); r <= bkExt->end(); r = bkExt->next()) {
        //cout << "scanning bucket " << r.str() << endl;
        string key = body_->rangeHelper->gen_key(dataset, r);
        string data;
        s = body_->db->get(coll, key, data);
        if (s != StatusCode::NOT_FOUND) {
            scan_bucket(dataset, key, data, hdr, query, records);
        }
    }

    return Status::OK();
}

// test whether a gVCF file is compatible for deposition into the database
bool gvcf_compatible(const MetadataCache& metadata, const bcf_hdr_t *hdr) {
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
Status validate_bcf(BCFBucketRange& rangeHelper, const bcf_hdr_t *hdr, bcf1_t *bcf) {
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
        if (bcf->rlen != strlen(bcf->d.allele[0])) {
            return Status::Invalid("gVCF rlen doesn't match strlen(REF) (and no END field)", range(bcf).str());
        }
    }

    if (bcf->rlen > rangeHelper.interval_len) {
        return Status::Invalid("gVCF has record that is longer than ",
                               std::to_string(rangeHelper.interval_len));
    }

    return Status::OK();
}

Status BCFKeyValueData::import_gvcf(MetadataCache& metadata,
                                    const string& dataset,
                                    const string& filename,
                                    set<string>& samples_out) {
    unique_ptr<vcfFile, void(*)(vcfFile*)> vcf(bcf_open(filename.c_str(), "r"),
                                               [](vcfFile* f) { bcf_close(f); });
    if (!vcf) return Status::IOError("opening gVCF file", filename);

    unique_ptr<bcf_hdr_t, void(*)(bcf_hdr_t*)> hdr(bcf_hdr_read(vcf.get()), &bcf_hdr_destroy);
    if (!hdr) return Status::IOError("reading gVCF header", filename);
    if (!gvcf_compatible(metadata, hdr.get())) {
        return Status::Invalid("Incompatible gVCF. The reference contigs must match the database configuration exactly.", filename);
    }

    vector<string> samples;
    unsigned n = bcf_hdr_nsamples(hdr);
    if (n == 0) {
        return Status::Invalid("gVCF contains no samples", dataset + " (" + filename + ")");
    }
    for (unsigned i = 0; i < n; i++) {
        string sample(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, i));
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

    // check uniqueness of data set and samples
    // TODO: design a thread-safe scheme for this
    Status s;
    std::shared_ptr<const bcf_hdr_t> ignored_hdr;
    s = dataset_header(dataset, ignored_hdr);
    if (s != StatusCode::NOT_FOUND) {
        return Status::Exists("data set already exists", dataset + " (" + filename + ")");
    }
    for (const auto& sample : samples) {
        string ignored_dataset;
        s = metadata.sample_dataset(sample, ignored_dataset);
        if (s != StatusCode::NOT_FOUND) {
            return Status::Exists("sample already exists", sample + " " + dataset + " (" + filename + ")");
        }
    }

    // This code assumes that the VCF is sorted by chromosome, and position.
    unique_ptr<BCFWriter> writer;
    unique_ptr<bcf1_t, void(*)(bcf1_t*)> vt(bcf_init(), &bcf_destroy);

    range bucket(-1, 0, body_->rangeHelper->interval_len);
    int c;
    for(c = bcf_read(vcf.get(), hdr.get(), vt.get());
        c == 0;
        c = bcf_read(vcf.get(), hdr.get(), vt.get())) {
        // Make sure a record is not longer than [BCFBucketRange::interval_len].
        // If this is not true, then we need to compensate in the search
        // routine.
        S(validate_bcf(*body_->rangeHelper, hdr.get(), vt.get()));

        // should we start a new bucket?
        if (vt->rid != bucket.rid || vt->pos >= bucket.end) {
            // write old bucket to DB
            if (writer != nullptr && writer->get_num_entries() > 0)
                S(write_bucket(*body_->rangeHelper, body_->db, writer.get(), dataset, bucket));

            // start a new in-memory chunk
            writer = nullptr;
            S(BCFWriter::Open(writer));

            // Move to a new genomic range. Round down the start address to
            // a natural multiple of the interval length.
            bucket = body_->rangeHelper->bucket(dataset, vt.get());
        }
        S(writer->write(vt.get()));
    }
    if (c != -1) return Status::IOError("reading from gVCF file", filename);

    // close last bucket
    if (writer != nullptr && writer->get_num_entries() > 0) {
        S(write_bucket(*body_->rangeHelper, body_->db, writer.get(), dataset, bucket));
    }
    writer = nullptr;  // clear the memory state

    // Serialize header into a string
    string hdr_data = BCFWriter::write_header(hdr.get());

    // Store header and metadata
    KeyValue::CollectionHandle coll_header, coll_sample_dataset, coll_sampleset;
    S(body_->db->collection("header", coll_header));
    S(body_->db->collection("sample_dataset", coll_sample_dataset));
    S(body_->db->collection("sampleset", coll_sampleset));
    unique_ptr<KeyValue::WriteBatch> wb;
    S(body_->db->begin_writes(wb));
    S(wb->put(coll_header, dataset, hdr_data));
    for (const auto& sample : samples) {
        S(wb->put(coll_sample_dataset, sample, dataset));
        string key = "*" + string(1,'\0') + sample;
        assert(key.size() == sample.size()+2);
        S(wb->put(coll_sampleset, key, string()));
    }
    // TODO: invalidate metadata cache for "*"
    return wb->commit();
}

} // namespace GLnexus

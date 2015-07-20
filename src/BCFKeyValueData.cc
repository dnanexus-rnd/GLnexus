#include "BCFKeyValueData.h"
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

// pImpl idiom
struct BCFKeyValueData::body {
    KeyValue::DB* db;
};

auto collections { "config", "sampleset", "sample_dataset", "header", "bcf" };

BCFKeyValueData::BCFKeyValueData() = default;
BCFKeyValueData::~BCFKeyValueData() = default;

Status BCFKeyValueData::InitializeDB(KeyValue::DB* db, const vector<pair<string,size_t>>& contigs) {
    Status s;

    // create collections
    for (const auto& coll : collections) {
        S(db->create_collection(coll));
    }

    // store contigs
    KeyValue::CollectionHandle config;
    S(db->collection("config", config));
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
    return db->put(config, "contigs", yaml.c_str());
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

/// Helper classes to de/serialize BCF records to memory buffers
class BCFReader {
    shared_ptr<const bcf_hdr_t> hdr_;
    vcfFile *bcf_ = nullptr;
    const char* buf_ = nullptr;
    size_t bufsz_;

    BCFReader(const char* buf) : buf_(buf) {}

public:
    // hdr should be null iff the data begins with the header.
    static Status Open(shared_ptr<const bcf_hdr_t> hdr, const char* buf, size_t bufsz, unique_ptr<BCFReader>& ans) {
        ans.reset(new BCFReader(buf));
        ans->bufsz_ = bufsz;

        char fn[4+sizeof(void**)+sizeof(size_t*)];
        char *pfn = fn;

        memcpy(pfn, "mem:", 4);
        *(char***)(pfn+4) = (char**) &(ans->buf_);
        *(size_t**)(pfn+4+sizeof(void**)) = &(ans->bufsz_);

        ans->bcf_ = bcf_open(pfn, "r");
        if (!ans->bcf_) return Status::Failure("BCFReader::Open");

        if (!hdr) {
            hdr = shared_ptr<bcf_hdr_t>(bcf_hdr_read(ans->bcf_), &bcf_hdr_destroy);
            if (!hdr) {
                return Status::IOError("BCFReader::Open bcf_hdr_read");
            }
        }
        ans->hdr_ = hdr;

        return Status::OK();
    }

    virtual ~BCFReader() {
        if (bcf_) {
            bcf_close(bcf_);
        }
    }

    shared_ptr<const bcf_hdr_t> header() const {
        assert(hdr_);
        return hdr_;
    }

    Status read(shared_ptr<bcf1_t>& ans) {
        if (!ans) {
            ans = shared_ptr<bcf1_t>(bcf_init(), &bcf_destroy);
        }
        // TODO use bcf_read; need to override htsFile format detection
        switch (vcf_read(bcf_, hdr_.get(), ans.get())) {
            case 0: return Status::OK();
            case -1: return Status::NotFound();
        }
        return Status::IOError("BCFReader::read bcf_read");
    }
};

class BCFWriter {
    bcf_hdr_t *hdr_;
    vcfFile *bcf_ = nullptr;
    char *buf_ = nullptr;
    size_t bufsz_;

    BCFWriter() = default;

public:

    static Status Open(bcf_hdr_t *hdr, bool write_header, unique_ptr<BCFWriter>& ans) {
        ans.reset(new BCFWriter);
        ans->hdr_ = hdr;

        char fn[4+sizeof(void**)+sizeof(size_t*)];
        char *pfn = fn;

        memcpy(pfn, "mem:", 4);
        *(char***)(pfn+4) = &(ans->buf_);
        *(size_t**)(pfn+4+sizeof(void**)) = &(ans->bufsz_);

        ans->bcf_ = bcf_open(pfn, "wu");
        if (!ans->bcf_) return Status::Failure("BCFWriter::Open");

        if (write_header) {
            if (bcf_hdr_write(ans->bcf_, hdr) != 0) {
                return Status::IOError("BCFWriter::Open bcf_hdr_write");
            }
        }

        return Status::OK();
    }

    virtual ~BCFWriter() {
        if (bcf_) {
            bcf_close(bcf_);
        }
        if (buf_) {
            free(buf_);
        }
    }

    Status write(bcf1_t* x) {
        if (bcf_write(bcf_, hdr_, x) != 0) {
            return Status::IOError("BCFWriter::write");
        }
        return Status::OK();
    }

    Status contents(string& ans) {
        if (hflush(bcf_->fp.hfile) != 0) {
            return Status::IOError("BCFWriter::buffer hflush");
        }
        size_t sz = htell(bcf_->fp.hfile);
        ans.clear();
        if (buf_) {
            ans = string(buf_, sz);
        }
        return Status::OK();
    }
};

Status BCFKeyValueData::dataset_bcf_header(const string& dataset,
                                           shared_ptr<const bcf_hdr_t>& hdr) const {
    // Retrieve the header
    Status s;
    KeyValue::CollectionHandle coll;
    S(body_->db->collection("header",coll));
    string data;
    S(body_->db->get(coll, dataset, data));

    // Parse the header
    unique_ptr<BCFReader> reader;
    S(BCFReader::Open(nullptr, data.c_str(), data.size(), reader));
    hdr = reader->header();
    assert(hdr);
    return Status::OK();
}

Status BCFKeyValueData::dataset_bcf(const string& dataset, const range& pos,
                                    shared_ptr<const bcf_hdr_t>& hdr,
                                    vector<shared_ptr<bcf1_t> >& records) const {
    // TODO cache header...
    Status s;
    S(dataset_bcf_header(dataset, hdr));

    // Retrieve the pertinent DB entries
    // Placeholder: one DB entry per dataset...
    KeyValue::CollectionHandle coll;
    S(body_->db->collection("bcf",coll));
    string data;
    S(body_->db->get(coll, dataset, data));

    // Parse the records and extract those overlapping pos
    unique_ptr<BCFReader> reader;
    S(BCFReader::Open(hdr, data.c_str(), data.size(), reader));

    records.clear();
    shared_ptr<bcf1_t> vt;
    while ((s = reader->read(vt)).ok()) {
        assert(vt);
        if (pos.overlaps(vt.get())) {
            if (bcf_unpack(vt.get(), BCF_UN_ALL) != 0) {
                return Status::IOError("BCFKeyValueData::dataset_bcf bcf_unpack", dataset + "@" + pos.str());
            }
            records.push_back(vt);
        }
        vt.reset(); // important! otherwise reader overwrites the stored copy.
    }
    return Status::OK();
}

// test whether a gVCF file is compatible for deposition into the database
bool gvcf_compatible(const DataCache *cache, const bcf_hdr_t *hdr) {
    Status s;
    auto& contigs = cache->contigs();

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

Status BCFKeyValueData::import_gvcf(const DataCache* cache, const string& dataset, const string& filename) {
    unique_ptr<vcfFile, void(*)(vcfFile*)> vcf(bcf_open(filename.c_str(), "r"),
                                               [](vcfFile* f) { bcf_close(f); });
    if (!vcf) return Status::IOError("opening gVCF file", filename);

    unique_ptr<bcf_hdr_t, void(*)(bcf_hdr_t*)> hdr(bcf_hdr_read(vcf.get()), &bcf_hdr_destroy);
    if (!hdr) return Status::IOError("reading gVCF header", filename);
    if (!gvcf_compatible(cache, hdr.get())) {
        return Status::Invalid("Incompatible gVCF. The reference contigs must match the database configuration exactly.", filename);
    }

    vector<string> samples;
    unsigned n = bcf_hdr_nsamples(hdr);
    for (unsigned i = 0; i < n; i++) {
        samples.push_back(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, i)));
    }
    // TODO check uniqueness of samples and dataset

    // Placeholder functionality: pile all the bcf records into one DB
    // entry...later we need to split it up by genomic range bucket, with
    // records overlapping bucket boundaries duplicated. hence some unique
    // record ID will also be needed for later deduplication...

    Status s;
    unique_ptr<BCFWriter> writer;
    S(BCFWriter::Open(hdr.get(), false, writer));

    unique_ptr<bcf1_t, void(*)(bcf1_t*)> vt(bcf_init(), &bcf_destroy);

    int c;
    for(c = bcf_read(vcf.get(), hdr.get(), vt.get()); c == 0; c = bcf_read(vcf.get(), hdr.get(), vt.get())) {
        S(writer->write(vt.get()));
    }
    if (c != -1) return Status::IOError("reading from gVCF file", filename);

    string data;
    S(writer->contents(data));

    KeyValue::CollectionHandle coll_bcf;
    S(body_->db->collection("bcf", coll_bcf));
    S(body_->db->put(coll_bcf, dataset, data));

    // Serialize header
    writer.release();
    S(BCFWriter::Open(hdr.get(), true, writer));
    S(writer->contents(data));

    // Store header and metadata
    KeyValue::CollectionHandle coll_header, coll_sample_dataset;
    S(body_->db->collection("header", coll_header));
    S(body_->db->collection("sample_dataset", coll_sample_dataset));
    unique_ptr<KeyValue::WriteBatch> wb;
    S(body_->db->begin_writes(wb));
    S(wb->put(coll_header, dataset, data));
    for (const auto& sample : samples) {
        S(wb->put(coll_sample_dataset, sample, dataset));
    }
    return wb->commit();
    
}

} // namespace GLnexus

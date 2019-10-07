#include <iostream>
#include <vcf.h>
#include "service.h"
#include "unifier.h"
#include "genotyper.h"
using namespace std;
using namespace GLnexus;
#include "BCFKeyValueData_utils.h"

// serves data from VCF files in the test/data directory
// x.vcf is loaded as data set "x" with a sample set "x", and a sample set for
// each sample containing just that sample. additionally, the sample set
// "<ALL>" designates all samples across the VCFs.
class VCFData : public Metadata, public BCFData {
    struct vcf_data_t {
        shared_ptr<bcf_hdr_t> header;
        shared_ptr<const vector<pair<string,size_t> > > contigs;
        shared_ptr<const vector<string> > samples;
        vector<shared_ptr<bcf1_t> > records;
    };
    map<string,vcf_data_t> datasets_;
    map<string,string> sample_datasets_;

    VCFData() {}

    static Status load_vcf(const string& path, vcf_data_t& ans, bool validate) {
        unique_ptr<vcfFile, void(*)(vcfFile*)> vcf(bcf_open(path.c_str(), "r"),
                                                   [](vcfFile* f) { bcf_close(f); });
        Status s;
        if (!vcf) {
            return Status::IOError("bcf_open failed", path);
        }

        shared_ptr<bcf_hdr_t> hdr(bcf_hdr_read(vcf.get()), &bcf_hdr_destroy);
        if (!hdr) {
            return Status::IOError("bcf_hdr_read failed", path);
        }

        int ncontigs = 0;
        const char **contignames = bcf_hdr_seqnames(hdr.get(), &ncontigs);
        if (!contignames) return Status::Failure("bcf_hdr_seqnames");
        assert(ncontigs == hdr->n[BCF_DT_CTG]);

        auto contigs = make_shared<vector<pair<string,size_t> > >();
        for (int i = 0; i < ncontigs; ++i)
        {
            size_t sz = 0;
            if (hdr->id[BCF_DT_CTG][i].val) {
                sz = hdr->id[BCF_DT_CTG][i].val->info[0];
            }
            contigs->push_back(make_pair(string(contignames[i]),sz));
        }
        free(contignames);

        int nsamples = bcf_hdr_nsamples(hdr);
        auto samples = make_shared<vector<string> >();
        for (int i = 0; i < nsamples; i++) {
            samples->push_back(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, i)));
        }

        vector<shared_ptr<bcf1_t>> records;
        while (true) {
            shared_ptr<bcf1_t> record(bcf_init(), &bcf_destroy);
            int ret = bcf_read(vcf.get(), hdr.get(), record.get());
            if (ret == -1) {
                break;
            } else if (ret != 0) {
                return Status::IOError("bcf_read", path);
            } else if (bcf_unpack(record.get(),BCF_UN_ALL) != 0) {
                return Status::IOError("bcf_unpack", path);
            }

            if (bcf_has_filter(hdr.get(), record.get(), "VRFromDeletion") == 1) {
                continue;
            }
            bool skip_ingestion = false;
            if (validate) {
                S(validate_bcf(*contigs, path, hdr.get(), record.get(), 0, 0, skip_ingestion));
            }
            if (!skip_ingestion) {
                records.push_back(move(record));
            }
        }

        ans.header = hdr;
        ans.contigs = contigs;
        ans.samples = move(samples);
        ans.records = move(records);

        return Status::OK();
    }

public:
    static Status Open(const set<string> names, unique_ptr<VCFData>& ans,
                       bool validate = true, const string path="test/data/") {
        Status s;
        map<string,vcf_data_t> datasets;
        for (const auto& nm : names) {
            vcf_data_t d;
            string f_path = path + nm;
            s = load_vcf(f_path, d, validate);
            if (s.bad()) return s;

            datasets[nm.substr(0,nm.find_last_of("."))] = d;
        }

        // TODO: verify all headers are mutually "compatible"
        // (contigs, info types etc.)

        map<string,string> sample_datasets;
        for (const auto& ds : datasets) {
            UNPAIR(ds,dataset,dts)
            for (const auto& sample : *dts.samples) {
                if (sample_datasets.find(sample) != sample_datasets.end()) {
                    return Status::Invalid("sample is duplicated", sample);
                }
                sample_datasets[sample] = dataset;
            }
        }
        ans.reset(new VCFData());
        ans->datasets_ = move(datasets);
        ans->sample_datasets_ = sample_datasets;

        return Status::OK();
    }

    Status contigs(vector<pair<string,size_t> >& ans) const override {
        ans = *(datasets_.begin()->second.contigs);
        return Status::OK();
    }

    Status sampleset_samples(const string& sampleset, shared_ptr<const set<string> >& ans) const override {
        auto p = datasets_.find(sampleset);
        if (p != datasets_.end()) {
            ans = make_shared<set<string>>(p->second.samples->begin(), p->second.samples->end());
            return Status::OK();
        }
        auto p2 = sample_datasets_.find(sampleset);
        if (p2 != sample_datasets_.end()) {
            auto sample = make_shared<set<string>>();
            sample->insert(sampleset);
            ans = move(sample);
            return Status::OK();
        }
        if (sampleset == "<ALL>") {
            auto s = make_shared<set<string>>();
            for (const auto& ds : datasets_) {
                s->insert(ds.second.samples->begin(), ds.second.samples->end());
            }
            ans = move(s);
            return Status::OK();
        }
        return Status::NotFound("unknown sample set", sampleset);
    }

    Status all_samples_sampleset(string& ans) override {
        return Status::NotImplemented();
    }

    Status sample_count(size_t& ans) const override {
        ans = sample_datasets_.size();
        return Status::OK();
    }

    Status sample_dataset(const string& sample, string& ans) const override {
        auto p = sample_datasets_.find(sample);
        if (p == sample_datasets_.end()) {
            return Status::NotFound("unknown sample", sample);
        }
        ans = p->second;
        return Status::OK();
    }

    Status dataset_header(const string& dataset, shared_ptr<const bcf_hdr_t>* hdr) override {
        auto p = datasets_.find(dataset);
        if (p == datasets_.end()) {
            return Status::NotFound("unknown data set", dataset);
        }
        *hdr = p->second.header;
        return Status::OK();
    }

    Status dataset_range(const string& dataset, const bcf_hdr_t *hdr,
                         const range& pos, bcf_predicate predicate,
                         vector<shared_ptr<bcf1_t>>* records) override {
        Status s;
        auto p = datasets_.find(dataset);
        if (p == datasets_.end()) {
            return Status::NotFound("unknown data set", dataset);
        }
        records->clear();
        for (const auto& bcf : p->second.records) {
            if (!range(bcf).overlaps(pos))
                continue;

            bool rec_ok = false;
            if (predicate == nullptr) {
                rec_ok = true;
            } else {
                S(predicate(hdr, bcf.get(), rec_ok));
            }
            if (rec_ok)
                records->push_back(bcf);
        }
        return Status::OK();
    }
};

struct TestUtils {
    // Load the text of a VCF file
    static Status load_vcf(const char* vcf, shared_ptr<bcf_hdr_t>& hdr, vector<shared_ptr<bcf1_t>>& records) {
        records.clear();
        hdr.reset();

        const char *tmpfn = "/tmp/GLnexus_unit_test_tmp.vcf";
        ofstream ofs;
        ofs.open(tmpfn, ofstream::out);
        if (!ofs.good()) return Status::IOError("test_vcf1 ofstream");
        ofs << vcf;
        ofs.close();

        shared_ptr<vcfFile> f(vcf_open(tmpfn, "r"), [](vcfFile* h) { vcf_close(h); });
        hdr = shared_ptr<bcf_hdr_t>(bcf_hdr_read(f.get()), &bcf_hdr_destroy);

        int c;
        auto record = shared_ptr<bcf1_t>(bcf_init1(), &bcf_destroy);
        while ((c = bcf_read1(f.get(), hdr.get(), record.get())) == 0) {
            if (bcf_unpack(record.get(), BCF_UN_ALL)) return Status::IOError("test_vcf1 bcf_unpack");
            records.push_back(record);
            record = shared_ptr<bcf1_t>(bcf_init1(), &bcf_destroy);
        }
        if (c != -1) {
            return Status::IOError("test_vcf1 bcf_read1");
        }
        return Status::OK();
    }

    // id. but assuming there's exactly one record
    static Status load_vcf1(const char* vcf, shared_ptr<bcf_hdr_t>& hdr, shared_ptr<bcf1_t>& record) {
        vector<shared_ptr<bcf1_t>> records;
        Status s;
        S(load_vcf(vcf, hdr, records));
        if (records.size() != 1) {
            return Status::Failure("load_test_vcf1: test VCF is supposed to have exactly one record");
        }
        record = records[0];
        return Status::OK();
    }
};

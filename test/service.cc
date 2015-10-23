#include <iostream>
#include <algorithm>
#include <vcf.h>
#include "service.h"
#include "unifier.h"
#include "genotyper.h"
#include "catch.hpp"
using namespace std;
using namespace GLnexus;

// serves data from VCF files in the test/data directory
// x.vcf is loaded as data set "x" with a sample set "x", and a sample set for
// each sample containing just that sample. additionally, the sample set
// "<ALL>" designates all samples across the VCFs.
class VCFData : public Metadata, public BCFData {
    struct vcf_data_t {
        shared_ptr<bcf_hdr_t> header;
        shared_ptr<const vector<string> > samples;
        vector<shared_ptr<bcf1_t> > records;
    };
    map<string,vcf_data_t> datasets_;
    map<string,string> sample_datasets_;

    VCFData() {}

    static Status load_vcf(const string& name, vcf_data_t& ans) {
        string path = "test/data/" + name;
        unique_ptr<vcfFile, void(*)(vcfFile*)> vcf(bcf_open(path.c_str(), "r"),
                                                   [](vcfFile* f) { bcf_close(f); });
        if (!vcf) {
            return Status::IOError("bcf_open failed", path);
        }

        shared_ptr<bcf_hdr_t> hdr(bcf_hdr_read(vcf.get()), &bcf_hdr_destroy);
        if (!hdr) {
            return Status::IOError("bcf_hdr_read failed", path);
        }

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
            records.push_back(move(record));
        }

        ans.header = hdr;
        ans.samples = move(samples);
        ans.records = move(records);

        return Status::OK();
    }

public:
    static Status Open(const set<string> names, unique_ptr<VCFData>& ans) {
        Status s;
        map<string,vcf_data_t> datasets;
        for (const auto& nm : names) {
            vcf_data_t d;
            s = load_vcf(nm, d);
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
        int ncontigs = 0;
        const auto& hdr = datasets_.begin()->second.header;
        const char **contignames = bcf_hdr_seqnames(hdr.get(), &ncontigs);
        if (!contignames) return Status::Failure("bcf_hdr_seqnames");
        assert(ncontigs == hdr->n[BCF_DT_CTG]);

        ans.clear();
        for (int i = 0; i < ncontigs; ++i)
        {
            size_t sz = 0;
            if (hdr->id[BCF_DT_CTG][i].val) {
                sz = hdr->id[BCF_DT_CTG][i].val->info[0];
            }
            ans.push_back(make_pair(string(contignames[i]),sz));
        }
        free(contignames);

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

    Status sample_dataset(const string& sample, string& ans) const override {
        auto p = sample_datasets_.find(sample);
        if (p == sample_datasets_.end()) {
            return Status::NotFound("unknown sample", sample);
        }
        ans = p->second;
        return Status::OK();
    }

    Status dataset_header(const string& dataset, shared_ptr<const bcf_hdr_t>& hdr) const override {
        auto p = datasets_.find(dataset);
        if (p == datasets_.end()) {
            return Status::NotFound("unknown data set", dataset);
        }
        hdr = p->second.header;
        return Status::OK();
    }

    Status dataset_range(const string& dataset, const bcf_hdr_t *hdr, const range& pos, vector<shared_ptr<bcf1_t> >& records) override {
        Status s;
        auto p = datasets_.find(dataset);
        if (p == datasets_.end()) {
            return Status::NotFound("unknown data set", dataset);
        }
        records.clear();
        for (const auto& bcf : p->second.records) {
            if (range(bcf).overlaps(pos)) {
                records.push_back(bcf);
            }
        }
        return Status::OK();
    }
};

// Wrap BCFData to simulate I/O errors for quasi-random testing of errors
class SimFailBCFData : public BCFData {
    BCFData& inner_;
    const size_t fail_every_;
    mutable size_t i_ = 0;
    mutable bool failed_once_ = false;

    SimFailBCFData(BCFData& inner, size_t fail_every)
        : inner_(inner), fail_every_(fail_every) {}

public:
    static Status Open(BCFData& inner, size_t fail_every, unique_ptr<SimFailBCFData>& ans) {
        ans.reset(new SimFailBCFData(inner, fail_every));
        return Status::OK();
    }

    Status dataset_header(const string& dataset, shared_ptr<const bcf_hdr_t>& hdr) const override {
        if (i_++ % fail_every_ == 0) {
            failed_once_ = true;
            return Status::IOError("SIM");
        }
        return inner_.dataset_header(dataset, hdr);
    }

    Status dataset_range(const string& dataset, const bcf_hdr_t *hdr, const range& pos, vector<shared_ptr<bcf1_t> >& records) override {
        if (i_++ % fail_every_ == 0) {
            failed_once_ = true;
            return Status::IOError("SIM");
        }
        return inner_.dataset_range(dataset, hdr, pos, records);
    }

    bool failed_once() { return failed_once_; }
};

TEST_CASE("service::discover_alleles") {
    unique_ptr<VCFData> data;
    Status s = VCFData::Open({"discover_alleles_trio1.vcf", "discover_alleles_trio2.vcf"}, data);
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(*data, *data, svc);
    REQUIRE(s.ok());

    discovered_alleles als;

    SECTION("nonexistent sampleset") {
        s = svc->discover_alleles("bogus", range(0, 0, 1000000), als);
        REQUIRE(s == StatusCode::NOT_FOUND);
    }

    SECTION("trio1") {
        s = svc->discover_alleles("discover_alleles_trio1", range(0, 0, 1099), als);
        REQUIRE(s.ok());

        REQUIRE(als.size() == 7);
        auto p = als.find(allele(range(0, 1000, 1001), "G"));
        REQUIRE(p != als.end());
        REQUIRE(p->second.observation_count == 4);
        p = als.find(allele(range(0, 1010, 1012), "CC"));
        REQUIRE(p != als.end());
        REQUIRE(p->second.observation_count == 3);
    }

    SECTION("trio1 partial") {
        s = svc->discover_alleles("discover_alleles_trio1", range(0, 1009, 1011), als);
        REQUIRE(s.ok());

        REQUIRE(als.size() == 0);

        s = svc->discover_alleles("discover_alleles_trio1", range(0, 1009, 1012), als);
        REQUIRE(s.ok());

        REQUIRE(als.size() == 2);
        REQUIRE(als.find(allele(range(0, 1010, 1012), "CC"))->second.observation_count == 3);
    }

    SECTION("2 trios") {
        s = svc->discover_alleles("<ALL>", range(0, 0, 1099), als);
        REQUIRE(s.ok());

        REQUIRE(als.size() == 10);

        REQUIRE(als.find(allele(range(0, 1000, 1001), "A"))->second.observation_count == 6);
        REQUIRE(als.find(allele(range(0, 1000, 1001), "G"))->second.observation_count == 6);
        REQUIRE(als.find(allele(range(0, 1001, 1002), "A"))->second.observation_count == 6);
        REQUIRE(als.find(allele(range(0, 1001, 1002), "C"))->second.observation_count == 2);
        REQUIRE(als.find(allele(range(0, 1001, 1002), "G"))->second.observation_count == 2);
        REQUIRE(als.find(allele(range(0, 1001, 1002), "T"))->second.observation_count == 2);
        REQUIRE(als.find(allele(range(0, 1010, 1012), "AG"))->second.observation_count == 3);
        REQUIRE(als.find(allele(range(0, 1010, 1012), "CC"))->second.observation_count == 3);
        REQUIRE(als.find(allele(range(0, 1010, 1013), "AGA"))->second.observation_count == 2);
        REQUIRE(als.find(allele(range(0, 1010, 1013), "CCC"))->second.observation_count == 4);
    }

    SECTION("spanning allele") {
        s = svc->discover_alleles("<ALL>", range(1, 1010, 1012), als);
        REQUIRE(s.ok());
        REQUIRE(als.size() == 2);
        REQUIRE(als.find(allele(range(1, 1001, 1016), "AAAAAAAAAAAAAAA")) == als.end());

        s = svc->discover_alleles("<ALL>", range(1, 1001, 1016), als);
        REQUIRE(s.ok());
        REQUIRE(als.size() == 6);
        REQUIRE(als.find(allele(range(1, 1001, 1016), "AAAAAAAAAAAAAAA"))->second.observation_count == 3);
    }

    SECTION("detect inconsistent reference alleles") {
        s = svc->discover_alleles("<ALL>", range(2, 1000, 1004), als);
        REQUIRE_FALSE(s.ok());
        REQUIRE(s == StatusCode::INVALID);
        REQUIRE(s.str().find("allele appears as both REF and ALT") != string::npos);
    }

    SECTION("detect missing reference alleles") {
        /*
        We used to have this test case when the test harness manually
        constructed allele lists, but it's not easy to simulate a missing
        reference allele with test input VCFs...

        als[allele(range(2, 1010, 1011), "A")] = make_pair(false, 1);
        als[allele(range(2, 1010, 1011), "A")] = make_pair(false, 1);

        s = svc->discover_alleles("2trios", range(2, 1010, 1014), als);
        REQUIRE_FALSE(s.ok());
        REQUIRE(s == StatusCode::INVALID);
        REQUIRE(s.str().find("no reference allele") != string::npos);
        */
    }

    SECTION("only return alleles observed in desired samples") {
        // looking at all three samples in the dataset, we should receive all
        // three alleles at this position
        s = svc->discover_alleles("discover_alleles_trio1", range(0, 1001, 1002), als);
        REQUIRE(s.ok());
        REQUIRE(als.size() == 3);

        // looking only at one homozygous ref individual, we should get
        // nothing
        s = svc->discover_alleles("trio1.fa", range(0, 1001, 1002), als);
        REQUIRE(s.ok());
        REQUIRE(als.size() == 0);

        // throw in a previous position as positive control
        s = svc->discover_alleles("trio1.fa", range(0, 1000, 1002), als);
        REQUIRE(s.ok());
        REQUIRE(als.size() == 2);
    }

    SECTION("multiple ranges") {
        vector<range> ranges;
        ranges.push_back(range(0, 1000, 1001));
        ranges.push_back(range(0, 1001, 1002));
        ranges.push_back(range(0, 1010, 1013));
        ranges.push_back(range(1, 1000, 1001));
        ranges.push_back(range(1, 1010, 1012));
        s = svc->discover_alleles("<ALL>", ranges, als);
        REQUIRE(s.ok());

        REQUIRE(als.size() == 14);

        REQUIRE(als.find(allele(range(0, 1000, 1001), "A"))->second.observation_count == 6);
        REQUIRE(als.find(allele(range(0, 1000, 1001), "G"))->second.observation_count == 6);
        REQUIRE(als.find(allele(range(0, 1001, 1002), "A"))->second.observation_count == 6);
        REQUIRE(als.find(allele(range(0, 1001, 1002), "C"))->second.observation_count == 2);
        REQUIRE(als.find(allele(range(0, 1001, 1002), "G"))->second.observation_count == 2);
        REQUIRE(als.find(allele(range(0, 1001, 1002), "T"))->second.observation_count == 2);
        REQUIRE(als.find(allele(range(0, 1010, 1012), "AG"))->second.observation_count == 3);
        REQUIRE(als.find(allele(range(0, 1010, 1012), "CC"))->second.observation_count == 3);
        REQUIRE(als.find(allele(range(0, 1010, 1013), "AGA"))->second.observation_count == 2);
        REQUIRE(als.find(allele(range(0, 1010, 1013), "CCC"))->second.observation_count == 4);

        REQUIRE(als.find(allele(range(1, 1000, 1001), "A"))->second.observation_count == 2);
        REQUIRE(als.find(allele(range(1, 1000, 1001), "AA"))->second.observation_count == 4);
        REQUIRE(als.find(allele(range(1, 1010, 1012), "AG"))->second.observation_count == 3);
        REQUIRE(als.find(allele(range(1, 1010, 1012), "CC"))->second.observation_count == 3);
    }

    SECTION("simulate I/O errors - single") {
        unique_ptr<SimFailBCFData> faildata;
        bool worked = false;

        for (size_t fail_every = 1; fail_every < 100; fail_every++) {
            s = SimFailBCFData::Open(*data, fail_every, faildata);
            REQUIRE(s.ok());

            s = Service::Start(*data, *faildata, svc);
            REQUIRE(s.ok());

            s = svc->discover_alleles("<ALL>", range(0, 0, 1099), als);
            if (faildata->failed_once()) {
                worked = true;
                REQUIRE(s == StatusCode::IO_ERROR);
                REQUIRE(s.str() == "IOError: SIM");
            } else {
                REQUIRE(s.ok());
            }
        }

        REQUIRE(worked);
    }

    SECTION("simulate I/O errors - multi") {
        unique_ptr<SimFailBCFData> faildata;
        bool worked = false;

        for (size_t fail_every = 1; fail_every < 100; fail_every++) {
            s = SimFailBCFData::Open(*data, fail_every, faildata);
            REQUIRE(s.ok());

            s = Service::Start(*data, *faildata, svc);
            REQUIRE(s.ok());

            vector<range> ranges;
            ranges.push_back(range(0, 1000, 1001));
            ranges.push_back(range(0, 1001, 1002));
            ranges.push_back(range(0, 1010, 1013));
            ranges.push_back(range(1, 1000, 1001));
            ranges.push_back(range(1, 1010, 1012));
            s = svc->discover_alleles("<ALL>", ranges, als);
            if (faildata->failed_once()) {
                worked = true;
                REQUIRE(s == StatusCode::IO_ERROR);
                REQUIRE(s.str() == "IOError: SIM");
            } else {
                REQUIRE(s.ok());
            }
        }

        REQUIRE(worked);
    }
}


TEST_CASE("service::discover_alleles gVCF") {
    unique_ptr<VCFData> data;
    Status s = VCFData::Open({"NA12878D_HiSeqX.21.10009462-10009469.gvcf", "NA12878D_HiSeqX.21.10009462-10009469.bogus.gvcf"}, data);
    cout << s.str() << endl;
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(*data, *data, svc);
    REQUIRE(s.ok());

    discovered_alleles als;

    SECTION("exclude symbolic alleles") {
        s = svc->discover_alleles("<ALL>", range(0, 10000000, 10010000), als);

        REQUIRE(als.size() == 2);
        auto p = als.find(allele(range(0, 10009463, 10009465), "TA"));
        REQUIRE(p != als.end());
        REQUIRE(p->second.observation_count == 1);
        p = als.find(allele(range(0, 10009463, 10009465), "T"));
        REQUIRE(p != als.end());
        REQUIRE(p->second.observation_count == 1);
    }

    SECTION("exclusion/detection of bogus alleles") {
        s = svc->discover_alleles("<ALL>", range(1, 10009463, 10009465), als);

        REQUIRE(als.size() == 2);
        auto p = als.find(allele(range(1, 10009463, 10009465), "TA"));
        REQUIRE(p != als.end());
        REQUIRE(p->second.observation_count == 1);
        p = als.find(allele(range(1, 10009463, 10009465), "T"));
        REQUIRE(p != als.end());
        REQUIRE(p->second.observation_count == 1);

        s = svc->discover_alleles("<ALL>", range(1, 10009465, 10009466), als);
        REQUIRE(s == StatusCode::INVALID);
    }
}

TEST_CASE("unified_sites placeholder") {
    unique_ptr<VCFData> data;
    Status s = VCFData::Open({"discover_alleles_trio1.vcf", "discover_alleles_trio2.vcf"}, data);
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(*data, *data, svc);
    REQUIRE(s.ok());

    discovered_alleles als;

    SECTION("trio1") {
        s = svc->discover_alleles("discover_alleles_trio1", range(0, 0, 1000000), als);
        REQUIRE(s.ok());

        vector<unified_site> sites;
        s = unified_sites(als, sites);
        REQUIRE(s.ok());

        vector<pair<string,size_t> > contigs;
        REQUIRE(data->contigs(contigs).ok());

        REQUIRE(sites.size() == 6);

        REQUIRE(sites[0].pos == range(0,1000,1001));
        REQUIRE(sites[0].alleles.size() == 2);
        REQUIRE(sites[0].alleles[0] == "A");
        REQUIRE(sites[0].alleles[1] == "G");
        REQUIRE(sites[0].unification[make_pair(range(0,1000,1001),string("A"))] == 0);
        REQUIRE(sites[0].unification[make_pair(range(0,1000,1001),string("G"))] == 1);
        REQUIRE(sites[0].observation_count.size() == 2);
        REQUIRE(sites[0].observation_count[0] == 2);
        REQUIRE(sites[0].observation_count[1] == 4);

        REQUIRE(sites[1].pos == range(0,1001,1002));
        REQUIRE(sites[1].alleles.size() == 3);
        REQUIRE(sites[1].alleles[0] == "C");
        REQUIRE(sites[1].alleles[1] == "G");
        REQUIRE(sites[1].alleles[2] == "T");
        REQUIRE(sites[1].unification[make_pair(range(0,1001,1002),string("C"))] == 0);
        REQUIRE(sites[1].unification[make_pair(range(0,1001,1002),string("G"))] == 1);
        REQUIRE(sites[1].unification[make_pair(range(0,1001,1002),string("T"))] == 2);
        REQUIRE(sites[1].observation_count.size() == 3);
        REQUIRE(sites[1].observation_count[0] == 2);
        REQUIRE(sites[1].observation_count[1] == 2);
        REQUIRE(sites[1].observation_count[2] == 2);

        REQUIRE(sites[2].pos == range(0,1010,1012));
        REQUIRE(sites[2].alleles.size() == 2);
        REQUIRE(sites[2].alleles[0] == "CC");
        REQUIRE(sites[2].alleles[1] == "AG");
        REQUIRE(sites[2].unification[make_pair(range(0,1010,1012),string("CC"))] == 0);
        REQUIRE(sites[2].unification[make_pair(range(0,1010,1012),string("AG"))] == 1);
        REQUIRE(sites[2].observation_count.size() == 2);
        REQUIRE(sites[2].observation_count[0] == 3);
        REQUIRE(sites[2].observation_count[1] == 3);

        REQUIRE(sites[3].pos == range(0,1100,1101));
        REQUIRE(sites[3].alleles.size() == 2);
        REQUIRE(sites[3].alleles[0] == "C");
        REQUIRE(sites[3].alleles[1] == "A");
        REQUIRE(sites[3].unification[make_pair(range(0,1100,1101),string("C"))] == 0);
        REQUIRE(sites[3].unification[make_pair(range(0,1100,1101),string("A"))] == 1);
        REQUIRE(sites[3].observation_count.size() == 2);
        REQUIRE(sites[3].observation_count[0] == 3);
        REQUIRE(sites[3].observation_count[1] == 3);

        REQUIRE(sites[4].pos == range(0,1102,1103));
        REQUIRE(sites[4].alleles.size() == 2);
        REQUIRE(sites[4].alleles[0] == "C");
        REQUIRE(sites[4].alleles[1] == "G");
        REQUIRE(sites[4].unification[make_pair(range(0,1102,1103),string("C"))] == 0);
        REQUIRE(sites[4].unification[make_pair(range(0,1102,1103),string("G"))] == 1);
        REQUIRE(sites[4].observation_count.size() == 2);
        REQUIRE(sites[4].observation_count[0] == 3);
        REQUIRE(sites[4].observation_count[1] == 3);

        REQUIRE(sites[5].pos == range(0,1200,1201));
        REQUIRE(sites[5].alleles.size() == 2);
        REQUIRE(sites[5].alleles[0] == "C");
        REQUIRE(sites[5].alleles[1] == "A");
        REQUIRE(sites[5].unification[make_pair(range(0,1200,1201),string("C"))] == 0);
        REQUIRE(sites[5].unification[make_pair(range(0,1200,1201),string("A"))] == 1);
        REQUIRE(sites[5].observation_count.size() == 2);
        REQUIRE(sites[5].observation_count[0] == 3);
        REQUIRE(sites[5].observation_count[1] == 3);
    }

    SECTION("2 trios") {
        s = svc->discover_alleles("<ALL>", range(0, 0, 1000000), als);
        REQUIRE(s.ok());

        vector<unified_site> sites;
        s = unified_sites(als, sites);
        REQUIRE(s.ok());

        vector<pair<string,size_t> > contigs;
        REQUIRE(data->contigs(contigs).ok());

        REQUIRE(sites.size() == 6);

        REQUIRE(sites[0].pos == range(0,1000,1001));
        REQUIRE(sites[0].alleles.size() == 2);
        REQUIRE(sites[0].alleles[0] == "A");
        REQUIRE(sites[0].alleles[1] == "G");
        REQUIRE(sites[0].unification[make_pair(range(0,1000,1001),string("A"))] == 0);
        REQUIRE(sites[0].unification[make_pair(range(0,1000,1001),string("G"))] == 1);
        REQUIRE(sites[0].observation_count.size() == 2);
        REQUIRE(sites[0].observation_count[0] == 6);
        REQUIRE(sites[0].observation_count[1] == 6);

        REQUIRE(sites[1].pos == range(0,1001,1002));
        REQUIRE(sites[1].alleles.size() == 4);
        REQUIRE(sites[1].alleles[0] == "C");
        REQUIRE(sites[1].alleles[1] == "A");
        REQUIRE(sites[1].alleles[2] == "G");
        REQUIRE(sites[1].alleles[3] == "T");
        REQUIRE(sites[1].unification[make_pair(range(0,1001,1002),string("C"))] == 0);
        REQUIRE(sites[1].unification[make_pair(range(0,1001,1002),string("A"))] == 1);
        REQUIRE(sites[1].unification[make_pair(range(0,1001,1002),string("G"))] == 2);
        REQUIRE(sites[1].unification[make_pair(range(0,1001,1002),string("T"))] == 3);
        REQUIRE(sites[1].observation_count.size() == 4);
        REQUIRE(sites[1].observation_count[0] == 2);
        REQUIRE(sites[1].observation_count[1] == 6);
        REQUIRE(sites[1].observation_count[2] == 2);
        REQUIRE(sites[1].observation_count[3] == 2);

        REQUIRE(sites[2].pos == range(0,1010,1012));
        REQUIRE(sites[2].alleles.size() == 2);
        REQUIRE(sites[2].alleles[0] == "CC");
        REQUIRE(sites[2].alleles[1] == "AG");
        REQUIRE(sites[2].unification[make_pair(range(0,1010,1012),string("CC"))] == 0);
        REQUIRE(sites[2].unification[make_pair(range(0,1010,1012),string("AG"))] == 1);
        REQUIRE(sites[2].observation_count.size() == 2);
        REQUIRE(sites[2].observation_count[0] == 3);
        REQUIRE(sites[2].observation_count[1] == 3);

        REQUIRE(sites[3].pos == range(0,1100,1101));
        REQUIRE(sites[3].alleles.size() == 2);
        REQUIRE(sites[3].alleles[0] == "C");
        REQUIRE(sites[3].alleles[1] == "A");
        REQUIRE(sites[3].unification[make_pair(range(0,1100,1101),string("C"))] == 0);
        REQUIRE(sites[3].unification[make_pair(range(0,1100,1101),string("A"))] == 1);
        REQUIRE(sites[3].observation_count.size() == 2);
        REQUIRE(sites[3].observation_count[0] == 3);
        REQUIRE(sites[3].observation_count[1] == 3);

        REQUIRE(sites[4].pos == range(0,1102,1103));
        REQUIRE(sites[4].alleles.size() == 2);
        REQUIRE(sites[4].alleles[0] == "C");
        REQUIRE(sites[4].alleles[1] == "G");
        REQUIRE(sites[4].unification[make_pair(range(0,1102,1103),string("C"))] == 0);
        REQUIRE(sites[4].unification[make_pair(range(0,1102,1103),string("G"))] == 1);
        REQUIRE(sites[4].observation_count.size() == 2);
        REQUIRE(sites[4].observation_count[0] == 3);
        REQUIRE(sites[4].observation_count[1] == 3);

        REQUIRE(sites[5].pos == range(0,1200,1201));
        REQUIRE(sites[5].alleles.size() == 2);
        REQUIRE(sites[5].alleles[0] == "C");
        REQUIRE(sites[5].alleles[1] == "A");
        REQUIRE(sites[5].unification[make_pair(range(0,1200,1201),string("C"))] == 0);
        REQUIRE(sites[5].unification[make_pair(range(0,1200,1201),string("A"))] == 1);
        REQUIRE(sites[5].observation_count.size() == 2);
        REQUIRE(sites[5].observation_count[0] == 3);
        REQUIRE(sites[5].observation_count[1] == 3);

        // An allele from trio2 that would collapse sites[3] and sites[4] was
        // pruned.

        /*
        for (const auto& site : sites) {
            cout << site.pos.str(contigs);
            for (unsigned i = 0; i < site.alleles.size(); i++) {
                cout << ' ' << site.alleles[i] << '[' << site.observation_count[i] << ']';
            }
            for (const auto& u : site.unification) {
                UNPAIR(u, p, i)
                UNPAIR(p, pos, suballele)
                cout << ' ' << suballele << '@' << pos << "->" << site.alleles[i];
            }
            cout << endl;
        }
        */
    }
}
TEST_CASE("genotyper with ref changes") {
    unique_ptr<VCFData> data;
    Status s = VCFData::Open({"joint_A.gvcf", "joint_B.gvcf", "joint_C.gvcf", "joint_D.gvcf"}, data);
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(*data, *data, svc);
    REQUIRE(s.ok());

    discovered_alleles als;

    const string tfn("/tmp/GLnexus_unit_tests.bcf");

    s = svc->discover_alleles("<ALL>", range(0, 0, 1000000), als);
    REQUIRE(s.ok());

    vector<unified_site> sites;
    s = unified_sites(als, sites);
    REQUIRE(s.ok());

    consolidated_loss losses;
    s = svc->genotype_sites(genotyper_config(), string("<ALL>"), sites, tfn, losses);
    REQUIRE(s.ok());

    unique_ptr<vcfFile, void(*)(vcfFile*)> vcf(bcf_open(tfn.c_str(), "r"), [](vcfFile* f) { bcf_close(f); });
    REQUIRE(bool(vcf));

    shared_ptr<bcf_hdr_t> hdr(bcf_hdr_read(vcf.get()), &bcf_hdr_destroy);
    REQUIRE(bool(hdr));

    int nsamples = bcf_hdr_nsamples(hdr);
    REQUIRE(nsamples == 4);
    REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 0)) == "JOINT_A");
    REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 1)) == "JOINT_B");
    REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 2)) == "JOINT_C");
    REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 3)) == "JOINT_D");

    unique_ptr<bcf1_t,void(*)(bcf1_t*)> record(bcf_init(), &bcf_destroy);
    REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
    REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

    REQUIRE(record->n_allele == 2);
    REQUIRE(string(record->d.allele[0]) == "ATG");
    REQUIRE(string(record->d.allele[1]) == "ATGTG");

    int *gt = nullptr, gtsz = 0;
    int nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
    REQUIRE(nGT == 8);

    // Because the REF has changed in JOINT_A, both ALT alleles
    // should be lost instead of being interpreted as is
    REQUIRE(bcf_gt_is_missing(gt[0]));
    REQUIRE(bcf_gt_is_missing(gt[1]));
    REQUIRE(bcf_gt_allele(gt[2]) == 0);
    REQUIRE(bcf_gt_allele(gt[3]) == 1);
    REQUIRE(bcf_gt_allele(gt[4]) == 0);
    REQUIRE(bcf_gt_allele(gt[5]) == 1);
    REQUIRE(bcf_gt_allele(gt[6]) == 0);
    REQUIRE(bcf_gt_allele(gt[7]) == 1);
}

TEST_CASE("genotyper placeholder") {
    unique_ptr<VCFData> data;
    Status s = VCFData::Open({"discover_alleles_trio1.vcf", "discover_alleles_trio2.vcf"}, data);
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(*data, *data, svc);
    REQUIRE(s.ok());

    discovered_alleles als;

    const string tfn("/tmp/GLnexus_unit_tests.bcf");

    SECTION("discover_alleles_trio1") {
        s = svc->discover_alleles("discover_alleles_trio1", range(0, 0, 1000000), als);
        REQUIRE(s.ok());

        vector<unified_site> sites;
        s = unified_sites(als, sites);
        REQUIRE(s.ok());

        consolidated_loss losses;
        s = svc->genotype_sites(genotyper_config(), string("discover_alleles_trio1"), sites, tfn, losses);
        REQUIRE(s.ok());

        unique_ptr<vcfFile, void(*)(vcfFile*)> vcf(bcf_open(tfn.c_str(), "r"), [](vcfFile* f) { bcf_close(f); });
        REQUIRE(bool(vcf));

        shared_ptr<bcf_hdr_t> hdr(bcf_hdr_read(vcf.get()), &bcf_hdr_destroy);
        REQUIRE(bool(hdr));

        int nsamples = bcf_hdr_nsamples(hdr);
        REQUIRE(nsamples == 3);
        REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 0)) == "trio1.ch");
        REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 1)) == "trio1.fa");
        REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 2)) == "trio1.mo");

        unique_ptr<bcf1_t,void(*)(bcf1_t*)> record(bcf_init(), &bcf_destroy);
        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "A");
        REQUIRE(string(record->d.allele[1]) == "G");

        int *gt = nullptr, gtsz = 0;
        int nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 6);
        REQUIRE(bcf_gt_allele(gt[0]) == 1);
        REQUIRE(bcf_gt_allele(gt[1]) == 1);
        REQUIRE(bcf_gt_allele(gt[2]) == 0);
        REQUIRE(bcf_gt_allele(gt[3]) == 1);
        REQUIRE(bcf_gt_allele(gt[4]) == 0);
        REQUIRE(bcf_gt_allele(gt[5]) == 1);

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 3);
        REQUIRE(string(record->d.allele[0]) == "C");
        REQUIRE(string(record->d.allele[1]) == "G");
        REQUIRE(string(record->d.allele[2]) == "T");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 6);
        REQUIRE(bcf_gt_allele(gt[0]) == 2);
        REQUIRE(bcf_gt_allele(gt[1]) == 2);
        REQUIRE(bcf_gt_allele(gt[2]) == 0);
        REQUIRE(bcf_gt_allele(gt[3]) == 0);
        REQUIRE(bcf_gt_allele(gt[4]) == 1);
        REQUIRE(bcf_gt_allele(gt[5]) == 1);

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "CC");
        REQUIRE(string(record->d.allele[1]) == "AG");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 6);
        REQUIRE(bcf_gt_allele(gt[0]) == 0);
        REQUIRE(bcf_gt_allele(gt[1]) == 1);
        REQUIRE(bcf_gt_allele(gt[2]) == 0);
        REQUIRE(bcf_gt_allele(gt[3]) == 1);
        REQUIRE(bcf_gt_allele(gt[4]) == 0);
        REQUIRE(bcf_gt_allele(gt[5]) == 1);

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "C");
        REQUIRE(string(record->d.allele[1]) == "A");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 6);
        REQUIRE(bcf_gt_allele(gt[0]) == 0);
        REQUIRE(bcf_gt_allele(gt[1]) == 1);
        REQUIRE(bcf_gt_allele(gt[2]) == 0);
        REQUIRE(bcf_gt_allele(gt[3]) == 1);
        REQUIRE(bcf_gt_allele(gt[4]) == 0);
        REQUIRE(bcf_gt_allele(gt[5]) == 1);

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "C");
        REQUIRE(string(record->d.allele[1]) == "G");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 6);
        REQUIRE(bcf_gt_allele(gt[0]) == 0);
        REQUIRE(bcf_gt_allele(gt[1]) == 1);
        REQUIRE(bcf_gt_allele(gt[2]) == 0);
        REQUIRE(bcf_gt_allele(gt[3]) == 1);
        REQUIRE(bcf_gt_allele(gt[4]) == 0);
        REQUIRE(bcf_gt_allele(gt[5]) == 1);

        // validate loss information
        auto loss1 = losses.find("trio1.ch");
        REQUIRE(loss1 != losses.end());
        loss_stats loss_ch = loss1->second;
        REQUIRE(loss_ch.n_no_calls_total == 0);
        REQUIRE(loss_ch.n_bp_lost == 0);

        auto loss2 = losses.find("trio1.fa");
        REQUIRE(loss2 != losses.end());
        loss_stats loss_fa = loss2->second;
        REQUIRE(loss_fa.n_no_calls_total == 0);
        REQUIRE(loss_fa.n_bp_lost == 0);

        auto loss3 = losses.find("trio1.mo");
        REQUIRE(loss3 != losses.end());
        loss_stats loss_mo = loss3->second;
        REQUIRE(loss_mo.n_no_calls_total == 0);
        REQUIRE(loss_mo.n_bp_lost == 0);

        free(gt);
    }

    SECTION("2 trios") {
        s = svc->discover_alleles("<ALL>", range(0, 0, 1000000), als);
        REQUIRE(s.ok());

        vector<unified_site> sites;
        s = unified_sites(als, sites);
        REQUIRE(s.ok());

        consolidated_loss losses;
        s = svc->genotype_sites(genotyper_config(), string("<ALL>"), sites, tfn, losses);
        REQUIRE(s.ok());

        unique_ptr<vcfFile, void(*)(vcfFile*)> vcf(bcf_open(tfn.c_str(), "r"), [](vcfFile* f) { bcf_close(f); });
        REQUIRE(bool(vcf));

        shared_ptr<bcf_hdr_t> hdr(bcf_hdr_read(vcf.get()), &bcf_hdr_destroy);
        REQUIRE(bool(hdr));

        int nsamples = bcf_hdr_nsamples(hdr);
        REQUIRE(nsamples == 6);
        REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 0)) == "trio1.ch");
        REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 1)) == "trio1.fa");
        REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 2)) == "trio1.mo");
        REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 3)) == "trio2.ch");
        REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 4)) == "trio2.fa");
        REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 5)) == "trio2.mo");

        unique_ptr<bcf1_t,void(*)(bcf1_t*)> record(bcf_init(), &bcf_destroy);
        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(range(record.get()) == range(0, 1000, 1001));
        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "A");
        REQUIRE(string(record->d.allele[1]) == "G");

        int *gt = nullptr, gtsz = 0;
        int nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 12);
        REQUIRE(bcf_gt_allele(gt[0]) == 1);
        REQUIRE(bcf_gt_allele(gt[1]) == 1);
        REQUIRE(bcf_gt_allele(gt[2]) == 0);
        REQUIRE(bcf_gt_allele(gt[3]) == 1);
        REQUIRE(bcf_gt_allele(gt[4]) == 0);
        REQUIRE(bcf_gt_allele(gt[5]) == 1);
        REQUIRE(bcf_gt_allele(gt[6]) == 0);
        REQUIRE(bcf_gt_allele(gt[7]) == 0);
        REQUIRE(bcf_gt_allele(gt[8]) == 0);
        REQUIRE(bcf_gt_allele(gt[9]) == 1);
        REQUIRE(bcf_gt_allele(gt[10]) == 0);
        REQUIRE(bcf_gt_allele(gt[11]) == 1);

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(range(record.get()) == range(0, 1001, 1002));
        REQUIRE(record->n_allele == 4);
        REQUIRE(string(record->d.allele[0]) == "C");
        REQUIRE(string(record->d.allele[1]) == "A");
        REQUIRE(string(record->d.allele[2]) == "G");
        REQUIRE(string(record->d.allele[3]) == "T");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 12);
        REQUIRE(bcf_gt_allele(gt[0]) == 3);
        REQUIRE(bcf_gt_allele(gt[1]) == 3);
        REQUIRE(bcf_gt_allele(gt[2]) == 0);
        REQUIRE(bcf_gt_allele(gt[3]) == 0);
        REQUIRE(bcf_gt_allele(gt[4]) == 2);
        REQUIRE(bcf_gt_allele(gt[5]) == 2);
        REQUIRE(bcf_gt_allele(gt[6]) == 1);
        REQUIRE(bcf_gt_allele(gt[7]) == 1);
        REQUIRE(bcf_gt_allele(gt[8]) == 1);
        REQUIRE(bcf_gt_allele(gt[9]) == 1);
        REQUIRE(bcf_gt_allele(gt[10]) == 1);
        REQUIRE(bcf_gt_allele(gt[11]) == 1);

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "CC");
        REQUIRE(string(record->d.allele[1]) == "AG");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 12);
        REQUIRE(bcf_gt_allele(gt[0]) == 0);
        REQUIRE(bcf_gt_allele(gt[1]) == 1);
        REQUIRE(bcf_gt_allele(gt[2]) == 0);
        REQUIRE(bcf_gt_allele(gt[3]) == 1);
        REQUIRE(bcf_gt_allele(gt[4]) == 0);
        REQUIRE(bcf_gt_allele(gt[5]) == 1);
        REQUIRE(bcf_gt_allele(gt[6]) == 0);
        REQUIRE(bcf_gt_allele(gt[7]) == 0);
        REQUIRE(bcf_gt_allele(gt[8]) == 0);
        REQUIRE(bcf_gt_is_missing(gt[9]));
        REQUIRE(bcf_gt_allele(gt[10]) == 0);
        REQUIRE(bcf_gt_is_missing(gt[11]));

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "C");
        REQUIRE(string(record->d.allele[1]) == "A");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 12);
        REQUIRE(bcf_gt_allele(gt[0]) == 0);
        REQUIRE(bcf_gt_allele(gt[1]) == 1);
        REQUIRE(bcf_gt_allele(gt[2]) == 0);
        REQUIRE(bcf_gt_allele(gt[3]) == 1);
        REQUIRE(bcf_gt_allele(gt[4]) == 0);
        REQUIRE(bcf_gt_allele(gt[5]) == 1);
        REQUIRE(bcf_gt_allele(gt[6]) == 0);
        REQUIRE(bcf_gt_allele(gt[7]) == 0);
        REQUIRE(bcf_gt_allele(gt[8]) == 0);
        REQUIRE(bcf_gt_is_missing(gt[9]));
        REQUIRE(bcf_gt_allele(gt[10]) == 0);
        REQUIRE(bcf_gt_is_missing(gt[11]));

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "C");
        REQUIRE(string(record->d.allele[1]) == "G");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 12);
        REQUIRE(bcf_gt_allele(gt[0]) == 0);
        REQUIRE(bcf_gt_allele(gt[1]) == 1);
        REQUIRE(bcf_gt_allele(gt[2]) == 0);
        REQUIRE(bcf_gt_allele(gt[3]) == 1);
        REQUIRE(bcf_gt_allele(gt[4]) == 0);
        REQUIRE(bcf_gt_allele(gt[5]) == 1);
        REQUIRE(bcf_gt_allele(gt[6]) == 0);
        REQUIRE(bcf_gt_allele(gt[7]) == 0);
        REQUIRE(bcf_gt_allele(gt[8]) == 0);
        REQUIRE(bcf_gt_is_missing(gt[9]));
        REQUIRE(bcf_gt_allele(gt[10]) == 0);
        REQUIRE(bcf_gt_is_missing(gt[11]));

        char* rnc = nullptr;
        int rncsz = 0;
        int nRNC = bcf_get_format_char(hdr.get(), record.get(), "RNC", &rnc, &rncsz);
        REQUIRE(nRNC == 12);
        REQUIRE(rncsz == 12);
        REQUIRE(string(rnc,12) == ".........L.L");

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "C");
        REQUIRE(string(record->d.allele[1]) == "A");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 12);
        REQUIRE(bcf_gt_allele(gt[0]) == 0);
        REQUIRE(bcf_gt_allele(gt[1]) == 1);
        REQUIRE(bcf_gt_allele(gt[2]) == 0);
        REQUIRE(bcf_gt_allele(gt[3]) == 1);
        REQUIRE(bcf_gt_allele(gt[4]) == 0);
        REQUIRE(bcf_gt_allele(gt[5]) == 1);
        REQUIRE(bcf_gt_is_missing(gt[6]));
        REQUIRE(bcf_gt_is_missing(gt[7]));
        REQUIRE(bcf_gt_is_missing(gt[8]));
        REQUIRE(bcf_gt_is_missing(gt[9]));
        REQUIRE(bcf_gt_is_missing(gt[10]));
        REQUIRE(bcf_gt_is_missing(gt[11]));

        nRNC = bcf_get_format_char(hdr.get(), record.get(), "RNC", &rnc, &rncsz);
        REQUIRE(nRNC == 12);
        REQUIRE(rncsz == 12);
        REQUIRE(string(rnc,12) == "......MMMMMM");

        free(gt);
        free(rnc);

        // validate loss information
        auto loss1 = losses.find("trio1.ch");
        REQUIRE(loss1 != losses.end());
        loss_stats loss_ch = loss1->second;
        REQUIRE(loss_ch.n_no_calls_total == 0);
        REQUIRE(loss_ch.n_bp_lost == 0);

        auto loss2 = losses.find("trio1.fa");
        REQUIRE(loss2 != losses.end());
        loss_stats loss_fa = loss2->second;
        REQUIRE(loss_fa.n_no_calls_total == 0);
        REQUIRE(loss_fa.n_bp_lost == 0);

        auto loss3 = losses.find("trio1.mo");
        REQUIRE(loss3 != losses.end());
        loss_stats loss_mo = loss3->second;
        REQUIRE(loss_mo.n_no_calls_total == 0);
        REQUIRE(loss_mo.n_bp_lost == 0);

        auto loss4 = losses.find("trio2.ch");
        REQUIRE(loss4 != losses.end());
        loss_stats loss_ch2 = loss4->second;
        REQUIRE(loss_ch2.n_no_calls_total == 2);
        REQUIRE(loss_ch2.n_bp_lost == 0);

        int expected_loss_bp = sites[2].pos.size() + sites[3].pos.size() + sites[4].pos.size();

        auto loss5 = losses.find("trio2.fa");
        REQUIRE(loss5 != losses.end());
        loss_stats loss_fa2 = loss5->second;
        REQUIRE(loss_fa2.n_no_calls_total == 5);
        REQUIRE(loss_fa2.n_bp_lost == expected_loss_bp);

        auto loss6 = losses.find("trio2.mo");
        REQUIRE(loss6 != losses.end());
        loss_stats loss_mo2 = loss6->second;
        REQUIRE(loss_mo2.n_no_calls_total == 5);
        REQUIRE(loss_mo2.n_bp_lost == expected_loss_bp);
    }

    // alleles discovered from all samples, genotyped only in trio2
    SECTION("trio2 with all alleles") {
        s = svc->discover_alleles("<ALL>", range(0, 0, 1000000), als);
        REQUIRE(s.ok());

        vector<unified_site> sites;
        s = unified_sites(als, sites);
        REQUIRE(s.ok());

        consolidated_loss losses;
        s = svc->genotype_sites(genotyper_config(), string("discover_alleles_trio2"), sites, tfn, losses);
        REQUIRE(s.ok());

        unique_ptr<vcfFile, void(*)(vcfFile*)> vcf(bcf_open(tfn.c_str(), "r"), [](vcfFile* f) { bcf_close(f); });
        REQUIRE(bool(vcf));

        shared_ptr<bcf_hdr_t> hdr(bcf_hdr_read(vcf.get()), &bcf_hdr_destroy);
        REQUIRE(bool(hdr));

        int nsamples = bcf_hdr_nsamples(hdr);
        REQUIRE(nsamples == 3);
        REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 0)) == "trio2.ch");
        REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 1)) == "trio2.fa");
        REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 2)) == "trio2.mo");

        unique_ptr<bcf1_t,void(*)(bcf1_t*)> record(bcf_init(), &bcf_destroy);
        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "A");
        REQUIRE(string(record->d.allele[1]) == "G");

        int *gt = nullptr, gtsz = 0;
        int nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 6);
        REQUIRE(bcf_gt_allele(gt[0]) == 0);
        REQUIRE(bcf_gt_allele(gt[1]) == 0);
        REQUIRE(bcf_gt_allele(gt[2]) == 0);
        REQUIRE(bcf_gt_allele(gt[3]) == 1);
        REQUIRE(bcf_gt_allele(gt[4]) == 0);
        REQUIRE(bcf_gt_allele(gt[5]) == 1);

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 4);
        REQUIRE(string(record->d.allele[0]) == "C");
        REQUIRE(string(record->d.allele[1]) == "A");
        REQUIRE(string(record->d.allele[2]) == "G");
        REQUIRE(string(record->d.allele[3]) == "T");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 6);
        REQUIRE(bcf_gt_allele(gt[0]) == 1);
        REQUIRE(bcf_gt_allele(gt[1]) == 1);
        REQUIRE(bcf_gt_allele(gt[2]) == 1);
        REQUIRE(bcf_gt_allele(gt[3]) == 1);
        REQUIRE(bcf_gt_allele(gt[4]) == 1);
        REQUIRE(bcf_gt_allele(gt[5]) == 1);

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "CC");
        REQUIRE(string(record->d.allele[1]) == "AG");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 6);
        REQUIRE(bcf_gt_allele(gt[0]) == 0);
        REQUIRE(bcf_gt_allele(gt[1]) == 0);
        REQUIRE(bcf_gt_allele(gt[2]) == 0);
        REQUIRE(bcf_gt_is_missing(gt[3]));
        REQUIRE(bcf_gt_allele(gt[4]) == 0);
        REQUIRE(bcf_gt_is_missing(gt[5]));

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "C");
        REQUIRE(string(record->d.allele[1]) == "A");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 6);
        REQUIRE(bcf_gt_allele(gt[0]) == 0);
        REQUIRE(bcf_gt_allele(gt[1]) == 0);
        REQUIRE(bcf_gt_allele(gt[2]) == 0);
        REQUIRE(bcf_gt_is_missing(gt[3]));
        REQUIRE(bcf_gt_allele(gt[4]) == 0);
        REQUIRE(bcf_gt_is_missing(gt[5]));

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "C");
        REQUIRE(string(record->d.allele[1]) == "G");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 6);
        REQUIRE(bcf_gt_allele(gt[0]) == 0);
        REQUIRE(bcf_gt_allele(gt[1]) == 0);
        REQUIRE(bcf_gt_allele(gt[2]) == 0);
        REQUIRE(bcf_gt_is_missing(gt[3]));
        REQUIRE(bcf_gt_allele(gt[4]) == 0);
        REQUIRE(bcf_gt_is_missing(gt[5]));

        auto loss4 = losses.find("trio2.ch");
        REQUIRE(loss4 != losses.end());
        loss_stats loss_ch2 = loss4->second;
        REQUIRE(loss_ch2.n_no_calls_total == 2);
        REQUIRE(loss_ch2.n_bp_lost == 0);

        int expected_loss_bp = sites[2].pos.size() + sites[3].pos.size() + sites[4].pos.size();

        auto loss5 = losses.find("trio2.fa");
        REQUIRE(loss5 != losses.end());
        loss_stats loss_fa2 = loss5->second;
        REQUIRE(loss_fa2.n_no_calls_total == 5);
        REQUIRE(loss_fa2.n_bp_lost == expected_loss_bp);

        auto loss6 = losses.find("trio2.mo");
        REQUIRE(loss6 != losses.end());
        loss_stats loss_mo2 = loss6->second;
        REQUIRE(loss_mo2.n_no_calls_total == 5);
        REQUIRE(loss_mo2.n_bp_lost == expected_loss_bp);

        free(gt);
    }

    SECTION("unification with multiple contigs") {
        discovered_alleles als0, als1;

        s = svc->discover_alleles("<ALL>", range(0, 0, 1000000), als0);
        REQUIRE(s.ok());

        s = svc->discover_alleles("<ALL>", range(1, 0, 1000000), als1);
        REQUIRE(s.ok());

        als.clear();
        REQUIRE(merge_discovered_alleles(als0, als).ok());
        REQUIRE(merge_discovered_alleles(als1, als).ok());

        vector<unified_site> sites;
        s = unified_sites(als, sites);
        REQUIRE(s.ok());

        REQUIRE(is_sorted(sites.begin(), sites.end()));
        REQUIRE(sites[0].pos.rid == 0);
        REQUIRE(sites[sites.size()-1].pos.rid == 1);
    }

    SECTION("simulate I/O errors") {
        s = Service::Start(*data, *data, svc);
        REQUIRE(s.ok());
        s = svc->discover_alleles("<ALL>", range(0, 0, 1000000), als);
        REQUIRE(s.ok());
        vector<unified_site> sites;
        s = unified_sites(als, sites);
        REQUIRE(s.ok());

        unique_ptr<SimFailBCFData> faildata;
        bool worked = false;

        for (size_t fail_every = 1; fail_every < 100; fail_every++) {
            s = SimFailBCFData::Open(*data, fail_every, faildata);
            REQUIRE(s.ok());

            s = Service::Start(*data, *faildata, svc);
            REQUIRE(s.ok());

            consolidated_loss losses;
            s = svc->genotype_sites(genotyper_config(), string("<ALL>"), sites, tfn, losses);
            if (faildata->failed_once()) {
                worked = true;
                REQUIRE(s == StatusCode::IO_ERROR);
                REQUIRE(s.str() == "IOError: SIM");
            } else {
                REQUIRE(s.ok());
            }
        }

        REQUIRE(worked);
    }
}

TEST_CASE("gVCF genotyper") {
    unique_ptr<VCFData> data;
    Status s = VCFData::Open({"NA12878D_HiSeqX.21.10009462-10009469.gvcf"}, data);
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(*data, *data, svc);
    REQUIRE(s.ok());

    discovered_alleles als;

    const string tfn("/tmp/GLnexus_unit_tests.bcf");

    SECTION("no depth requirement") {
        vector<unified_site> sites;

        unified_site us(range(0,10009461,10009462));
        us.alleles.push_back("T");
        us.alleles.push_back("A");
        us.unification[make_pair(range(0,10009461,10009462),"T")] = 0;
        us.unification[make_pair(range(0,10009461,10009462),"A")] = 1;
        us.observation_count = { 1, 1 };
        sites.push_back(us);

        us.pos = range(0,10009462,10009463);
        us.alleles.clear();
        us.alleles.push_back("C");
        us.alleles.push_back("G");
        us.unification.clear();
        us.unification[make_pair(range(0,10009462,10009463),"C")] = 0;
        us.unification[make_pair(range(0,10009462,10009463),"G")] = 1;
        us.observation_count = { 1, 1 };
        sites.push_back(us);

        // this site spans two gVCF records and will not be called by the current algorithm.
        us.pos = range(0,10009465,10009467);
        us.alleles.clear();
        us.alleles.push_back("AA");
        us.alleles.push_back("GT");
        us.unification.clear();
        us.unification[make_pair(range(0,10009465,10009467),"AA")] = 0;
        us.unification[make_pair(range(0,10009465,10009467),"GT")] = 1;
        us.observation_count = { 1, 1 };
        sites.push_back(us);

        consolidated_loss losses;
        s = svc->genotype_sites(genotyper_config(), string("NA12878D_HiSeqX.21.10009462-10009469"), sites, tfn, losses);
        REQUIRE(s.ok());

        unique_ptr<vcfFile, void(*)(vcfFile*)> vcf(bcf_open(tfn.c_str(), "r"), [](vcfFile* f) { bcf_close(f); });
        REQUIRE(bool(vcf));
        shared_ptr<bcf_hdr_t> hdr(bcf_hdr_read(vcf.get()), &bcf_hdr_destroy);
        REQUIRE(bool(hdr));

        int nsamples = bcf_hdr_nsamples(hdr);
        REQUIRE(nsamples == 1);
        REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 0)) == "NA12878");

        unique_ptr<bcf1_t,void(*)(bcf1_t*)> record(bcf_init(), &bcf_destroy);
        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "T");
        REQUIRE(string(record->d.allele[1]) == "A");

        // With no depth requirement, we should just get the hard-called genotypes back.

        int *gt = nullptr, gtsz = 0;
        int nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 2);
        REQUIRE(bcf_gt_allele(gt[0]) == 0);
        REQUIRE(bcf_gt_allele(gt[1]) == 0);

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "C");
        REQUIRE(string(record->d.allele[1]) == "G");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 2);
        REQUIRE(bcf_gt_allele(gt[0]) == 0);
        REQUIRE(bcf_gt_allele(gt[1]) == 0);

        // this site spans two gVCF records and will not be called by the current algorithm.
        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "AA");
        REQUIRE(string(record->d.allele[1]) == "GT");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 2);
        REQUIRE(bcf_gt_is_missing(gt[0]));
        REQUIRE(bcf_gt_is_missing(gt[1]));

        auto loss = losses.find("NA12878");
        REQUIRE(loss != losses.end());
        loss_stats loss_na = loss->second;
        REQUIRE(loss_na.n_no_calls_total == 2);
        REQUIRE(loss_na.n_bp_lost == 4);
        REQUIRE(loss_na.n_calls_lost == 4);
        REQUIRE(loss_na.n_gvcf_calls_lost == 4);
        REQUIRE(loss_na.n_gvcf_bp_lost == 4);
    }

    SECTION("require depth > 12") {
        vector<unified_site> sites;

        unified_site us(range(0,10009463,10009465));
        us.alleles.push_back("TA");
        us.alleles.push_back("T");
        us.unification[make_pair(range(0,10009463,10009465),"TA")] = 0;
        us.unification[make_pair(range(0,10009463,10009465),"T")] = 1;
        us.observation_count = { 1, 1 };
        sites.push_back(us);

        us.pos = range(0,10009465,10009466);
        us.alleles.clear();
        us.alleles.push_back("A");
        us.alleles.push_back("G");
        us.unification.clear();
        us.unification[make_pair(range(0,10009465,10009466),"A")] = 0;
        us.unification[make_pair(range(0,10009465,10009466),"G")] = 1;
        us.observation_count = { 1, 1 };
        sites.push_back(us);

        us.pos = range(0,10009466,10009467);
        us.alleles.clear();
        us.alleles.push_back("A");
        us.alleles.push_back("C");
        us.unification.clear();
        us.unification[make_pair(range(0,10009466,10009467),"A")] = 0;
        us.unification[make_pair(range(0,10009466,10009467),"C")] = 1;
        us.observation_count = { 1, 1 };
        sites.push_back(us);

        genotyper_config cfg;
        cfg.required_dp = 13;

        consolidated_loss losses;
        s = svc->genotype_sites(cfg, string("NA12878D_HiSeqX.21.10009462-10009469"), sites, tfn, losses);
        REQUIRE(s.ok());

        unique_ptr<vcfFile, void(*)(vcfFile*)> vcf(bcf_open(tfn.c_str(), "r"), [](vcfFile* f) { bcf_close(f); });
        REQUIRE(bool(vcf));
        shared_ptr<bcf_hdr_t> hdr(bcf_hdr_read(vcf.get()), &bcf_hdr_destroy);
        REQUIRE(bool(hdr));

        int nsamples = bcf_hdr_nsamples(hdr);
        REQUIRE(nsamples == 1);
        REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 0)) == "NA12878");

        // Some calls will be missing due to the depth requirement.
        unique_ptr<bcf1_t,void(*)(bcf1_t*)> record(bcf_init(), &bcf_destroy);
        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "TA");
        REQUIRE(string(record->d.allele[1]) == "T");

        int *gt = nullptr, gtsz = 0;
        int nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 2);
        REQUIRE(bcf_gt_is_missing(gt[0]));
        REQUIRE(bcf_gt_is_missing(gt[1]));


        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "A");
        REQUIRE(string(record->d.allele[1]) == "G");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 2);
        REQUIRE(bcf_gt_allele(gt[0]) == 0);
        REQUIRE(bcf_gt_allele(gt[1]) == 0);


        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "A");
        REQUIRE(string(record->d.allele[1]) == "C");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 2);
        REQUIRE(bcf_gt_is_missing(gt[0]));
        REQUIRE(bcf_gt_is_missing(gt[1]));

        auto loss = losses.find("NA12878");
        REQUIRE(loss != losses.end());
        loss_stats loss_d12 = loss->second;
        REQUIRE(loss_d12.n_no_calls_total == 4);
        REQUIRE(loss_d12.n_bp_lost == 6);
        REQUIRE(loss_d12.n_calls_lost == 4);
        REQUIRE(loss_d12.n_gvcf_bp_lost == 2);
        REQUIRE(loss_d12.n_gvcf_calls_lost == 2);
    }


    SECTION("require depth > 9") {
        vector<unified_site> sites;

        unified_site us(range(0,10009463,10009465));
        us.alleles.push_back("TA");
        us.alleles.push_back("T");
        us.unification[make_pair(range(0,10009463,10009465),"TA")] = 0;
        us.unification[make_pair(range(0,10009463,10009465),"T")] = 1;
        us.observation_count = { 1, 1 };
        sites.push_back(us);

        genotyper_config cfg;
        cfg.required_dp = 9;
        consolidated_loss losses;
        s = svc->genotype_sites(cfg, string("NA12878D_HiSeqX.21.10009462-10009469"), sites, tfn, losses);
        REQUIRE(s.ok());

        unique_ptr<vcfFile, void(*)(vcfFile*)> vcf(bcf_open(tfn.c_str(), "r"), [](vcfFile* f) { bcf_close(f); });
        REQUIRE(bool(vcf));
        shared_ptr<bcf_hdr_t> hdr(bcf_hdr_read(vcf.get()), &bcf_hdr_destroy);
        REQUIRE(bool(hdr));

        int nsamples = bcf_hdr_nsamples(hdr);
        REQUIRE(nsamples == 1);
        REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 0)) == "NA12878");

        // One of the allele calls will be missing due to the depth requirement.
        unique_ptr<bcf1_t,void(*)(bcf1_t*)> record(bcf_init(), &bcf_destroy);
        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "TA");
        REQUIRE(string(record->d.allele[1]) == "T");

        int *gt = nullptr, gtsz = 0;
        int nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 2);
        REQUIRE(bcf_gt_allele(gt[0]) == 0);
        REQUIRE(bcf_gt_is_missing(gt[1]));
        // TODO consider normalizing half-calls so that the missing allele is always first (or second)

        auto loss = losses.find("NA12878");
        REQUIRE(loss != losses.end());
        loss_stats loss_d9 = loss->second;
        REQUIRE(loss_d9.n_no_calls_total == 1);
        REQUIRE(loss_d9.n_bp_lost == 2);
        REQUIRE(loss_d9.n_calls_lost == 1);
        REQUIRE(loss_d9.n_gvcf_calls_lost == 0);
        REQUIRE(loss_d9.n_gvcf_bp_lost == 0);
    }
}


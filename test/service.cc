#include <iostream>
#include <algorithm>
#include <vcf.h>
#include "service.h"
#include "alleles.h"
#include "catch.hpp"
using namespace std;
using namespace GLnexus;

// serves data from VCF files in the test/data directory
// x.vcf is loaded as data set "x" with one sample set "x"
// additionally, the sample set "<ALL>" designates all samples across the VCFs.
class VCFData : public Data {
    struct vcf_data_t {
        shared_ptr<bcf_hdr_t> header;
        shared_ptr<const vector<string> > samples;
        vector<shared_ptr<bcf1_t> > records;
    };
    map<string,vcf_data_t> datasets_;
    map<string,string> sample_datasets_;

    VCFData() {}

    static Status load_vcf(const string& name, vcf_data_t& ans) {
        string path = "test/data/" + name + ".vcf";
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
            datasets[nm] = d;
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
        } else if (sampleset == "<ALL>") {
            auto s = make_shared<set<string>>();
            for (const auto& ds : datasets_) {
                s->insert(ds.second.samples->begin(), ds.second.samples->end());
            }
            ans = move(s);
            return Status::OK();
        }
        return Status::NotFound("unknown sample set", sampleset);
    }

    Status sample_dataset(const string& sample, string& ans) const override {
        auto p = sample_datasets_.find(sample);
        if (p == sample_datasets_.end()) {
            return Status::NotFound("unknown sample", sample);
        }
        ans = p->second;
        return Status::OK();
    }

    Status dataset_bcf_header(const string& dataset, shared_ptr<const bcf_hdr_t>& hdr) const override {
        auto p = datasets_.find(dataset);
        if (p == datasets_.end()) {
            return Status::NotFound("unknown data set", dataset);
        }
        hdr = p->second.header;
        return Status::OK();
    }

    Status dataset_bcf(const string& dataset, const range& pos, shared_ptr<const bcf_hdr_t>& hdr, vector<shared_ptr<bcf1_t> >& records) const override {
        auto p = datasets_.find(dataset);
        if (p == datasets_.end()) {
            return Status::NotFound("unknown data set", dataset);
        }
        hdr = p->second.header;
        records.clear();
        for (const auto& bcf : p->second.records) {
            if (range(bcf).overlaps(pos)) {
                records.push_back(bcf);
            }
        }
        return Status::OK();
    }
};

TEST_CASE("service::discover_alleles") {
    unique_ptr<VCFData> data;
    Status s = VCFData::Open({"discover_alleles_trio1", "discover_alleles_trio2"}, data);
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(data.get(), svc);
    REQUIRE(s.ok());

    discovered_alleles als;

    SECTION("nonexistent sampleset") {
        s = svc->discover_alleles("bogus", range(0, 0, 1000000), als);
        REQUIRE(s == StatusCode::NOT_FOUND);
    }

    SECTION("trio1") {
        s = svc->discover_alleles("discover_alleles_trio1", range(0, 0, 1000000), als);
        REQUIRE(s.ok());

        REQUIRE(als.size() == 7);
        auto p = als.find(allele(range(0, 1000, 1001), "G"));
        REQUIRE(p != als.end());
        REQUIRE(p->second.second == 4);
        p = als.find(allele(range(0, 1010, 1012), "CC"));
        REQUIRE(p != als.end());
        REQUIRE(p->second.second == 3);
    }

    SECTION("trio1 partial") {
        s = svc->discover_alleles("discover_alleles_trio1", range(0, 1009, 1011), als);
        REQUIRE(s.ok());

        REQUIRE(als.size() == 2);
        REQUIRE(als.find(allele(range(0, 1010, 1012), "CC"))->second.second == 3);
    }

    SECTION("2 trios") {
        s = svc->discover_alleles("<ALL>", range(0, 0, 1000000), als);
        REQUIRE(s.ok());

        REQUIRE(als.size() == 10);
        REQUIRE(is_sorted(als.begin(), als.end()));

        REQUIRE(als.find(allele(range(0, 1000, 1001), "A"))->second.second == 6);
        REQUIRE(als.find(allele(range(0, 1000, 1001), "G"))->second.second == 6);
        REQUIRE(als.find(allele(range(0, 1001, 1002), "A"))->second.second == 6);
        REQUIRE(als.find(allele(range(0, 1001, 1002), "C"))->second.second == 2);
        REQUIRE(als.find(allele(range(0, 1001, 1002), "G"))->second.second == 2);
        REQUIRE(als.find(allele(range(0, 1001, 1002), "T"))->second.second == 2);
        REQUIRE(als.find(allele(range(0, 1010, 1012), "AG"))->second.second == 3);
        REQUIRE(als.find(allele(range(0, 1010, 1012), "CC"))->second.second == 3);
        REQUIRE(als.find(allele(range(0, 1010, 1013), "AGA"))->second.second == 2);
        REQUIRE(als.find(allele(range(0, 1010, 1013), "CCC"))->second.second == 4);
    }

    SECTION("spanning allele") {
        s = svc->discover_alleles("<ALL>", range(1, 1010, 1012), als);
        REQUIRE(s.ok());

        REQUIRE(als.size() == 4);
        REQUIRE(is_sorted(als.begin(), als.end()));
     
        REQUIRE(als.find(allele(range(1, 1001, 1016), "AAAAAAAAAAAAAAA"))->second.second == 3);
        // tentative: by default alleles that extend outside of the range of interest get unified to no-call
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

    SECTION("unify_alleles placeholder") {
        s = svc->discover_alleles("<ALL>", range(0, 0, 1000000), als);
        REQUIRE(s.ok());

        vector<unified_site> sites;
        s = unify_alleles(als, sites);
        REQUIRE(s.ok());

        vector<pair<string,size_t> > contigs;
        REQUIRE(data->contigs(contigs).ok());

        REQUIRE(sites.size() > 0);

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
    }
}

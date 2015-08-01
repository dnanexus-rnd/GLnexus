#include <iostream>
#include <algorithm>
#include <vcf.h>
#include "service.h"
#include "alleles.h"
#include "genotyper.h"
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

    Status dataset_bcf(const string& dataset, const range& pos, vector<shared_ptr<bcf1_t> >& records) const override {
        auto p = datasets_.find(dataset);
        if (p == datasets_.end()) {
            return Status::NotFound("unknown data set", dataset);
        }
        //hdr = p->second.header;
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

        REQUIRE(als.size() == 4);
     
        REQUIRE(als.find(allele(range(1, 1001, 1016), "AAAAAAAAAAAAAAA"))->second.observation_count == 3);
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
}

TEST_CASE("unified_sites placeholder") {
    unique_ptr<VCFData> data;
    Status s = VCFData::Open({"discover_alleles_trio1", "discover_alleles_trio2"}, data);
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(data.get(), svc);
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

        REQUIRE(sites.size() == 5);

        REQUIRE(sites[0].pos == range(0,1000,1001));
        REQUIRE(sites[0].alleles.size() == 2);
        REQUIRE(sites[0].alleles[0] == "A");
        REQUIRE(sites[0].alleles[1] == "G");
        REQUIRE(sites[0].unification[make_pair(1000,string("A"))] == 0);
        REQUIRE(sites[0].unification[make_pair(1000,string("G"))] == 1);
        REQUIRE(sites[0].observation_count.size() == 2);
        REQUIRE(sites[0].observation_count[0] == 2);
        REQUIRE(sites[0].observation_count[1] == 4);

        REQUIRE(sites[1].pos == range(0,1001,1002));
        REQUIRE(sites[1].alleles.size() == 3);
        REQUIRE(sites[1].alleles[0] == "C");
        REQUIRE(sites[1].alleles[1] == "G");
        REQUIRE(sites[1].alleles[2] == "T");
        REQUIRE(sites[1].unification[make_pair(1001,string("C"))] == 0);
        REQUIRE(sites[1].unification[make_pair(1001,string("G"))] == 1);
        REQUIRE(sites[1].unification[make_pair(1001,string("T"))] == 2);
        REQUIRE(sites[1].observation_count.size() == 3);
        REQUIRE(sites[1].observation_count[0] == 2);
        REQUIRE(sites[1].observation_count[1] == 2);
        REQUIRE(sites[1].observation_count[2] == 2);

        REQUIRE(sites[2].pos == range(0,1010,1012));
        REQUIRE(sites[2].alleles.size() == 2);
        REQUIRE(sites[2].alleles[0] == "CC");
        REQUIRE(sites[2].alleles[1] == "AG");
        REQUIRE(sites[2].unification[make_pair(1010,string("CC"))] == 0);
        REQUIRE(sites[2].unification[make_pair(1010,string("AG"))] == 1);
        REQUIRE(sites[2].observation_count.size() == 2);
        REQUIRE(sites[2].observation_count[0] == 3);
        REQUIRE(sites[2].observation_count[1] == 3);

        REQUIRE(sites[3].pos == range(0,1100,1101));
        REQUIRE(sites[3].alleles.size() == 2);
        REQUIRE(sites[3].alleles[0] == "C");
        REQUIRE(sites[3].alleles[1] == "A");
        REQUIRE(sites[3].unification[make_pair(1100,string("C"))] == 0);
        REQUIRE(sites[3].unification[make_pair(1100,string("A"))] == 1);
        REQUIRE(sites[3].observation_count.size() == 2);
        REQUIRE(sites[3].observation_count[0] == 3);
        REQUIRE(sites[3].observation_count[1] == 3);

        REQUIRE(sites[4].pos == range(0,1102,1103));
        REQUIRE(sites[4].alleles.size() == 2);
        REQUIRE(sites[4].alleles[0] == "C");
        REQUIRE(sites[4].alleles[1] == "G");
        REQUIRE(sites[4].unification[make_pair(1102,string("C"))] == 0);
        REQUIRE(sites[4].unification[make_pair(1102,string("G"))] == 1);
        REQUIRE(sites[4].observation_count.size() == 2);
        REQUIRE(sites[4].observation_count[0] == 3);
        REQUIRE(sites[4].observation_count[1] == 3);
    }

    SECTION("2 trios") {
        s = svc->discover_alleles("<ALL>", range(0, 0, 1000000), als);
        REQUIRE(s.ok());

        vector<unified_site> sites;
        s = unified_sites(als, sites);
        REQUIRE(s.ok());

        vector<pair<string,size_t> > contigs;
        REQUIRE(data->contigs(contigs).ok());

        REQUIRE(sites.size() == 5);

        REQUIRE(sites[0].pos == range(0,1000,1001));
        REQUIRE(sites[0].alleles.size() == 2);
        REQUIRE(sites[0].alleles[0] == "A");
        REQUIRE(sites[0].alleles[1] == "G");
        REQUIRE(sites[0].unification[make_pair(1000,string("A"))] == 0);
        REQUIRE(sites[0].unification[make_pair(1000,string("G"))] == 1);
        REQUIRE(sites[0].observation_count.size() == 2);
        REQUIRE(sites[0].observation_count[0] == 6);
        REQUIRE(sites[0].observation_count[1] == 6);

        REQUIRE(sites[1].pos == range(0,1001,1002));
        REQUIRE(sites[1].alleles.size() == 4);
        REQUIRE(sites[1].alleles[0] == "C");
        REQUIRE(sites[1].alleles[1] == "A");
        REQUIRE(sites[1].alleles[2] == "G");
        REQUIRE(sites[1].alleles[3] == "T");
        REQUIRE(sites[1].unification[make_pair(1001,string("C"))] == 0);
        REQUIRE(sites[1].unification[make_pair(1001,string("A"))] == 1);
        REQUIRE(sites[1].unification[make_pair(1001,string("G"))] == 2);
        REQUIRE(sites[1].unification[make_pair(1001,string("T"))] == 3);
        REQUIRE(sites[1].observation_count.size() == 4);
        REQUIRE(sites[1].observation_count[0] == 2);
        REQUIRE(sites[1].observation_count[1] == 6);
        REQUIRE(sites[1].observation_count[2] == 2);
        REQUIRE(sites[1].observation_count[3] == 2);

        REQUIRE(sites[2].pos == range(0,1010,1012));
        REQUIRE(sites[2].alleles.size() == 2);
        REQUIRE(sites[2].alleles[0] == "CC");
        REQUIRE(sites[2].alleles[1] == "AG");
        REQUIRE(sites[2].unification[make_pair(1010,string("CC"))] == 0);
        REQUIRE(sites[2].unification[make_pair(1010,string("AG"))] == 1);
        REQUIRE(sites[2].observation_count.size() == 2);
        REQUIRE(sites[2].observation_count[0] == 3);
        REQUIRE(sites[2].observation_count[1] == 3);

        REQUIRE(sites[3].pos == range(0,1100,1101));
        REQUIRE(sites[3].alleles.size() == 2);
        REQUIRE(sites[3].alleles[0] == "C");
        REQUIRE(sites[3].alleles[1] == "A");
        REQUIRE(sites[3].unification[make_pair(1100,string("C"))] == 0);
        REQUIRE(sites[3].unification[make_pair(1100,string("A"))] == 1);
        REQUIRE(sites[3].observation_count.size() == 2);
        REQUIRE(sites[3].observation_count[0] == 3);
        REQUIRE(sites[3].observation_count[1] == 3);

        REQUIRE(sites[4].pos == range(0,1102,1103));
        REQUIRE(sites[4].alleles.size() == 2);
        REQUIRE(sites[4].alleles[0] == "C");
        REQUIRE(sites[4].alleles[1] == "G");
        REQUIRE(sites[4].unification[make_pair(1102,string("C"))] == 0);
        REQUIRE(sites[4].unification[make_pair(1102,string("G"))] == 1);
        REQUIRE(sites[4].observation_count.size() == 2);
        REQUIRE(sites[4].observation_count[0] == 3);
        REQUIRE(sites[4].observation_count[1] == 3);

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

TEST_CASE("genotyper placeholder") {
    unique_ptr<VCFData> data;
    Status s = VCFData::Open({"discover_alleles_trio1", "discover_alleles_trio2"}, data);
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(data.get(), svc);
    REQUIRE(s.ok());

    discovered_alleles als;

    const string tfn("/tmp/GLnexus_unit_tests.bcf");

    SECTION("discover_alleles_trio1") {
        s = svc->discover_alleles("discover_alleles_trio1", range(0, 0, 1000000), als);
        REQUIRE(s.ok());

        vector<unified_site> sites;
        s = unified_sites(als, sites);
        REQUIRE(s.ok());

        s = svc->genotype_sites(string("discover_alleles_trio1"), sites, tfn);
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
        REQUIRE(bcf_gt_allele(gt[0] == 1));
        REQUIRE(bcf_gt_allele(gt[1] == 1));
        REQUIRE(bcf_gt_allele(gt[2] == 0));
        REQUIRE(bcf_gt_allele(gt[3] == 1));
        REQUIRE(bcf_gt_allele(gt[4] == 0));
        REQUIRE(bcf_gt_allele(gt[5] == 1));

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 3);
        REQUIRE(string(record->d.allele[0]) == "C");
        REQUIRE(string(record->d.allele[1]) == "G");
        REQUIRE(string(record->d.allele[2]) == "T");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 6);
        REQUIRE(bcf_gt_allele(gt[0] == 2));
        REQUIRE(bcf_gt_allele(gt[1] == 2));
        REQUIRE(bcf_gt_allele(gt[2] == 0));
        REQUIRE(bcf_gt_allele(gt[3] == 0));
        REQUIRE(bcf_gt_allele(gt[4] == 1));
        REQUIRE(bcf_gt_allele(gt[5] == 1));

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "CC");
        REQUIRE(string(record->d.allele[1]) == "AG");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 6);
        REQUIRE(bcf_gt_allele(gt[0] == 0));
        REQUIRE(bcf_gt_allele(gt[1] == 1));
        REQUIRE(bcf_gt_allele(gt[2] == 0));
        REQUIRE(bcf_gt_allele(gt[3] == 1));
        REQUIRE(bcf_gt_allele(gt[4] == 0));
        REQUIRE(bcf_gt_allele(gt[5] == 1));

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "C");
        REQUIRE(string(record->d.allele[1]) == "A");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 6);
        REQUIRE(bcf_gt_allele(gt[0] == 0));
        REQUIRE(bcf_gt_allele(gt[1] == 1));
        REQUIRE(bcf_gt_allele(gt[2] == 0));
        REQUIRE(bcf_gt_allele(gt[3] == 1));
        REQUIRE(bcf_gt_allele(gt[4] == 0));
        REQUIRE(bcf_gt_allele(gt[5] == 1));

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "C");
        REQUIRE(string(record->d.allele[1]) == "G");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 6);
        REQUIRE(bcf_gt_allele(gt[0] == 0));
        REQUIRE(bcf_gt_allele(gt[1] == 1));
        REQUIRE(bcf_gt_allele(gt[2] == 0));
        REQUIRE(bcf_gt_allele(gt[3] == 1));
        REQUIRE(bcf_gt_allele(gt[4] == 0));
        REQUIRE(bcf_gt_allele(gt[5] == 1));

        free(gt);
    }

    SECTION("2 trios") {
        s = svc->discover_alleles("<ALL>", range(0, 0, 1000000), als);
        REQUIRE(s.ok());

        vector<unified_site> sites;
        s = unified_sites(als, sites);
        REQUIRE(s.ok());

        s = svc->genotype_sites(string("<ALL>"), sites, tfn);
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

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "A");
        REQUIRE(string(record->d.allele[1]) == "G");

        int *gt = nullptr, gtsz = 0;
        int nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 12);
        REQUIRE(bcf_gt_allele(gt[0] == 1));
        REQUIRE(bcf_gt_allele(gt[1] == 1));
        REQUIRE(bcf_gt_allele(gt[2] == 0));
        REQUIRE(bcf_gt_allele(gt[3] == 1));
        REQUIRE(bcf_gt_allele(gt[4] == 0));
        REQUIRE(bcf_gt_allele(gt[5] == 1));
        REQUIRE(bcf_gt_allele(gt[6] == 0));
        REQUIRE(bcf_gt_allele(gt[7] == 0));
        REQUIRE(bcf_gt_allele(gt[8] == 0));
        REQUIRE(bcf_gt_allele(gt[9] == 1));
        REQUIRE(bcf_gt_allele(gt[10] == 0));
        REQUIRE(bcf_gt_allele(gt[11] == 1));

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 4);
        REQUIRE(string(record->d.allele[0]) == "C");
        REQUIRE(string(record->d.allele[1]) == "A");
        REQUIRE(string(record->d.allele[2]) == "G");
        REQUIRE(string(record->d.allele[3]) == "T");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 12);
        REQUIRE(bcf_gt_allele(gt[0] == 3));
        REQUIRE(bcf_gt_allele(gt[1] == 3));
        REQUIRE(bcf_gt_allele(gt[2] == 0));
        REQUIRE(bcf_gt_allele(gt[3] == 0));
        REQUIRE(bcf_gt_allele(gt[4] == 2));
        REQUIRE(bcf_gt_allele(gt[5] == 2));
        REQUIRE(bcf_gt_allele(gt[6] == 1));
        REQUIRE(bcf_gt_allele(gt[7] == 1));
        REQUIRE(bcf_gt_allele(gt[8] == 1));
        REQUIRE(bcf_gt_allele(gt[9] == 1));
        REQUIRE(bcf_gt_allele(gt[10] == 1));
        REQUIRE(bcf_gt_allele(gt[11] == 1));

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "CC");
        REQUIRE(string(record->d.allele[1]) == "AG");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 12);
        REQUIRE(bcf_gt_allele(gt[0] == 0));
        REQUIRE(bcf_gt_allele(gt[1] == 1));
        REQUIRE(bcf_gt_allele(gt[2] == 0));
        REQUIRE(bcf_gt_allele(gt[3] == 1));
        REQUIRE(bcf_gt_allele(gt[4] == 0));
        REQUIRE(bcf_gt_allele(gt[5] == 1));
        REQUIRE(bcf_gt_is_missing(gt[6]));
        REQUIRE(bcf_gt_is_missing(gt[7]));
        REQUIRE(bcf_gt_is_missing(gt[8]));
        REQUIRE(bcf_gt_is_missing(gt[9]));
        REQUIRE(bcf_gt_is_missing(gt[10]));
        REQUIRE(bcf_gt_is_missing(gt[11]));

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "C");
        REQUIRE(string(record->d.allele[1]) == "A");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 12);
        REQUIRE(bcf_gt_allele(gt[0] == 0));
        REQUIRE(bcf_gt_allele(gt[1] == 1));
        REQUIRE(bcf_gt_allele(gt[2] == 0));
        REQUIRE(bcf_gt_allele(gt[3] == 1));
        REQUIRE(bcf_gt_allele(gt[4] == 0));
        REQUIRE(bcf_gt_allele(gt[5] == 1));
        // TODO: ought to be able to make reference half-calls in trio2
        REQUIRE(bcf_gt_is_missing(gt[6]));
        REQUIRE(bcf_gt_is_missing(gt[7]));
        REQUIRE(bcf_gt_is_missing(gt[8]));
        REQUIRE(bcf_gt_is_missing(gt[9]));
        REQUIRE(bcf_gt_is_missing(gt[10]));
        REQUIRE(bcf_gt_is_missing(gt[11]));

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "C");
        REQUIRE(string(record->d.allele[1]) == "G");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 12);
        REQUIRE(bcf_gt_allele(gt[0] == 0));
        REQUIRE(bcf_gt_allele(gt[1] == 1));
        REQUIRE(bcf_gt_allele(gt[2] == 0));
        REQUIRE(bcf_gt_allele(gt[3] == 1));
        REQUIRE(bcf_gt_allele(gt[4] == 0));
        REQUIRE(bcf_gt_allele(gt[5] == 1));
        // TODO: ought to be able to make reference half-calls in trio2
        REQUIRE(bcf_gt_is_missing(gt[6]));
        REQUIRE(bcf_gt_is_missing(gt[7]));
        REQUIRE(bcf_gt_is_missing(gt[8]));
        REQUIRE(bcf_gt_is_missing(gt[9]));
        REQUIRE(bcf_gt_is_missing(gt[10]));
        REQUIRE(bcf_gt_is_missing(gt[11]));

        free(gt);
    }

    SECTION("trio2 with all alleles") {
        s = svc->discover_alleles("<ALL>", range(0, 0, 1000000), als);
        REQUIRE(s.ok());

        vector<unified_site> sites;
        s = unified_sites(als, sites);
        REQUIRE(s.ok());

        s = svc->genotype_sites(string("discover_alleles_trio2"), sites, tfn);
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
        REQUIRE(bcf_gt_allele(gt[0] == 0));
        REQUIRE(bcf_gt_allele(gt[1] == 0));
        REQUIRE(bcf_gt_allele(gt[2] == 0));
        REQUIRE(bcf_gt_allele(gt[3] == 1));
        REQUIRE(bcf_gt_allele(gt[4] == 0));
        REQUIRE(bcf_gt_allele(gt[5] == 1));

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 4);
        REQUIRE(string(record->d.allele[0]) == "C");
        REQUIRE(string(record->d.allele[1]) == "A");
        REQUIRE(string(record->d.allele[2]) == "G");
        REQUIRE(string(record->d.allele[3]) == "T");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 6);
        REQUIRE(bcf_gt_allele(gt[0] == 1));
        REQUIRE(bcf_gt_allele(gt[1] == 1));
        REQUIRE(bcf_gt_allele(gt[2] == 1));
        REQUIRE(bcf_gt_allele(gt[3] == 1));
        REQUIRE(bcf_gt_allele(gt[4] == 1));
        REQUIRE(bcf_gt_allele(gt[5] == 1));

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "CC");
        REQUIRE(string(record->d.allele[1]) == "AG");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 6);
        REQUIRE(bcf_gt_is_missing(gt[0]));
        REQUIRE(bcf_gt_is_missing(gt[1]));
        REQUIRE(bcf_gt_is_missing(gt[2]));
        REQUIRE(bcf_gt_is_missing(gt[3]));
        REQUIRE(bcf_gt_is_missing(gt[4]));
        REQUIRE(bcf_gt_is_missing(gt[5]));

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "C");
        REQUIRE(string(record->d.allele[1]) == "A");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 6);
        REQUIRE(bcf_gt_is_missing(gt[0]));
        REQUIRE(bcf_gt_is_missing(gt[1]));
        REQUIRE(bcf_gt_is_missing(gt[2]));
        REQUIRE(bcf_gt_is_missing(gt[3]));
        REQUIRE(bcf_gt_is_missing(gt[4]));
        REQUIRE(bcf_gt_is_missing(gt[5]));

        REQUIRE(bcf_read(vcf.get(), hdr.get(), record.get()) == 0);
        REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);

        REQUIRE(record->n_allele == 2);
        REQUIRE(string(record->d.allele[0]) == "C");
        REQUIRE(string(record->d.allele[1]) == "G");

        nGT = bcf_get_genotypes(hdr.get(), record.get(), &gt, &gtsz);
        REQUIRE(nGT == 6);
        REQUIRE(bcf_gt_is_missing(gt[0]));
        REQUIRE(bcf_gt_is_missing(gt[1]));
        REQUIRE(bcf_gt_is_missing(gt[2]));
        REQUIRE(bcf_gt_is_missing(gt[3]));
        REQUIRE(bcf_gt_is_missing(gt[4]));
        REQUIRE(bcf_gt_is_missing(gt[5]));

        free(gt);
    }
}

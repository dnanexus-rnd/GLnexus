#include <iostream>
#include <algorithm>
#include <vcf.h>
#include "service.h"
#include "unifier.h"
#include "genotyper.h"
#include "utils.cc"
#include "catch.hpp"
using namespace std;
using namespace GLnexus;

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

    Status dataset_header(const string& dataset, shared_ptr<const bcf_hdr_t>* hdr) override {
        if (i_++ % fail_every_ == 0) {
            failed_once_ = true;
            return Status::IOError("SIM");
        }
        return inner_.dataset_header(dataset, hdr);
    }

    Status dataset_range(const string& dataset, const bcf_hdr_t *hdr, const range& pos,
                         bcf_predicate predicate, vector<shared_ptr<bcf1_t>>* records) override {
        if (i_++ % fail_every_ == 0) {
            failed_once_ = true;
            return Status::IOError("SIM");
        }
        return inner_.dataset_range(dataset, hdr, pos, predicate, records);
    }

    bool failed_once() { return failed_once_; }
};

TEST_CASE("service::discover_alleles") {
    unique_ptr<VCFData> data;
    Status s = VCFData::Open({"discover_alleles_trio1.vcf", "discover_alleles_trio2.vcf"}, data);
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(service_config(), *data, *data, svc);
    REQUIRE(s.ok());

    discovered_alleles als;
    unsigned N;

    SECTION("multiple ranges") {
        vector<range> ranges;
        ranges.push_back(range(0, 1000, 1001));
        ranges.push_back(range(0, 1001, 1002));
        ranges.push_back(range(0, 1010, 1013));
        ranges.push_back(range(1, 1000, 1001));
        ranges.push_back(range(1, 1010, 1012));
        ranges.push_back(range(1, 2000, 2100));
        ranges.push_back(range(2, 1001, 1002));
        vector<discovered_alleles> mals;
        s = svc->discover_alleles("<ALL>", ranges, N, mals);
        REQUIRE(s.ok());
        REQUIRE(mals.size() == ranges.size());
        REQUIRE(N == 6);

        REQUIRE(mals[0].size() == 2);
        REQUIRE(mals[0].find(allele(range(0, 1000, 1001), "A"))->second.zGQ.copy_number() == 6);
        REQUIRE(mals[0].find(allele(range(0, 1000, 1001), "G"))->second.zGQ.copy_number() == 6);
        REQUIRE(mals[0].find(allele(range(0, 1000, 1001), "A"))->second.in_target == range(0, 1000, 1001));
        REQUIRE(mals[0].find(allele(range(0, 1000, 1001), "G"))->second.in_target == range(0, 1000, 1001));

        REQUIRE(mals[1].size() == 4);
        REQUIRE(mals[1].find(allele(range(0, 1001, 1002), "A"))->second.zGQ.copy_number() == 6);
        REQUIRE(mals[1].find(allele(range(0, 1001, 1002), "C"))->second.zGQ.copy_number() == 2);
        REQUIRE(mals[1].find(allele(range(0, 1001, 1002), "G"))->second.zGQ.copy_number() == 2);
        REQUIRE(mals[1].find(allele(range(0, 1001, 1002), "T"))->second.zGQ.copy_number() == 2);

        REQUIRE(mals[2].size() == 4);
        REQUIRE(mals[2].find(allele(range(0, 1010, 1012), "AG"))->second.zGQ.copy_number() == 3);
        REQUIRE(mals[2].find(allele(range(0, 1010, 1012), "CC"))->second.zGQ.copy_number() == 3);
        REQUIRE(mals[2].find(allele(range(0, 1010, 1013), "AGA"))->second.zGQ.copy_number() == 2);
        REQUIRE(mals[2].find(allele(range(0, 1010, 1013), "CCC"))->second.zGQ.copy_number() == 4);

        REQUIRE(mals[3].size() == 2);
        REQUIRE(mals[3].find(allele(range(1, 1000, 1001), "A"))->second.zGQ.copy_number() == 2);
        REQUIRE(mals[3].find(allele(range(1, 1000, 1001), "AA"))->second.zGQ.copy_number() == 4);

        REQUIRE(mals[4].size() == 4);
        REQUIRE(mals[4].find(allele(range(1, 1010, 1012), "AG"))->second.zGQ.copy_number() == 3);
        REQUIRE(mals[4].find(allele(range(1, 1010, 1012), "CC"))->second.zGQ.copy_number() == 3);
        REQUIRE(mals[4].find(allele(range(1, 1010, 1012), "AG"))->second.in_target == range(1, 1010, 1012));
        REQUIRE(mals[4].find(allele(range(1, 1010, 1012), "CC"))->second.in_target == range(1, 1010, 1012));
        REQUIRE(mals[2].find(allele(range(0, 1010, 1013), "AGA"))->second.zGQ.copy_number() == 2);
        REQUIRE(mals[2].find(allele(range(0, 1010, 1013), "CCC"))->second.zGQ.copy_number() == 4);

        REQUIRE(mals[5].empty());

        REQUIRE(mals[6].empty());
    }

    SECTION("allele overlaps two ranges") {
        vector<range> ranges;
        ranges.push_back(range(1, 1000, 1002));
        ranges.push_back(range(1, 1010, 1014));
        vector<discovered_alleles> mals;
        s = svc->discover_alleles("<ALL>", ranges, N, mals);
        REQUIRE(s.ok());
        REQUIRE(mals.size() == ranges.size());
        REQUIRE(N == 6);

        // a long allele is found in both ranges
        REQUIRE(mals[0].size() == 6);
        REQUIRE(mals[0].find(allele(range(1, 1001, 1016), "AAAAAAAAAAAAAAA"))->second.zGQ.copy_number() == 3);
        REQUIRE(mals[0].find(allele(range(1, 1001, 1016), "AAAAAAAAAAAAAAA"))->second.in_target == range(1, 1000, 1002));
        REQUIRE(mals[1].size() == 4);
        REQUIRE(mals[1].find(allele(range(1, 1001, 1016), "AAAAAAAAAAAAAAA"))->second.zGQ.copy_number() == 3);
        REQUIRE(mals[1].find(allele(range(1, 1001, 1016), "AAAAAAAAAAAAAAA"))->second.in_target == range(1, 1010, 1014));

        // merge_discovered_alleles should collapse them without double-counting
        discovered_alleles merged;
        REQUIRE(merge_discovered_alleles(mals[0], merged).ok());
        REQUIRE(merge_discovered_alleles(mals[1], merged).ok());
        REQUIRE(merged.find(allele(range(1, 1001, 1016), "AAAAAAAAAAAAAAA"))->second.zGQ.copy_number() == 3);
        REQUIRE(merged.find(allele(range(1, 1001, 1016), "AAAAAAAAAAAAAAA"))->second.in_target == range(1, 1010, 1014));
    }

    SECTION("simulate I/O errors - single") {
        unique_ptr<SimFailBCFData> faildata;
        bool worked = false;

        for (size_t fail_every = 1; fail_every < 100; fail_every++) {
            s = SimFailBCFData::Open(*data, fail_every, faildata);
            REQUIRE(s.ok());

            s = Service::Start(service_config(), *data, *faildata, svc);
            REQUIRE(s.ok());

            s = svc->discover_alleles("<ALL>", range(0, 0, 1099), N, als);
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

            s = Service::Start(service_config(), *data, *faildata, svc);
            REQUIRE(s.ok());

            vector<range> ranges;
            ranges.push_back(range(0, 1000, 1001));
            ranges.push_back(range(0, 1001, 1002));
            ranges.push_back(range(0, 1010, 1013));
            ranges.push_back(range(1, 1000, 1001));
            ranges.push_back(range(1, 1010, 1012));
            vector<discovered_alleles> mals;
            s = svc->discover_alleles("<ALL>", ranges, N, mals);
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
    Status s = VCFData::Open({"NA12878D_HiSeqX.21.10009462-10009469.gvcf", "NA12878D_HiSeqX.21.10009462-10009469.bogus.gvcf"}, data, false);
    cout << s.str() << endl;
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(service_config(), *data, *data, svc);
    REQUIRE(s.ok());

    discovered_alleles als;
    unsigned N;

    SECTION("exclude symbolic alleles") {
        s = svc->discover_alleles("<ALL>", range(0, 10000000, 10010000), N, als);

        REQUIRE(N == 2);
        REQUIRE(als.size() == 2);
        auto p = als.find(allele(range(0, 10009463, 10009465), "TA"));
        REQUIRE(p != als.end());
        REQUIRE(p->second.zGQ.copy_number() == 1);
        p = als.find(allele(range(0, 10009463, 10009465), "T"));
        REQUIRE(p != als.end());
        REQUIRE(p->second.zGQ.copy_number() == 1);
        REQUIRE(p->second.in_target == range(0, 10000000, 10010000));
    }

    SECTION("exclusion/detection of bogus alleles") {
        s = svc->discover_alleles("<ALL>", range(1, 10009463, 10009465), N, als);

        REQUIRE(N == 2);
        REQUIRE(als.size() == 2);
        auto p = als.find(allele(range(1, 10009463, 10009465), "TA"));
        REQUIRE(p != als.end());
        REQUIRE(p->second.zGQ.copy_number() == 1);
        p = als.find(allele(range(1, 10009463, 10009465), "T"));
        REQUIRE(p != als.end());
        REQUIRE(p->second.zGQ.copy_number() == 1);

        s = svc->discover_alleles("<ALL>", range(1, 10009465, 10009466), N, als);
        REQUIRE(s == StatusCode::INVALID);
    }
}

TEST_CASE("unified_sites") {
    unique_ptr<VCFData> data;
    Status s = VCFData::Open({"discover_alleles_trio1.vcf", "discover_alleles_trio2.vcf"}, data);
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(service_config(), *data, *data, svc);
    REQUIRE(s.ok());

    discovered_alleles als;
    unsigned N;

    SECTION("trio1") {
        s = svc->discover_alleles("discover_alleles_trio1", range(0, 0, 1000000), N, als);
        REQUIRE(s.ok());
        REQUIRE(N == 3);

        vector<unified_site> sites;
        unifier_stats stats;
        s = unified_sites(unifier_config(), N, als, sites, stats);
        REQUIRE(s.ok());

        vector<pair<string,size_t> > contigs;
        REQUIRE(data->contigs(contigs).ok());

        REQUIRE(sites.size() == 6);

        REQUIRE(sites[0].pos == range(0,1000,1001));
        REQUIRE(sites[0].alleles.size() == 2);
        REQUIRE(sites[0].alleles[0].dna == "A");
        REQUIRE(sites[0].alleles[1].dna == "G");
        REQUIRE(sites[0].unification[allele(range(0,1000,1001),string("A"))] == 0);
        REQUIRE(sites[0].unification[allele(range(0,1000,1001),string("G"))] == 1);
        REQUIRE(sites[0].alleles[0].frequency != sites[0].alleles[0].frequency);
        REQUIRE(sites[0].alleles[1].frequency == 0.666667f);

        REQUIRE(sites[1].pos == range(0,1001,1002));
        REQUIRE(sites[1].alleles.size() == 3);
        REQUIRE(sites[1].alleles[0].dna == "C");
        REQUIRE(sites[1].alleles[1].dna == "G");
        REQUIRE(sites[1].alleles[2].dna == "T");
        REQUIRE(sites[1].unification[allele(range(0,1001,1002),string("C"))] == 0);
        REQUIRE(sites[1].unification[allele(range(0,1001,1002),string("G"))] == 1);
        REQUIRE(sites[1].unification[allele(range(0,1001,1002),string("T"))] == 2);
        REQUIRE(sites[1].alleles[0].frequency != sites[1].alleles[0].frequency);
        REQUIRE(sites[1].alleles[1].frequency == 0.333333f);
        REQUIRE(sites[1].alleles[2].frequency == 0.333333f);

        REQUIRE(sites[2].pos == range(0,1010,1012));
        REQUIRE(sites[2].alleles.size() == 2);
        REQUIRE(sites[2].alleles[0].dna == "CC");
        REQUIRE(sites[2].alleles[1].dna == "AG");
        REQUIRE(sites[2].unification[allele(range(0,1010,1012),string("CC"))] == 0);
        REQUIRE(sites[2].unification[allele(range(0,1010,1012),string("AG"))] == 1);
        REQUIRE(sites[2].alleles[0].frequency != sites[2].alleles[0].frequency);
        REQUIRE(sites[2].alleles[1].frequency == 0.5f);

        REQUIRE(sites[3].pos == range(0,1100,1101));
        REQUIRE(sites[3].alleles.size() == 2);
        REQUIRE(sites[3].alleles[0].dna == "C");
        REQUIRE(sites[3].alleles[1].dna == "A");
        REQUIRE(sites[3].unification[allele(range(0,1100,1101),string("C"))] == 0);
        REQUIRE(sites[3].unification[allele(range(0,1100,1101),string("A"))] == 1);
        REQUIRE(sites[3].alleles[0].frequency != sites[3].alleles[0].frequency);
        REQUIRE(sites[3].alleles[1].frequency == 0.5f);

        REQUIRE(sites[4].pos == range(0,1102,1103));
        REQUIRE(sites[4].alleles.size() == 2);
        REQUIRE(sites[4].alleles[0].dna == "C");
        REQUIRE(sites[4].alleles[1].dna == "G");
        REQUIRE(sites[4].unification[allele(range(0,1102,1103),string("C"))] == 0);
        REQUIRE(sites[4].unification[allele(range(0,1102,1103),string("G"))] == 1);
        REQUIRE(sites[4].alleles[0].frequency != sites[4].alleles[0].frequency);
        REQUIRE(sites[4].alleles[1].frequency == 0.5f);

        REQUIRE(sites[5].pos == range(0,1200,1201));
        REQUIRE(sites[5].alleles.size() == 2);
        REQUIRE(sites[5].alleles[0].dna == "C");
        REQUIRE(sites[5].alleles[1].dna == "A");
        REQUIRE(sites[5].unification[allele(range(0,1200,1201),string("C"))] == 0);
        REQUIRE(sites[5].unification[allele(range(0,1200,1201),string("A"))] == 1);
        REQUIRE(sites[5].alleles[0].frequency != sites[5].alleles[0].frequency);
        REQUIRE(sites[5].alleles[1].frequency == 0.5f);
    }

    SECTION("2 trios") {
        s = svc->discover_alleles("<ALL>", range(0, 0, 1000000), N, als);
        REQUIRE(s.ok());
        REQUIRE(N == 6);

        vector<unified_site> sites;
        unifier_stats stats;
        s = unified_sites(unifier_config(), N, als, sites, stats);
        REQUIRE(s.ok());

        vector<pair<string,size_t> > contigs;
        REQUIRE(data->contigs(contigs).ok());

        REQUIRE(sites.size() == 6);

        REQUIRE(sites[0].pos == range(0,1000,1001));
        REQUIRE(sites[0].alleles.size() == 2);
        REQUIRE(sites[0].alleles[0].dna == "A");
        REQUIRE(sites[0].alleles[1].dna == "G");
        REQUIRE(sites[0].unification[allele(range(0,1000,1001),string("A"))] == 0);
        REQUIRE(sites[0].unification[allele(range(0,1000,1001),string("G"))] == 1);
        REQUIRE(sites[0].alleles[0].frequency != sites[0].alleles[0].frequency);
        REQUIRE(sites[0].alleles[1].frequency == 0.5f);

        REQUIRE(sites[1].pos == range(0,1001,1002));
        REQUIRE(sites[1].alleles.size() == 4);
        REQUIRE(sites[1].alleles[0].dna == "C");
        REQUIRE(sites[1].alleles[1].dna == "A");
        REQUIRE(sites[1].alleles[2].dna == "G");
        REQUIRE(sites[1].alleles[3].dna == "T");
        REQUIRE(sites[1].unification[allele(range(0,1001,1002),string("C"))] == 0);
        REQUIRE(sites[1].unification[allele(range(0,1001,1002),string("A"))] == 1);
        REQUIRE(sites[1].unification[allele(range(0,1001,1002),string("G"))] == 2);
        REQUIRE(sites[1].unification[allele(range(0,1001,1002),string("T"))] == 3);
        REQUIRE(sites[1].alleles[0].frequency != sites[1].alleles[0].frequency);
        REQUIRE(sites[1].alleles[1].frequency == 0.5f);
        REQUIRE(sites[1].alleles[2].frequency == 0.166667f);
        REQUIRE(sites[1].alleles[3].frequency == 0.166667f);

        REQUIRE(sites[2].pos == range(0,1010,1013));
        REQUIRE(sites[2].alleles.size() == 3);
        REQUIRE(sites[2].alleles[0].dna == "CCC");
        REQUIRE(sites[2].alleles[1].dna == "AGC");
        REQUIRE(sites[2].alleles[2].dna == "AGA");
        REQUIRE(sites[2].unification[allele(range(0,1010,1012),string("CC"))] == 0);
        REQUIRE(sites[2].unification[allele(range(0,1010,1013),string("CCC"))] == 0);
        REQUIRE(sites[2].unification[allele(range(0,1010,1012),string("AG"))] == 1);
        REQUIRE(sites[2].unification[allele(range(0,1010,1013),string("AGA"))] == 2);
        REQUIRE(sites[2].alleles[0].frequency != sites[2].alleles[0].frequency);
        REQUIRE(sites[2].alleles[1].frequency == 0.25f);
        REQUIRE(sites[2].alleles[2].frequency == 0.166667f);

        REQUIRE(sites[3].pos == range(0,1100,1101));
        REQUIRE(sites[3].alleles.size() == 2);
        REQUIRE(sites[3].alleles[0].dna == "C");
        REQUIRE(sites[3].alleles[1].dna == "A");
        REQUIRE(sites[3].unification[allele(range(0,1100,1101),string("C"))] == 0);
        REQUIRE(sites[3].unification[allele(range(0,1100,1101),string("A"))] == 1);
        REQUIRE(sites[3].alleles[0].frequency != sites[3].alleles[0].frequency);
        REQUIRE(sites[3].alleles[1].frequency == 0.25f);

        REQUIRE(sites[4].pos == range(0,1102,1103));
        REQUIRE(sites[4].alleles.size() == 2);
        REQUIRE(sites[4].alleles[0].dna == "C");
        REQUIRE(sites[4].alleles[1].dna == "G");
        REQUIRE(sites[4].unification[allele(range(0,1102,1103),string("C"))] == 0);
        REQUIRE(sites[4].unification[allele(range(0,1102,1103),string("G"))] == 1);
        REQUIRE(sites[4].alleles[0].frequency != sites[4].alleles[0].frequency);
        REQUIRE(sites[4].alleles[1].frequency == 0.25f);

        REQUIRE(sites[5].pos == range(0,1200,1201));
        REQUIRE(sites[5].alleles.size() == 2);
        REQUIRE(sites[5].alleles[0].dna == "C");
        REQUIRE(sites[5].alleles[1].dna == "A");
        REQUIRE(sites[5].unification[allele(range(0,1200,1201),string("C"))] == 0);
        REQUIRE(sites[5].unification[allele(range(0,1200,1201),string("A"))] == 1);
        REQUIRE(sites[5].alleles[0].frequency != sites[5].alleles[0].frequency);
        REQUIRE(sites[5].alleles[1].frequency == 0.25f);

        // An allele from trio2 that would collapse sites[3] and sites[4] was
        // pruned.

        /*
        for (const auto& site : sites) {
            cout << site.pos.str(contigs);
            for (unsigned i = 0; i < site.alleles.size(); i++) {
                cout << ' ' << site.alleles[i] << '[' << site.allele_frequencies[i] << ']';
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
    Status s = VCFData::Open({"discover_alleles_trio1.vcf", "discover_alleles_trio2.vcf"}, data);
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(service_config(), *data, *data, svc);
    REQUIRE(s.ok());

    discovered_alleles als;
    unsigned N;

    const string tfn("/tmp/GLnexus_unit_tests.bcf");

    SECTION("simulate I/O errors") {
        s = Service::Start(service_config(), *data, *data, svc);
        REQUIRE(s.ok());
        s = svc->discover_alleles("<ALL>", range(0, 0, 1000000), N, als);
        REQUIRE(s.ok());
        REQUIRE(N == 6);
        vector<unified_site> sites;
        unifier_stats stats;
        s = unified_sites(unifier_config(), N, als, sites, stats);
        REQUIRE(s.ok());

        unique_ptr<SimFailBCFData> faildata;
        bool worked = false;

        for (size_t fail_every = 1; fail_every < 100; fail_every++) {
            s = SimFailBCFData::Open(*data, fail_every, faildata);
            REQUIRE(s.ok());

            s = Service::Start(service_config(), *data, *faildata, svc);
            REQUIRE(s.ok());

            s = svc->genotype_sites(genotyper_config(), string("<ALL>"), sites, tfn);
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

    SECTION("unification with multiple contigs") {
        discovered_alleles als0, als1;
        unsigned N;

        s = svc->discover_alleles("<ALL>", range(0, 0, 1000000), N, als0);
        REQUIRE(s.ok());

        s = svc->discover_alleles("<ALL>", range(1, 0, 1000000), N, als1);
        REQUIRE(s.ok());

        als.clear();
        REQUIRE(merge_discovered_alleles(als0, als).ok());
        REQUIRE(merge_discovered_alleles(als1, als).ok());

        vector<unified_site> sites;
        unifier_stats stats;
        s = unified_sites(unifier_config(), N, als, sites, stats);
        //cout << s.str() << endl;
        REQUIRE(s.ok());

        REQUIRE(is_sorted(sites.begin(), sites.end()));
        REQUIRE(sites[0].pos.rid == 0);
        REQUIRE(sites[sites.size()-1].pos.rid == 1);
    }
}

TEST_CASE("gVCF genotyper") {
    unique_ptr<VCFData> data;
    Status s = VCFData::Open({"NA12878D_HiSeqX.21.10009462-10009469.gvcf"}, data);
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(service_config(), *data, *data, svc);
    REQUIRE(s.ok());

    discovered_alleles als;

    const string tfn("/tmp/GLnexus_unit_tests.bcf");

    SECTION("no depth requirement") {
        vector<unified_site> sites;

        unified_site us(range(0,10009461,10009462));
        us.unification[allele(range(0,10009461,10009462),"T")] = 0;
        us.unification[allele(range(0,10009461,10009462),"A")] = 1;
        us.alleles.push_back(unified_allele(range(0,10009461,10009462),"T"));
        us.alleles.push_back(unified_allele(range(0,10009461,10009462),"A"));
        us.alleles[1].frequency = 0.1;
        sites.push_back(us);

        us.pos = range(0,10009462,10009463);
        us.alleles.clear();
        us.unification.clear();
        us.unification[allele(range(0,10009462,10009463),"C")] = 0;
        us.unification[allele(range(0,10009462,10009463),"G")] = 1;
        us.alleles.push_back(unified_allele(range(0,10009462,10009463),"C"));
        us.alleles.push_back(unified_allele(range(0,10009462,10009463),"G"));
        us.alleles[1].frequency = 0.1;
        sites.push_back(us);

        // this site spans two gVCF records and will not be called by the current algorithm.
        us.pos = range(0,10009465,10009467);
        us.alleles.clear();
        us.unification.clear();
        us.unification[allele(range(0,10009465,10009467),"AA")] = 0;
        us.unification[allele(range(0,10009465,10009467),"GT")] = 1;
        us.alleles.push_back(unified_allele(range(0,10009465,10009467),"AA"));
        us.alleles.push_back(unified_allele(range(0,10009465,10009467),"GT"));
        us.alleles[1].frequency = 0.1;
        sites.push_back(us);

        s = svc->genotype_sites(genotyper_config(), string("NA12878D_HiSeqX.21.10009462-10009469"), sites, tfn);
        REQUIRE(s.ok());
    }

    SECTION("require depth > 12") {
        vector<unified_site> sites;

        unified_site us(range(0,10009463,10009465));
        us.unification[allele(range(0,10009463,10009465),"TA")] = 0;
        us.unification[allele(range(0,10009463,10009465),"T")] = 1;
        us.alleles.push_back(unified_allele(range(0,10009463,10009465),"TA"));
        us.alleles.push_back(unified_allele(range(0,10009463,10009465),"T"));
        us.alleles[1].frequency = 0.1;
        sites.push_back(us);

        us.pos = range(0,10009465,10009466);
        us.alleles.clear();
        us.unification.clear();
        us.unification[allele(range(0,10009465,10009466),"A")] = 0;
        us.unification[allele(range(0,10009465,10009466),"G")] = 1;
        us.alleles.push_back(unified_allele(range(0,10009465,10009466),"A"));
        us.alleles.push_back(unified_allele(range(0,10009465,10009466),"G"));
        us.alleles[1].frequency = 0.1;
        sites.push_back(us);

        us.pos = range(0,10009466,10009467);
        us.alleles.clear();
        us.unification.clear();
        us.unification[allele(range(0,10009466,10009467),"A")] = 0;
        us.unification[allele(range(0,10009466,10009467),"C")] = 1;
        us.alleles.push_back(unified_allele(range(0,10009466,10009467),"A"));
        us.alleles.push_back(unified_allele(range(0,10009466,10009467),"C"));
        us.alleles[1].frequency = 0.1;
        sites.push_back(us);

        genotyper_config cfg;
        cfg.required_dp = 13;

        s = svc->genotype_sites(cfg, string("NA12878D_HiSeqX.21.10009462-10009469"), sites, tfn);
        REQUIRE(s.ok());
    }


    SECTION("require depth > 9") {
        vector<unified_site> sites;

        unified_site us(range(0,10009463,10009465));
        us.unification[allele(range(0,10009463,10009465),"TA")] = 0;
        us.unification[allele(range(0,10009463,10009465),"T")] = 1;
        us.alleles.push_back(unified_allele(range(0,10009463,10009465),"TA"));
        us.alleles.push_back(unified_allele(range(0,10009463,10009465),"T"));
        us.alleles[1].frequency = 0.1;
        sites.push_back(us);

        genotyper_config cfg;
        cfg.required_dp = 9;
        s = svc->genotype_sites(cfg, string("NA12878D_HiSeqX.21.10009462-10009469"), sites, tfn);
        REQUIRE(s.ok());
    }
}

TEST_CASE("genotype residuals") {
    unique_ptr<VCFData> data;
    Status s = VCFData::Open({"discover_alleles_trio1.vcf", "discover_alleles_trio2.vcf"}, data);
    REQUIRE(s.ok());
    unique_ptr<Service> svc;
    s = Service::Start(service_config(), *data, *data, svc);
    REQUIRE(s.ok());

    discovered_alleles als;
    unsigned N;
    s = svc->discover_alleles("<ALL>", range(0, 0, 1000000), N, als);
    REQUIRE(s.ok());

    vector<unified_site> sites;
    unifier_stats stats;
    s = unified_sites(unifier_config(), N, als, sites, stats);
    REQUIRE(s.ok());

    const string tfn("/tmp/GLnexus_unit_tests.bcf");
    genotyper_config cfg;
    cfg.output_residuals = true;
    s = svc->genotype_sites(cfg, string("<ALL>"), sites, tfn);
    REQUIRE(s.ok());

    // verify that the residuals file is in correct yaml format
    YAML::Node resFile = YAML::LoadFile("/tmp/GLnexus_unit_tests.residuals.yml");
    REQUIRE(resFile.size() == 3);

    // It seems that a list of documents split by --- symbols
    // are parsed as a yaml map.
    REQUIRE(resFile.IsMap());
}

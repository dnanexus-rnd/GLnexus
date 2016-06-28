#include <iostream>
#include "unifier.h"
#include "types.h"
#include "catch.hpp"
using namespace std;
using namespace GLnexus;

TEST_CASE("unifier max_alleles_per_site") {
    discovered_alleles dal;
    discovered_allele_info dai;

    dai.is_ref = true; dai.maxAQ = 99; dai.zGQ = zygosity_by_GQ(1,0,100);
    dal[allele(range(0, 100, 101), "A")] = dai;
    dai.is_ref = false; dai.maxAQ = 99; dai.zGQ = zygosity_by_GQ(1,0,100);
    dal[allele(range(0, 100, 101), "G")] = dai;
    dai.is_ref = false; dai.maxAQ = 99; dai.zGQ = zygosity_by_GQ(1,0,10);
    dal[allele(range(0, 100, 101), "C")] = dai;
    dai.is_ref = false; dai.maxAQ = 99; dai.zGQ = zygosity_by_GQ(1,0,10);
    dal[allele(range(0, 100, 101), "T")] = dai;
    dai.is_ref = true; dai.maxAQ = 99; dai.zGQ = zygosity_by_GQ(1,0,100);
    dal[allele(range(0, 100, 102), "AA")] = dai;
    dai.is_ref = false; dai.maxAQ = 99; dai.zGQ = zygosity_by_GQ(1,0,100);
    dal[allele(range(0, 100, 102), "GG")] = dai;
    dai.is_ref = false; dai.maxAQ = 99; dai.zGQ = zygosity_by_GQ(1,0,1);
    dal[allele(range(0, 100, 102), "GC")] = dai;

    unifier_config cfg;
    vector<unified_site> sites;

    SECTION("0") {
        // No limit on # of max alleles
        Status s = unified_sites(cfg, dal, sites);
        REQUIRE(s.ok());
        REQUIRE(sites.size() == 1);
        REQUIRE(sites[0].alleles.size() == 6);
        REQUIRE(sites[0].alleles[5] == "GC");
    }

    SECTION("5") {
        cfg.max_alleles_per_site = 5;
        Status s = unified_sites(cfg, dal, sites);
        REQUIRE(s.ok());
        REQUIRE(sites.size() == 1);
        REQUIRE(sites[0].alleles.size() == 5);
        REQUIRE(sites[0].alleles[4] == "TA");
    }

    SECTION("4") {
        cfg.max_alleles_per_site = 4;
        Status s = unified_sites(cfg, dal, sites);
        REQUIRE(s.ok());
        REQUIRE(sites.size() == 1);
        REQUIRE(sites[0].alleles.size() == 3);
        REQUIRE(sites[0].alleles[2] == "GG");
        // both CA and TA were pruned as they have the same observation count
    }

    SECTION("3") {
        cfg.max_alleles_per_site = 3;
        Status s = unified_sites(cfg, dal, sites);
        REQUIRE(s.ok());
        REQUIRE(sites.size() == 1);
        REQUIRE(sites[0].alleles.size() == 3);
        REQUIRE(sites[0].alleles[2] == "GG");
    }

    SECTION("2") {
        cfg.max_alleles_per_site = 2;
        Status s = unified_sites(cfg, dal, sites);
        REQUIRE(s.ok());
        REQUIRE(sites.size() == 1);
        REQUIRE(sites[0].alleles.size() == 2);
        REQUIRE(sites[0].alleles[1] == "GA");
        // GA is kept even though it has the same observation count as GG
        // because we have to keep at least one alt allele...
    }
}

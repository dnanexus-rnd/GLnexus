#include <iostream>
#include "unifier.h"
#include "types.h"
#include "catch.hpp"
using namespace std;
using namespace GLnexus;

TEST_CASE("unifier max_alleles_per_site") {
    discovered_alleles dal;

    dal[allele(range(0, 100, 101), "A")] = discovered_allele_info({true, 99, zygosity_by_GQ(1,0,100)});
    dal[allele(range(0, 100, 101), "G")] = discovered_allele_info({false, 99, zygosity_by_GQ(1,0,100)});
    dal[allele(range(0, 100, 101), "C")] = discovered_allele_info({false, 99, zygosity_by_GQ(1,0,10)});
    dal[allele(range(0, 100, 101), "T")] = discovered_allele_info({false, 99, zygosity_by_GQ(1,0,10)});
    dal[allele(range(0, 100, 102), "AA")] = discovered_allele_info({true, 99, zygosity_by_GQ(1,0,100)});
    dal[allele(range(0, 100, 102), "GG")] = discovered_allele_info({false, 99, zygosity_by_GQ(1,0,100)});
    dal[allele(range(0, 100, 102), "GC")] = discovered_allele_info({false, 99, zygosity_by_GQ(1,0,1)});

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

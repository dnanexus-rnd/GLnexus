#include <iostream>
#include <memory>
#include "diploid.h"
#include "catch.hpp"
using namespace std;

TEST_CASE("diploid::genotypes") {
    REQUIRE(diploid::genotypes(1) == 1);
    REQUIRE(diploid::genotypes(2) == 3);
    REQUIRE(diploid::genotypes(3) == 6);
}

TEST_CASE("diploid::gt_alleles") {
    for (int n_allele=0; n_allele<16; n_allele++) {
        int x=0;
        for (int i=0; i<n_allele; i++) {
            for (int j=i; j<n_allele; j++) {
                auto p = diploid::gt_alleles(n_allele, x++);
                REQUIRE(p.first == i);
                REQUIRE(p.second == j);
            }
        }
    }
}

TEST_CASE("diploid::alleles_gt") {
    for (int n_allele=0; n_allele<16; n_allele++) {
        int x=0;
        for (int i=0; i<n_allele; i++) {
            for (int j=i; j<n_allele; j++) {
                REQUIRE(diploid::alleles_gt(n_allele, i, j) == x);
                REQUIRE(diploid::alleles_gt(n_allele, j, i) == x);
                x++;
            }
        }
    }
}

TEST_CASE("diploid::trio::mendelian_inconsistencies") {
    using namespace diploid;

    REQUIRE(trio::mendelian_inconsistencies(2, alleles_gt(2, 0, 0), alleles_gt(2, 0, 0), alleles_gt(2, 0, 0)) == 0);
    REQUIRE(trio::mendelian_inconsistencies(2, alleles_gt(2, 0, 0), alleles_gt(2, 0, 0), alleles_gt(2, 0, 1)) == 1);
    REQUIRE(trio::mendelian_inconsistencies(2, alleles_gt(2, 0, 0), alleles_gt(2, 0, 0), alleles_gt(2, 1, 1)) == 2);
    REQUIRE(trio::mendelian_inconsistencies(2, alleles_gt(2, 0, 0), alleles_gt(2, 0, 1), alleles_gt(2, 0, 1)) == 0);
    REQUIRE(trio::mendelian_inconsistencies(2, alleles_gt(2, 0, 0), alleles_gt(2, 0, 1), alleles_gt(2, 0, 1)) == 0);

    REQUIRE(trio::mendelian_inconsistencies(2, alleles_gt(2, 0, 0), alleles_gt(2, 0, 1), alleles_gt(2, 1, 1)) == 1);
    REQUIRE(trio::mendelian_inconsistencies(2, alleles_gt(2, 0, 0), alleles_gt(2, 1, 1), alleles_gt(2, 1, 1)) == 1);
    REQUIRE(trio::mendelian_inconsistencies(2, alleles_gt(2, 1, 0), alleles_gt(2, 0, 1), alleles_gt(2, 0, 1)) == 0);
    REQUIRE(trio::mendelian_inconsistencies(2, alleles_gt(2, 1, 1), alleles_gt(2, 0, 1), alleles_gt(2, 0, 1)) == 0);

    REQUIRE(trio::mendelian_inconsistencies(3, alleles_gt(3, 0, 1), alleles_gt(3, 2, 2), alleles_gt(3, 2, 2)) == 1);
    REQUIRE(trio::mendelian_inconsistencies(3, alleles_gt(3, 2, 2), alleles_gt(3, 2, 2), alleles_gt(3, 0, 1)) == 2);
    REQUIRE(trio::mendelian_inconsistencies(3, alleles_gt(3, 0, 0), alleles_gt(3, 1, 2), alleles_gt(3, 0, 2)) == 0);
    REQUIRE(trio::mendelian_inconsistencies(3, alleles_gt(3, 0, 0), alleles_gt(3, 1, 2), alleles_gt(3, 1, 2)) == 1);
    REQUIRE(trio::mendelian_inconsistencies(3, alleles_gt(3, 0, 0), alleles_gt(3, 1, 2), alleles_gt(3, 2, 2)) == 1);
    REQUIRE(trio::mendelian_inconsistencies(3, alleles_gt(3, 1, 2), alleles_gt(3, 0, 0), alleles_gt(3, 2, 2)) == 1);

    REQUIRE(trio::mendelian_inconsistencies(3, alleles_gt(3, 0, 0), alleles_gt(3, 1, 2), alleles_gt(3, 0, 2)) == 0);
    REQUIRE(trio::mendelian_inconsistencies(3, alleles_gt(3, 1, 2), alleles_gt(3, 0, 0), alleles_gt(3, 0, 2)) == 0);
    REQUIRE(trio::mendelian_inconsistencies(3, alleles_gt(3, 1, 2), alleles_gt(3, 0, 0), alleles_gt(3, 1, 2)) == 1);
    REQUIRE(trio::mendelian_inconsistencies(3, alleles_gt(3, 1, 2), alleles_gt(3, 0, 0), alleles_gt(3, 2, 2)) == 1);
}

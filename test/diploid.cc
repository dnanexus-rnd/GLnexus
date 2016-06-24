#include <iostream>
#include <memory>
#include "diploid.h"
#include "catch.hpp"
using namespace std;
using namespace GLnexus;

TEST_CASE("diploid::genotypes") {
    REQUIRE(diploid::genotypes(1) == 1);
    REQUIRE(diploid::genotypes(2) == 3);
    REQUIRE(diploid::genotypes(3) == 6);
}

TEST_CASE("diploid::gt_alleles") {
    for (int n_allele=0; n_allele<16; n_allele++) {
        int x=0;
        for (int j=0; j<n_allele; j++) {
            for (int i=0; i<=j; i++) {
                auto p = diploid::gt_alleles(x++);
                REQUIRE(p.first == i);
                REQUIRE(p.second == j);
            }
        }
    }
}

TEST_CASE("diploid::alleles_gt") {
    for (int n_allele=0; n_allele<16; n_allele++) {
        int x=0;
        for (int j=0; j<n_allele; j++) {
            for (int i=0; i<=j; i++) {
                REQUIRE(diploid::alleles_gt(i, j) == x);
                REQUIRE(diploid::alleles_gt(j, i) == x);
                x++;
            }
        }
    }
}

TEST_CASE("diploid::alleles_maxAQ") {
    // tests with one overwhelmingly likely genotype
    for (int n_allele=2; n_allele<16; n_allele++) {
        auto nGT = diploid::genotypes(n_allele);
        for (int gt=0; gt<nGT; gt++) {
            auto al1 = diploid::gt_alleles(gt).first;
            auto al2 = diploid::gt_alleles(gt).second;

            vector<double> gl(diploid::genotypes(n_allele), 0.01);
            gl[gt] = 1.0;
            vector<int> AQ;
            Status s = diploid::alleles_maxAQ(n_allele, 1, {0}, gl, AQ);
            REQUIRE(AQ.size() == n_allele);
            for (int al = 0; al<n_allele; al++) {
                if (al == al1 || al == al2) {
                    REQUIRE(AQ[al] == 20);
                } else {
                    REQUIRE(AQ[al] == 0);
                }
            }
        }
    }

    // tests with a second-rank genotype
    for (int n_allele=2; n_allele<16; n_allele++) {
        auto nGT = diploid::genotypes(n_allele);
        for (int gt=0; gt<nGT; gt++) {
            auto al1 = diploid::gt_alleles(gt).first;
            auto al2 = diploid::gt_alleles(gt).second;

            for (int gt2=0; gt2<nGT; gt2++) {
                if (gt != gt2) {
                    auto al3 = diploid::gt_alleles(gt2).first;
                    auto al4 = diploid::gt_alleles(gt2).second;

                    vector<double> gl(diploid::genotypes(n_allele), 0.001);
                    gl[gt] = 1.0;
                    gl[gt2] = 0.1;
                    vector<int> AQ;
                    Status s = diploid::alleles_maxAQ(n_allele, 1, {0}, gl, AQ);
                    REQUIRE(AQ.size() == n_allele);
                    for (int al = 0; al<n_allele; al++) {
                         if (al == al1 || al == al2) {
                             if (al == al3 || al == al4) {
                                 REQUIRE(AQ[al] == 30);
                             } else {
                                 REQUIRE(AQ[al] == 10);
                             }
                        } else {
                            REQUIRE(AQ[al] == 0);
                        }
                    }
                }
            }
        }
    }

    // TODO: multi-sample tests
}

TEST_CASE("diploid::trio::mendelian_inconsistencies") {
    using namespace diploid;

    REQUIRE(trio::mendelian_inconsistencies(alleles_gt(0, 0), alleles_gt(0, 0), alleles_gt(0, 0)) == 0);
    REQUIRE(trio::mendelian_inconsistencies(alleles_gt(0, 0), alleles_gt(0, 0), alleles_gt(0, 1)) == 1);
    REQUIRE(trio::mendelian_inconsistencies(alleles_gt(0, 0), alleles_gt(0, 0), alleles_gt(1, 1)) == 2);
    REQUIRE(trio::mendelian_inconsistencies(alleles_gt(0, 0), alleles_gt(0, 1), alleles_gt(0, 1)) == 0);
    REQUIRE(trio::mendelian_inconsistencies(alleles_gt(0, 0), alleles_gt(0, 1), alleles_gt(0, 1)) == 0);

    REQUIRE(trio::mendelian_inconsistencies(alleles_gt(0, 0), alleles_gt(0, 1), alleles_gt(1, 1)) == 1);
    REQUIRE(trio::mendelian_inconsistencies(alleles_gt(0, 0), alleles_gt(1, 1), alleles_gt(1, 1)) == 1);
    REQUIRE(trio::mendelian_inconsistencies(alleles_gt(1, 0), alleles_gt(0, 1), alleles_gt(0, 1)) == 0);
    REQUIRE(trio::mendelian_inconsistencies(alleles_gt(1, 1), alleles_gt(0, 1), alleles_gt(0, 1)) == 0);

    REQUIRE(trio::mendelian_inconsistencies(alleles_gt(0, 1), alleles_gt(2, 2), alleles_gt(2, 2)) == 1);
    REQUIRE(trio::mendelian_inconsistencies(alleles_gt(2, 2), alleles_gt(2, 2), alleles_gt(0, 1)) == 2);
    REQUIRE(trio::mendelian_inconsistencies(alleles_gt(0, 0), alleles_gt(1, 2), alleles_gt(0, 2)) == 0);
    REQUIRE(trio::mendelian_inconsistencies(alleles_gt(0, 0), alleles_gt(1, 2), alleles_gt(1, 2)) == 1);
    REQUIRE(trio::mendelian_inconsistencies(alleles_gt(0, 0), alleles_gt(1, 2), alleles_gt(2, 2)) == 1);
    REQUIRE(trio::mendelian_inconsistencies(alleles_gt(1, 2), alleles_gt(0, 0), alleles_gt(2, 2)) == 1);

    REQUIRE(trio::mendelian_inconsistencies(alleles_gt(0, 0), alleles_gt(1, 2), alleles_gt(0, 2)) == 0);
    REQUIRE(trio::mendelian_inconsistencies(alleles_gt(1, 2), alleles_gt(0, 0), alleles_gt(0, 2)) == 0);
    REQUIRE(trio::mendelian_inconsistencies(alleles_gt(1, 2), alleles_gt(0, 0), alleles_gt(1, 2)) == 1);
    REQUIRE(trio::mendelian_inconsistencies(alleles_gt(1, 2), alleles_gt(0, 0), alleles_gt(2, 2)) == 1);
}

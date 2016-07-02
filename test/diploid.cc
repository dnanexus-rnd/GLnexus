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

TEST_CASE("diploid::alleles_topAQ") {
    // tests with one overwhelmingly likely genotype
    for (int n_allele=2; n_allele<16; n_allele++) {
        auto nGT = diploid::genotypes(n_allele);
        for (int gt=0; gt<nGT; gt++) {
            auto al1 = diploid::gt_alleles(gt).first;
            auto al2 = diploid::gt_alleles(gt).second;

            vector<double> gll(diploid::genotypes(n_allele), log(0.01));
            gll[gt] = log(1.0);
            vector<top_AQ> AQ;
            Status s = diploid::alleles_topAQ(n_allele, 1, {0}, gll, AQ);
            REQUIRE(s.ok());
            REQUIRE(AQ.size() == n_allele);
            for (int al = 0; al<n_allele; al++) {
                if (al == al1 || al == al2) {
                    REQUIRE(AQ[al].V[0] == 20);
                } else {
                    REQUIRE(AQ[al].V[0] == 0);
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

                    vector<double> gll(diploid::genotypes(n_allele), log(0.001));
                    gll[gt] = log(1.0);
                    gll[gt2] = log(0.1);
                    vector<top_AQ> AQ;
                    Status s = diploid::alleles_topAQ(n_allele, 1, {0}, gll, AQ);
                    REQUIRE(s.ok());
                    REQUIRE(AQ.size() == n_allele);
                    for (int al = 0; al<n_allele; al++) {
                         if (al == al1 || al == al2) {
                             if (al == al3 || al == al4) {
                                 REQUIRE(AQ[al].V[0] == 30);
                             } else {
                                 REQUIRE(AQ[al].V[0] == 10);
                             }
                        } else {
                            REQUIRE(AQ[al].V[0] == 0);
                        }
                    }
                }
            }
        }
    }

    SECTION("top_AQ ranking") {
        top_AQ topAQ, topAQ2;
        auto &V = topAQ.V, &V2 = topAQ2.V;
        vector<int> v {1, 5, 4};
        topAQ += v;
        REQUIRE(V[0] == 5);
        REQUIRE(V[1] == 4);
        REQUIRE(V[2] == 1);
        REQUIRE(V[3] == 0);
        v = {2, 8, 3}; topAQ += v;
        REQUIRE(V[0] == 8);
        REQUIRE(V[1] == 5);
        REQUIRE(V[2] == 4);
        REQUIRE(V[3] == 3);
        REQUIRE(V[4] == 2);
        REQUIRE(V[5] == 1);
        REQUIRE(V[6] == 0);
        v = {8, 2, 3}; topAQ += v;
        REQUIRE(V[0] == 8);
        REQUIRE(V[1] == 8);
        REQUIRE(V[2] == 5);
        REQUIRE(V[3] == 4);
        REQUIRE(V[4] == 3);
        REQUIRE(V[5] == 3);
        REQUIRE(V[6] == 2);
        REQUIRE(V[7] == 2);
        REQUIRE(V[8] == 1);
        REQUIRE(V[9] == 0);
        v = {6, 7, 9, 0}; topAQ += v;
        REQUIRE(V[0] == 9);
        REQUIRE(V[1] == 8);
        REQUIRE(V[2] == 8);
        REQUIRE(V[3] == 7);
        REQUIRE(V[4] == 6);
        REQUIRE(V[5] == 5);
        REQUIRE(V[6] == 4);
        REQUIRE(V[7] == 3);
        REQUIRE(V[8] == 3);
        REQUIRE(V[9] == 2);
        topAQ2 = topAQ;
        topAQ2 += topAQ;
        REQUIRE(V2[0] == 9);
        REQUIRE(V2[1] == 9);
        REQUIRE(V2[2] == 8);
        REQUIRE(V2[3] == 8);
        REQUIRE(V2[4] == 8);
        REQUIRE(V2[5] == 8);
        REQUIRE(V2[6] == 7);
        REQUIRE(V2[7] == 7);
        REQUIRE(V2[8] == 6);
        REQUIRE(V2[9] == 6);
    }

    SECTION("potential phred underflow") {
        vector<int32_t> pl {4364,0,4983};
        vector<double> gll;
        for (auto pl_i : pl) {
            gll.push_back(double(pl_i)/(-10.0)/log10(exp(1.0)));
        }
        vector<top_AQ> AQ;
        Status s = diploid::alleles_topAQ(2, 1, {0}, gll, AQ);
        REQUIRE(s.ok());
        REQUIRE(AQ[0].V[0] == 4983);
        REQUIRE(AQ[1].V[0] == 4364);
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

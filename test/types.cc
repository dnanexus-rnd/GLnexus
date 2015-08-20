#include <iostream>
#include "types.h"
#include "catch.hpp"
using namespace std;
using namespace GLnexus;

TEST_CASE("range::overlaps") {
    range r(0, 1000, 1100);
    REQUIRE(r.overlaps(range(0, 1000, 1100)));
    REQUIRE(r.overlaps(range(0, 1001, 1100)));
    REQUIRE(r.overlaps(range(0, 1000, 1099)));
    REQUIRE(r.overlaps(range(0, 1050, 1090)));
    REQUIRE(r.overlaps(range(0, 999, 1001)));
    REQUIRE_FALSE(r.overlaps(range(0, 666, 999)));
    REQUIRE_FALSE(r.overlaps(range(0, 1101, 1200)));
    REQUIRE_FALSE(r.overlaps(range(1, 1000, 1100)));

    REQUIRE_FALSE(r.overlaps(range(0, 999, 1000)));
    REQUIRE(r.overlaps(range(0, 999, 1001)));

    REQUIRE_FALSE(r.overlaps(range(0, 1100, 1101)));
    REQUIRE(r.overlaps(range(0, 1099, 1100)));
}

TEST_CASE("range sort order") {
    REQUIRE(range(0, 1000, 1100) < range(1, 1000, 1100));
    REQUIRE(range(0, 1000, 1100) < range(1, 999, 1100));
    REQUIRE(range(0, 1000, 1100) < range(0, 1001, 1100));
    REQUIRE(range(0, 1000, 1100) < range(0, 1000, 1101));
    REQUIRE_FALSE(range(0, 1000, 1100) < range(0, 1000, 1100));
    REQUIRE_FALSE(range(0, 1000, 1100) < range(0, 999, 1100));
    REQUIRE_FALSE(range(1, 1000, 1100) < range(0, 999, 1100));
    REQUIRE_FALSE(range(1, 999, 1100) < range(0, 1000, 1100));
    REQUIRE_FALSE(range(0, 1000, 1099) < range(0, 999, 1100));

    vector<range> v { range(1, 1000, 1100), range(0, 1000, 1100), range(0, 1000, 1100), range(1, 999, 1100), range(0, 1001, 1100), range(1, 999, 1100) };
    sort(v.begin(), v.end());

    for (unsigned i = 0; i < v.size()-1; i++) {
        REQUIRE(v[i] <= v[i+1]);
    }
}

TEST_CASE("range_of_bcf") {
    #define UPD(T,name,ini,del) std::unique_ptr<T, void(*)(T*)> up_##name((ini), (del)); auto name = up_##name.get();
    UPD(vcfFile, vcf, bcf_open("test/data/NA12878D_HiSeqX.21.10009462-10009469.gvcf", "r"), [](vcfFile* f) { bcf_close(f); });
    UPD(bcf_hdr_t, hdr, bcf_hdr_read(vcf), &bcf_hdr_destroy);
    shared_ptr<bcf1_t> vt;
    vector<shared_ptr<bcf1_t>> records;

    do {
        if (vt) {
            REQUIRE(bcf_unpack(vt.get(), BCF_UN_ALL) == 0);
            records.push_back(vt);
        }
        vt = shared_ptr<bcf1_t>(bcf_init(), &bcf_destroy);
    } while (bcf_read(vcf, hdr, vt.get()) == 0);

    REQUIRE(records.size() == 5);

    range rng;
    Status s = range_of_bcf(hdr, records[0], rng);
    REQUIRE(s.ok());
    REQUIRE(rng.rid == 0);
    REQUIRE(rng.beg == 10009461);
    REQUIRE(rng.end == 10009463);
    REQUIRE(rng.size() == 2);

    s = range_of_bcf(hdr, records[1], rng);
    REQUIRE(s.ok());
    REQUIRE(rng.rid == 0);
    REQUIRE(rng.beg == 10009463);
    REQUIRE(rng.end == 10009465);
    REQUIRE(rng.size() == 2);

    s = range_of_bcf(hdr, records[2], rng);
    REQUIRE(s.ok());
    REQUIRE(rng.rid == 0);
    REQUIRE(rng.beg == 10009465);
    REQUIRE(rng.end == 10009466);
    REQUIRE(rng.size() == 1);
}

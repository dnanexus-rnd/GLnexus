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

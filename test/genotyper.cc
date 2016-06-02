#include <iostream>
#include "genotyper.h"
#include "types.h"
#include "catch.hpp"
using namespace std;
using namespace GLnexus;

TEST_CASE("One_call_ordering") {

    SECTION("Numerical calls") {
        one_call call1=one_call(1, NoCallReason::N_A);
        one_call call2=one_call(2, NoCallReason::N_A);
        REQUIRE(call2 > call1);
    }

    SECTION("No call"){
        one_call call1 = one_call();
        one_call call2 = one_call(0, NoCallReason::N_A);
        REQUIRE(call2 < call1);
    }

}

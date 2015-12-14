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

TEST_CASE("LossTracker_full_overlap") {
    unified_site site = unified_site(range(0, 1, 1000));
    LossTracker loss_tracker = LossTracker(site.pos);

    loss_tracker.add_call_for_site(range(0, 10, 15), 2, false);
    loss_tracker.finalize_loss_for_site(0);

    loss_stats loss;
    Status s = loss_tracker.get(loss);

    REQUIRE(s.ok());
    REQUIRE(loss.n_calls_total == 2);
    REQUIRE(loss.n_bp_total == 10);
    REQUIRE(loss.n_no_calls_total == 0);
    REQUIRE(loss.n_bp_lost == 0);
    REQUIRE(loss.n_calls_lost == 0);
}

TEST_CASE("LossTracker_invalid_operations") {
    unified_site site = unified_site(range(0, 1, 1000));
    LossTracker loss_tracker = LossTracker(site.pos);

    loss_tracker.add_call_for_site(range(0, 10, 15), 2, false);

    loss_stats loss;
    Status s = loss_tracker.get(loss);
    REQUIRE(!s.ok());

    loss_tracker.finalize_loss_for_site(0);

    s = loss_tracker.add_call_for_site(range(0, 20, 25), 2, false);
    REQUIRE(!s.ok());

    s = loss_tracker.finalize_loss_for_site(2);
    REQUIRE(!s.ok());

}


TEST_CASE("LossTracker_partial_overlap") {
    unified_site site = unified_site(range(0, 1, 1000));
    LossTracker loss_tracker = LossTracker(site.pos);

    loss_tracker.add_call_for_site(range(0, 995, 1005), 2, false);
    loss_tracker.finalize_loss_for_site(0);

    loss_stats loss;
    Status s = loss_tracker.get(loss);

    REQUIRE(s.ok());
    REQUIRE(loss.n_calls_total == 2);
    REQUIRE(loss.n_bp_total == 10);
    REQUIRE(loss.n_no_calls_total == 0);
    REQUIRE(loss.n_bp_lost == 0);
    REQUIRE(loss.n_calls_lost == 0);
}


TEST_CASE("LossTracker_1_missing_call") {
    unified_site site = unified_site(range(0, 1, 1000));
    LossTracker loss_tracker = LossTracker(site.pos);

    loss_tracker.add_call_for_site(range(0, 995, 1005), 2, false);
    loss_tracker.finalize_loss_for_site(1);

    loss_stats loss;
    Status s = loss_tracker.get(loss);

    REQUIRE(s.ok());
    REQUIRE(loss.n_calls_total == 2);
    REQUIRE(loss.n_bp_total == 10);
    REQUIRE(loss.n_no_calls_total == 1);
    REQUIRE(loss.n_bp_lost == 5);
    REQUIRE(loss.n_calls_lost == 1);
}

TEST_CASE("LossTracker_2_missing_call") {
    unified_site site = unified_site(range(0, 1, 1000));
    LossTracker loss_tracker = LossTracker(site.pos);

    loss_tracker.add_call_for_site(range(0, 995, 1005), 2, false);
    loss_tracker.finalize_loss_for_site(2);

    loss_stats loss;
    Status s = loss_tracker.get(loss);

    REQUIRE(s.ok());
    REQUIRE(loss.n_calls_total == 2);
    REQUIRE(loss.n_bp_total == 10);
    REQUIRE(loss.n_no_calls_total == 2);
    REQUIRE(loss.n_bp_lost == 10);
    REQUIRE(loss.n_calls_lost == 2);
}


TEST_CASE("LossTracker_1_missing_call_on_1_orig_call") {
    unified_site site = unified_site(range(0, 1, 1000));
    LossTracker loss_tracker = LossTracker(site.pos);

    loss_tracker.add_call_for_site(range(0, 995, 1005), 1, false);
    loss_tracker.finalize_loss_for_site(1);

    loss_stats loss;
    Status s = loss_tracker.get(loss);

    REQUIRE(s.ok());
    REQUIRE(loss.n_calls_total == 1);
    REQUIRE(loss.n_bp_total == 5);
    REQUIRE(loss.n_no_calls_total == 1);
    REQUIRE(loss.n_bp_lost == 0);
    REQUIRE(loss.n_calls_lost == 0);
}


TEST_CASE("LossTracker_2_missing_call_on_1_orig_call") {
    unified_site site = unified_site(range(0, 1, 1000));
    LossTracker loss_tracker = LossTracker(site.pos);

    loss_tracker.add_call_for_site(range(0, 995, 1005), 1, false);
    loss_tracker.finalize_loss_for_site(2);

    loss_stats loss;
    Status s = loss_tracker.get(loss);

    REQUIRE(s.ok());
    REQUIRE(loss.n_calls_total == 1);
    REQUIRE(loss.n_bp_total == 5);
    REQUIRE(loss.n_no_calls_total == 2);
    REQUIRE(loss.n_bp_lost == 5);
    REQUIRE(loss.n_calls_lost == 1);
}

TEST_CASE("loss_summary_add_stats") {
    consolidated_loss summary;

    unified_site site1 = unified_site(range(0, 1, 1000));
    LossTracker loss_tracker = LossTracker(site1.pos);
    loss_tracker.add_call_for_site(range(0, 995, 1005), 2, false);
    loss_tracker.finalize_loss_for_site(2);

    loss_stats loss1;
    Status s = loss_tracker.get(loss1);
    REQUIRE(s.ok());

    summary.insert(make_pair("sample_1", loss1));

    consolidated_loss summary2;
    unified_site site2 = unified_site(range(0, 500, 1500));
    LossTracker loss_tracker2 = LossTracker(site2.pos);
    loss_tracker2.add_call_for_site(range(0, 900, 910), 2, false);
    loss_tracker2.finalize_loss_for_site(1);

    loss_stats loss2;
    s = loss_tracker2.get(loss2);
    REQUIRE(s.ok());

    summary2.insert(make_pair("sample_1", loss2));

    s =merge_loss_stats(summary2, summary);
    REQUIRE(s.ok());

    loss_stats updated_loss = summary.find("sample_1")->second;
    REQUIRE(updated_loss.n_calls_total == 4);
    REQUIRE(updated_loss.n_calls_lost == 3);
    REQUIRE(updated_loss.n_no_calls_total == 3);
    REQUIRE(updated_loss.n_bp_lost == 20);
    REQUIRE(updated_loss.n_bp_total == 30);
}
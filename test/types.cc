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

    range rng(records[0]);
    REQUIRE(rng.rid == 0);
    REQUIRE(rng.beg == 10009461);
    REQUIRE(rng.end == 10009463);
    REQUIRE(rng.size() == 2);

    rng = range(records[1]);
    REQUIRE(rng.rid == 0);
    REQUIRE(rng.beg == 10009463);
    REQUIRE(rng.end == 10009465);
    REQUIRE(rng.size() == 2);

    rng = range(records[2]);
    REQUIRE(rng.rid == 0);
    REQUIRE(rng.beg == 10009465);
    REQUIRE(rng.end == 10009466);
    REQUIRE(rng.size() == 1);
}

TEST_CASE("loss_stats_full_overlap") {
    unified_site site = unified_site(range(0, 1, 1000));
    loss_stats loss = loss_stats(site);

    loss.add_call_for_site(range(0, 10, 15), 2);
    REQUIRE(loss.get_n_calls_total() == 2);
    REQUIRE(loss.get_n_bp_total() == 10);

    loss.add_call_for_site(range(0, 20, 25), 1);
    REQUIRE(loss.get_n_calls_total() == 3);
    REQUIRE(loss.get_n_bp_total() == 15);

    loss.add_call_for_site(range(0, 30, 35), 0);
    REQUIRE(loss.get_n_calls_total() == 3);
    REQUIRE(loss.get_n_bp_total() == 15);
}

TEST_CASE("loss_stats_partial_overlap") {
    unified_site site = unified_site(range(0, 1, 1000));
    loss_stats loss = loss_stats(site);

    loss.add_call_for_site(range(0, 995, 1005), 2);
    REQUIRE(loss.get_n_calls_total() == 2);
    REQUIRE(loss.get_n_bp_total() == 10);
}

TEST_CASE("loss_stats_no_missing_calls") {
    unified_site site = unified_site(range(0, 1, 1000));
    loss_stats loss = loss_stats(site);

    loss.add_call_for_site(range(0, 995, 1005), 2);
    loss.finalize_loss_for_site(0);

    REQUIRE(loss.get_n_no_calls_total() == 0);
    REQUIRE(loss.get_n_bp_lost() == 0);
    REQUIRE(loss.get_n_calls_lost() == 0);
}

TEST_CASE("loss_stats_1_missing_call") {
    unified_site site = unified_site(range(0, 1, 1000));
    loss_stats loss = loss_stats(site);

    loss.add_call_for_site(range(0, 995, 1005), 2);
    loss.finalize_loss_for_site(1);

    REQUIRE(loss.get_n_no_calls_total() == 1);
    REQUIRE(loss.get_n_bp_lost() == 5);
    REQUIRE(loss.get_n_calls_lost() == 1);
}

TEST_CASE("loss_stats_2_missing_call") {
    unified_site site = unified_site(range(0, 1, 1000));
    loss_stats loss = loss_stats(site);

    loss.add_call_for_site(range(0, 995, 1005), 2);
    loss.finalize_loss_for_site(2);

    REQUIRE(loss.get_n_no_calls_total() == 2);
    REQUIRE(loss.get_n_bp_lost() == 10);
    REQUIRE(loss.get_n_calls_lost() == 2);
}

TEST_CASE("loss_stats_1_missing_call_on_1_orig_call") {
    unified_site site = unified_site(range(0, 1, 1000));
    loss_stats loss = loss_stats( site);

    loss.add_call_for_site(range(0, 995, 1005), 1);
    loss.finalize_loss_for_site(1);

    REQUIRE(loss.get_n_no_calls_total() == 1);
    REQUIRE(loss.get_n_bp_lost() == 0);
    REQUIRE(loss.get_n_calls_lost() == 0);
}

TEST_CASE("loss_stats_2_missing_call_on_1_orig_call") {
    unified_site site = unified_site(range(0, 1, 1000));
    loss_stats loss = loss_stats(site);

    loss.add_call_for_site(range(0, 995, 1005), 1);
    loss.finalize_loss_for_site(2);

    REQUIRE(loss.get_n_no_calls_total() == 2);
    REQUIRE(loss.get_n_bp_lost() == 5);
    REQUIRE(loss.get_n_calls_lost() == 1);
}

TEST_CASE("loss_summary_add_stats") {
    consolidated_loss summary;

    unified_site site1 = unified_site(range(0, 1, 1000));
    loss_stats loss = loss_stats(site1);
    loss.add_call_for_site(range(0, 995, 1005), 2);
    loss.finalize_loss_for_site(2);

    summary.insert(make_pair("sample_1", loss));

    unified_site site2 = unified_site(range(0, 500, 1500));
    loss_stats loss2 = loss_stats(site2);
    loss2.add_call_for_site(range(0, 900, 910), 2);
    loss2.finalize_loss_for_site(1);

    consolidated_loss summary2;
    summary2.insert(make_pair("sample_1", loss2));

    Status s =merge_loss_stats(summary2, summary);
    REQUIRE(s.ok());

    loss_stats updated_loss = summary.find("sample_1")->second;
    REQUIRE(updated_loss.get_n_calls_total() == 4);
    REQUIRE(updated_loss.get_n_calls_lost() == 3);
    REQUIRE(updated_loss.get_n_no_calls_total() == 3);
    REQUIRE(updated_loss.get_n_bp_lost() == 20);
    REQUIRE(updated_loss.get_n_bp_total() == 30);
}
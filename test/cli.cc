#include <iostream>
#include <fstream>
#include <map>
#include "spdlog/spdlog.h"
#include "BCFKeyValueData.h"
#include "BCFSerialize.h"
#include "cli_utils.h"
#include "catch.hpp"

using namespace std;
using namespace GLnexus;

static auto console = spdlog::stderr_logger_mt("GLnexus_cli");

TEST_CASE("cli") {
    Status s;

    // setup database directory
    string dbdir = "/tmp/cli";
    REQUIRE(system(("rm -rf " + dbdir).c_str()) == 0);
    REQUIRE(system(("mkdir -p " + dbdir).c_str()) == 0);
    string dbpath = dbdir + "/DB";

    string basedir = "test/data/cli";
    string exemplar_gvcf = basedir + "/" + "F1.gvcf";
    vector<pair<string,size_t>> contigs;
    s = cli::utils::db_init(console, dbpath, exemplar_gvcf, contigs);
    REQUIRE(s.ok());
    REQUIRE(contigs.size() >= 1);

    vector<string> gvcfs;
    for (auto fname : {"F1.gvcf", "F2.gvcf"}) {
         gvcfs.push_back(basedir + "/" + fname);
    }
    vector<range> ranges;
    s = cli::utils::db_bulk_load(console, gvcfs, dbpath, ranges, 8, contigs);
    REQUIRE(s.ok());
    REQUIRE(contigs.size() >= 1);

    // construct ranges that cover the entire contig set
    ranges.clear();
    for (int rid=0; rid < contigs.size(); rid++) {
        size_t len = contigs[rid].second;
        ranges.push_back(range(rid, 0, len));
    }

    discovered_alleles dsals;
    unsigned sample_count;
    s = cli::utils::discover_alleles(console, dbpath, ranges, contigs, 8, dsals, sample_count);
    REQUIRE(s.ok());

    string filename = dbdir + "/dsals.yml";
    s = cli::utils::yaml_write_discovered_alleles_to_file(dsals, contigs, sample_count, filename);
    REQUIRE(s.ok());
}

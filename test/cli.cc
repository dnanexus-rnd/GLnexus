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

static auto console = spdlog::stderr_logger_mt("cli_test");
static int nr_threads = 2;
static string DB_DIR = "/tmp/cli";
static string DB_PATH = DB_DIR + "/DB";

TEST_CASE("cli") {
    SECTION("Everything") {
        Status s;

        // setup database directory
        REQUIRE(system(("rm -rf " + DB_DIR).c_str()) == 0);
        REQUIRE(system(("mkdir -p " + DB_DIR).c_str()) == 0);

        string basedir = "test/data/cli";
        string exemplar_gvcf = basedir + "/" + "F1.gvcf.gz";
        vector<pair<string,size_t>> contigs;
        s = cli::utils::db_init(console, DB_PATH, exemplar_gvcf, contigs);
        REQUIRE(s.ok());
        REQUIRE(contigs.size() >= 1);

        // DB directory already exists -- make sure we get an error
        s = cli::utils::db_init(console, DB_PATH, exemplar_gvcf, contigs);
        REQUIRE(s.bad());

        vector<string> gvcfs;
        for (auto fname : {"F1.gvcf.gz", "F2.gvcf.gz", "F3.gvcf.gz", "F4.gvcf.gz"}) {
            gvcfs.push_back(basedir + "/" + fname);
        }

        // empty range, so that all records will be loaded into the database
        vector<range> ranges;
        s = cli::utils::db_bulk_load(console, nr_threads, gvcfs, DB_PATH, ranges, contigs);
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
        s = cli::utils::discover_alleles(console, nr_threads, DB_PATH, ranges, contigs, dsals, sample_count);
        REQUIRE(s.ok());

        string filename = DB_DIR + "/dsals.yml";
        s = cli::utils::yaml_write_discovered_alleles_to_file(dsals, contigs, sample_count, filename);
        REQUIRE(s.ok());

        GLnexus::unifier_config unifier_cfg;
        GLnexus::genotyper_config genotyper_cfg;
        string config_preset = "test";
        s = cli::utils::load_config_preset(console, config_preset, unifier_cfg, genotyper_cfg);
        REQUIRE(s.ok());
        unifier_cfg.min_AQ1 = 0;
        unifier_cfg.min_AQ2 = 0;

        // unify sites
        vector<GLnexus::unified_site> sites;
        unifier_stats stats;
        s = cli::utils::unify_sites(console, unifier_cfg,
                                    contigs, dsals, sample_count, sites, stats);
        filename = DB_DIR + "/sites.yml";
        s = cli::utils::write_unified_sites_to_file(sites, contigs, filename);
        REQUIRE(s.ok());

        filename = DB_DIR + "/results.bcf";
        s = cli::utils::genotype(console, nr_threads, DB_PATH, genotyper_cfg, sites, {}, filename);
        REQUIRE(s.ok());
    }

    SECTION("read contigs") {
        Status s;

        vector<pair<string,size_t>> contigs;
        s = cli::utils::db_get_contigs(console, DB_PATH, contigs);
        REQUIRE(s.ok());

        // check that the contigs are correct
        // ##contig=<ID=1,length=3000>
        // ##contig=<ID=2,length=4000>
        // ##contig=<ID=3,length=1000>
        // ##contig=<ID=4,length=4000>
        REQUIRE(contigs[0].first == "1");
        REQUIRE(contigs[0].second == 3000);
        REQUIRE(contigs[1].first == "2");
        REQUIRE(contigs[1].second == 4000);
        REQUIRE(contigs[2].first == "3");
        REQUIRE(contigs[2].second == 1000);
        REQUIRE(contigs[3].first == "4");
        REQUIRE(contigs[3].second == 4000);
    }
}

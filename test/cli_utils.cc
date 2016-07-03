#include <iostream>
#include <fstream>
#include <cstdio>
#include <memory>
#include "cli_utils.h"
#include "catch.hpp"

using namespace std;
using namespace GLnexus;
using namespace GLnexus::cli;

        // create two yaml files in the right format
static void yaml_file_of_discover_alleles(const string &filename,
                                          const vector<pair<string,size_t>> &contigs,
                                          const vector<range> ranges,
                                          const char* buf) {
    std::remove(filename.c_str());

    vector<discovered_alleles> valleles;
    YAML::Node n = YAML::Load(buf);
    discovered_alleles dal1;
    Status s = discovered_alleles_of_yaml(n, contigs, dal1);
    REQUIRE(s.ok());
    valleles.push_back(dal1);

    YAML::Emitter yaml;
    s = utils::yaml_of_contigs_alleles_ranges(contigs, ranges, valleles, yaml);
    REQUIRE(s.ok());

    ofstream fos;
    fos.open(filename);
    fos << yaml.c_str();
    fos.close();
}

static auto console = spdlog::stderr_logger_mt("cli_utils_test");

TEST_CASE("cli_utils") {
    vector<pair<string,size_t>> contigs;
    contigs.push_back(make_pair("16",12345));
    contigs.push_back(make_pair("17",23456));

    const char* da_yaml1 = 1 + R"(
- range: {ref: '16', beg: 100, end: 100}
  dna: A
  is_ref: true
  top_AQ: [99]
  zygosity_by_GQ: [[100,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
- range: {ref: '16', beg: 113, end: 120}
  dna: G
  is_ref: false
  top_AQ: [99]
  zygosity_by_GQ: [[0,0],[10,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,5]]
)";

    const char* da_yaml2 = 1 + R"(
- range: {ref: '17', beg: 100, end: 100}
  dna: A
  is_ref: true
  top_AQ: [99]
  zygosity_by_GQ: [[100,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
- range: {ref: '17', beg: 200, end: 310}
  dna: G
  is_ref: false
  top_AQ: [99]
  zygosity_by_GQ: [[0,0],[10,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,5]]
)";

    const char* da_yaml3 = 1 + R"(
- range: {ref: '16', beg: 107, end: 109}
  dna: A
  is_ref: true
  maxAQ: 99
  zygosity_by_GQ: [[100,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
- range: {ref: '17', beg: 220, end: 330}
  dna: G
  is_ref: true
  maxAQ: 99
  zygosity_by_GQ: [[0,0],[10,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,5]]
)";

    SECTION("parse_range") {
        GLnexus::range query(-1,-1,-1);
        string range_txt = "17:100-2000";
        REQUIRE(utils::parse_range(contigs, range_txt, query));

        range_txt = "20:10000-30000";
        REQUIRE(!utils::parse_range(contigs, range_txt, query));

        range_txt = "17:0-10000000";
        REQUIRE(!utils::parse_range(contigs, range_txt, query));

        string cmdline("xxx yyy");
        vector<range> ranges;
        REQUIRE(!utils::parse_ranges(contigs, cmdline, ranges));

        cmdline = string("16:10-37,17:901-3220");
        REQUIRE(utils::parse_ranges(contigs, cmdline, ranges));
    }

    SECTION("yaml_discovered_alleles") {
        vector<range> ranges;
        ranges.push_back(range(0, 1000, 1100));
        ranges.push_back(range(1, 40, 70));

        vector<discovered_alleles> valleles;
        {
            YAML::Node n = YAML::Load(da_yaml1);
            discovered_alleles dal1;
            Status s = discovered_alleles_of_yaml(n, contigs, dal1);
            REQUIRE(s.ok());
            valleles.push_back(dal1);

            n = YAML::Load(da_yaml2);
            discovered_alleles dal2;
            s = discovered_alleles_of_yaml(n, contigs, dal2);
            REQUIRE(s.ok());
            valleles.push_back(dal2);
        }

        YAML::Emitter yaml;
        Status s = utils::yaml_of_contigs_alleles_ranges(contigs, ranges, valleles, yaml);
        REQUIRE(s.ok());

        YAML::Node n = YAML::Load(yaml.c_str());
        std::vector<std::pair<std::string,size_t> > contigs2;
        std::vector<range> ranges2;
        std::vector<discovered_alleles> valleles2;

        s = utils::contigs_alleles_ranges_of_yaml(n, contigs2, ranges2, valleles2);
        REQUIRE(valleles.size() == valleles2.size());
        for (int i=0; i < valleles.size(); i++) {
            REQUIRE(valleles[i] == valleles2[i]);
        }
    }

    SECTION("yaml_discovered_alleles, #ranges does not match #valleles") {
        vector<range> ranges;
        ranges.push_back(range(0, 1000, 1100));
        ranges.push_back(range(1, 40, 70));

        vector<discovered_alleles> valleles;
        YAML::Emitter yaml;
        Status s = utils::yaml_of_contigs_alleles_ranges(contigs, ranges, valleles, yaml);
        REQUIRE(s.bad());
    }


    const char* snp = 1 + R"(
range: {ref: '17', beg: 100, end: 100}
alleles: [A, G]
copy_number: [100, 51]
unification:
  - range: {ref: '17', beg: 100, end: 100}
    alt: A
    to: 0
  - range: {ref: '17', beg: 100, end: 100}
    alt: G
    to: 1
)";


    const char* del = 1 + R"(
range: {ref: '17', beg: 1000, end: 1001}
alleles: [AG, AC, C]
copy_number: [100, 50, 1]
unification:
  - range: {ref: '17', beg: 1000, end: 1001}
    alt: AG
    to: 0
  - range: {ref: '17', beg: 1000, end: 1001}
    alt: AC
    to: 1
  - range: {ref: '17', beg: 1000, end: 1001}
    alt: C
    to: 2
  - range: {ref: '17', beg: 1001, end: 1001}
    alt: C
    to: 1
)";

    SECTION("yaml_of_unified_sites") {
        // generate a vector of sites
        vector<unified_site> sites;
        {
            YAML::Node n = YAML::Load(snp);
            unified_site us(range(-1,-1,-1));
            Status s = unified_site::of_yaml(n, contigs, us);
            REQUIRE(s.ok());
            sites.push_back(us);

            n = YAML::Load(del);
            s = unified_site::of_yaml(n, contigs, us);
            REQUIRE(s.ok());
            sites.push_back(us);
        }

        // convert to yaml
        YAML::Emitter yaml;
        {
            Status s = utils::yaml_of_unified_sites(sites, contigs, yaml);
            REQUIRE(s.ok());
        }

        // convert back and compare
        {
            YAML::Node n = YAML::Load(yaml.c_str());
            vector<unified_site> sites2;
            Status s = utils::unified_sites_of_yaml(n, contigs, sites2);
            REQUIRE(s.ok());
            REQUIRE(sites.size() == sites2.size());
            for (int i=0; i < sites.size(); i++) {
                REQUIRE(sites[i] == sites2[i]);
            }
        }
    }

    SECTION("LoadYAMLFile") {
        string tmp_file_name = "/tmp/xxx.yml";
        std::remove(tmp_file_name.c_str());

        {
            YAML::Node node;
            string emptyname = "";
            Status s = utils::LoadYAMLFile(emptyname, node);
            REQUIRE(s.bad());
        }

        // Create a trivial YAML file
        {
            YAML::Emitter yaml;

            yaml << YAML::BeginSeq;
            yaml << "x";
            yaml << "y";
            yaml << "z";
            yaml << YAML::EndSeq;

            ofstream fos;
            fos.open(tmp_file_name);
            fos << yaml.c_str();
            fos.close();
        }

        // verify the file
        {
            YAML::Node node;
            Status s = utils::LoadYAMLFile(tmp_file_name, node);
            REQUIRE(s.ok());
        }

        // Create an invalid YAML file
        std::remove(tmp_file_name.c_str());
        {
            ofstream fos;
            fos.open(tmp_file_name);
            fos << "hello: false:" << endl;
            fos << "world: oyster" << endl;
            fos << "movies: Dune Airplane Princess Bride" << endl;
            fos.close();
        }

        // verify the file does not load
        {
            YAML::Node node;
            Status s = utils::LoadYAMLFile(tmp_file_name, node);
            REQUIRE(s.bad());
        }

        // Create an invalid YAML file
        std::remove(tmp_file_name.c_str());
        {
            ofstream fos;
            fos.open(tmp_file_name);
            fos << "" << endl;
            fos.close();
        }

        // verify the file does not load
        {
            YAML::Node node;
            Status s = utils::LoadYAMLFile(tmp_file_name, node);
            REQUIRE(s.bad());
        }

    }


    SECTION("merging discovered allele files, special case, only one file") {
        vector<range> ranges;
        ranges.push_back(range(0, 1, 1100));

        string tmp_file_name1 = "/tmp/xxx_1.yml";
        yaml_file_of_discover_alleles(tmp_file_name1, contigs, ranges, da_yaml1);

        vector<string> filenames;
        filenames.push_back(tmp_file_name1);

        vector<pair<string,size_t>> contigs2;
        vector<range> ranges2;
        vector<discovered_alleles> valleles2;
        Status s = utils::merge_discovered_allele_files(console, filenames, contigs2, ranges2, valleles2);
        REQUIRE(s.ok());
        REQUIRE(contigs2.size() == contigs.size());
    }

    SECTION("merging discovered allele files, 2 files") {
        // create two files, with different ranges
        string tmp_file_name1 = "/tmp/xxx_1.yml";
        {
            vector<range> ranges;
            ranges.push_back(range(0, 1, 1100));
            yaml_file_of_discover_alleles(tmp_file_name1, contigs, ranges, da_yaml1);
        }

        string tmp_file_name2 = "/tmp/xxx_2.yml";
        {
            vector<range> ranges2;
            ranges2.push_back(range(0, 2002, 2208));
            yaml_file_of_discover_alleles(tmp_file_name2, contigs, ranges2, da_yaml2);
        }

        vector<string> filenames;
        filenames.push_back(tmp_file_name1);
        filenames.push_back(tmp_file_name2);

        vector<pair<string,size_t>> contigs2;
        vector<range> ranges2;
        vector<discovered_alleles> valleles;
        Status s = utils::merge_discovered_allele_files(console, filenames, contigs2, ranges2, valleles);
        REQUIRE(s.ok());
        REQUIRE(contigs2.size() == contigs.size());
        REQUIRE(valleles.size() == 2);
        REQUIRE(valleles[0].size() == 2);
        REQUIRE(valleles[1].size() == 2);
    }

    SECTION("merging discovered allele files, 3 files") {
        vector<range> ranges;
        ranges.push_back(range(0, 1, 1100));

        vector<string> filenames;
        const char* yamls[] = {da_yaml1, da_yaml2, da_yaml3};
        const char* i_filenames[3] = {"/tmp/xxx_1.yml", "/tmp/xxx_2.yml", "/tmp/xxx_3.yml"};
        for (int i=0; i < 3; i++) {
            string fname = string(i_filenames[i]);
            yaml_file_of_discover_alleles(fname, contigs, ranges, yamls[i]);
            filenames.push_back(fname);
        }

        vector<pair<string,size_t>> contigs2;
        vector<range> ranges2;
        vector<discovered_alleles> valleles;
        Status s = utils::merge_discovered_allele_files(console, filenames, contigs2, ranges2, valleles);
        REQUIRE(s.ok());
        REQUIRE(contigs2.size() == contigs.size());
        REQUIRE(valleles.size() == 1);
        REQUIRE(valleles[0].size() == 6);
    }

    SECTION("merging discovered allele files, error, no files provided") {
        vector<range> ranges;
        ranges.push_back(range(0, 1, 1100));
        vector<string> filenames;
        vector<discovered_alleles> valleles;

        Status s = utils::merge_discovered_allele_files(console, filenames, contigs, ranges, valleles);
    }

}

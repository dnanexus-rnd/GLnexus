#include <iostream>
#include <memory>
#include "applet_utils.h"
#include "catch.hpp"

using namespace std;
using namespace GLnexus;
using namespace GLnexus::applet;

TEST_CASE("applet_utils") {
    vector<pair<string,size_t>> contigs;
    contigs.push_back(make_pair("16",12345));
    contigs.push_back(make_pair("17",23456));

    const char* da_yaml1 = 1 + R"(
- range: {ref: '16', beg: 100, end: 100}
  dna: A
  is_ref: true
  copy_number: 100
- range: {ref: '16', beg: 113, end: 120}
  dna: G
  is_ref: false
  copy_number: 10.5
)";

    const char* da_yaml2 = 1 + R"(
- range: {ref: '17', beg: 100, end: 100}
  dna: A
  is_ref: true
  copy_number: 100
- range: {ref: '17', beg: 200, end: 310}
  dna: G
  is_ref: false
  copy_number: 10.5
)";

    SECTION("parse_range") {
        GLnexus::range query(-1,-1,-1);
        string range_txt = "17:100-2000";
        REQUIRE(utils::parse_range(contigs, range_txt, query));

        range_txt = "20:10000-30000";
        REQUIRE(!utils::parse_range(contigs, range_txt, query));
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
        cout << s.str() << endl;
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
}

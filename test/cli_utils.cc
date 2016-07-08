#include <iostream>
#include <fstream>
#include <cstdio>
#include <memory>
#include "cli_utils.h"
#include "catch.hpp"

using namespace std;
using namespace GLnexus;
using namespace GLnexus::cli;

// create a yaml file in the right format
static void yaml_file_of_discovered_alleles(const string &filename,
                                            const vector<pair<string,size_t>> &contigs,
                                            const char* buf) {
    std::remove(filename.c_str());
    stringstream iss(buf, ios_base::in);

    YAML::Node node = YAML::Load(buf);
    discovered_alleles dals;
    Status s = discovered_alleles_of_yaml(node, contigs, dals);
    REQUIRE(s.ok());

    ofstream fos;
    fos.open(filename);
    s = utils::yaml_stream_of_discovered_alleles(contigs, dals, fos);
    REQUIRE(s.ok());
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
  top_AQ: [99]
  zygosity_by_GQ: [[100,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
- range: {ref: '17', beg: 220, end: 330}
  dna: G
  is_ref: true
  top_AQ: [99]
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

    SECTION("yaml_discovered_alleles, an empty list") {
        // an empty list
        std::stringstream ss;
        {
            GLnexus::discovered_alleles dsals;
            Status s = utils::yaml_stream_of_discovered_alleles(contigs, dsals, ss);
            REQUIRE(s.ok());
        }
        discovered_alleles dal_empty;
        vector<pair<string,size_t>> contigs2;
        Status s = utils::discovered_alleles_of_yaml_stream(ss, contigs2, dal_empty);
        REQUIRE(s.bad());
        REQUIRE(contigs == contigs2);
    }

    SECTION("yaml_discovered_alleles, an empty contigs") {
        Status s;

        //YAML::Node n = YAML::Load(da_yaml1);
        discovered_alleles dal1;
        // s = discovered_alleles_of_yaml(n, contigs, dal1);
        //REQUIRE(s.ok());

        vector<pair<string,size_t>> contigs_empty;
        std::stringstream ss;
        s = utils::yaml_stream_of_discovered_alleles(contigs_empty, dal1, ss);
        REQUIRE(s.ok());

        discovered_alleles dal_empty;
        vector<pair<string,size_t>> contigs2;
        s = utils::discovered_alleles_of_yaml_stream(ss, contigs2, dal_empty);
        REQUIRE(s.bad());
    }

    SECTION("yaml_discovered_alleles") {
        discovered_alleles dsals;
        {
            YAML::Node n = YAML::Load(da_yaml1);
            discovered_alleles dal1;
            Status s = discovered_alleles_of_yaml(n, contigs, dal1);
            REQUIRE(s.ok());
            merge_discovered_alleles(dal1, dsals);

            n = YAML::Load(da_yaml2);
            discovered_alleles dal2;
            s = discovered_alleles_of_yaml(n, contigs, dal2);
            REQUIRE(s.ok());
            merge_discovered_alleles(dal2, dsals);
        }

        std::stringstream ss;
        Status s = utils::yaml_stream_of_discovered_alleles(contigs, dsals, ss);
        REQUIRE(s.ok());

        std::vector<std::pair<std::string,size_t> > contigs2;
        discovered_alleles dsals2;
        s = utils::discovered_alleles_of_yaml_stream(ss, contigs2, dsals2);
        REQUIRE(contigs == contigs2);
        REQUIRE(dsals.size() == dsals2.size());
    }

    const char* bad_yaml_1 = 1 + R"(
contigs: xxx
alleles: yyy
)";

    const char* bad_yaml_2 = 1 + R"(
contigs:
  - name: A1
    size: 1000
  - name: A2
    size: 300000
alleles: yyy
)";


    SECTION("bad yaml inputs for discovered_alleles_of_yaml_stream") {
        discovered_alleles dsals;
        {
            stringstream ss("aaa");
            Status s = utils::discovered_alleles_of_yaml_stream(ss, contigs, dsals);
            REQUIRE(s.bad());
        }

        {
            stringstream ss(bad_yaml_1);
            Status s = utils::discovered_alleles_of_yaml_stream(ss, contigs, dsals);
            REQUIRE(s.bad());
        }

        {
            stringstream ss(bad_yaml_2);
            Status s = utils::discovered_alleles_of_yaml_stream(ss, contigs, dsals);
            REQUIRE(s.bad());
        }

        {
            std::stringstream ss;
            ss.write("---\n", 4);
            ss.write("xxx\n", 4);
            ss.write("yyy\n", 4);
            ss.write("zzz\n", 4);
            ss.write("----\n", 5);
            ss.write("zzz\n", 4);
            ss.write("...\n", 5);

            discovered_alleles dal_empty;
            vector<pair<string,size_t>> contigs2;
            Status s = utils::discovered_alleles_of_yaml_stream(ss, contigs2, dal_empty);
            REQUIRE(s.bad());
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
        string tmp_file_name1 = "/tmp/xxx_1.yml";
        yaml_file_of_discovered_alleles(tmp_file_name1, contigs, da_yaml1);

        vector<string> filenames;
        filenames.push_back(tmp_file_name1);

        vector<pair<string,size_t>> contigs2;
        discovered_alleles dsals2;
        Status s = utils::merge_discovered_allele_files(console, filenames, contigs2, dsals2);
        REQUIRE(s.ok());
        REQUIRE(contigs2.size() == contigs.size());
    }

    SECTION("merging discovered allele files, 2 files") {
        // create two files, with different ranges
        string tmp_file_name1 = "/tmp/xxx_1.yml";
        yaml_file_of_discovered_alleles(tmp_file_name1, contigs, da_yaml1);

        string tmp_file_name2 = "/tmp/xxx_2.yml";
        yaml_file_of_discovered_alleles(tmp_file_name2, contigs, da_yaml2);

        vector<string> filenames;
        filenames.push_back(tmp_file_name1);
        filenames.push_back(tmp_file_name2);

        vector<pair<string,size_t>> contigs2;
        discovered_alleles dsals;
        Status s = utils::merge_discovered_allele_files(console, filenames, contigs2, dsals);
        REQUIRE(s.ok());
        REQUIRE(contigs2.size() == contigs.size());
        REQUIRE(dsals.size() == 4);
    }

    SECTION("merging discovered allele files, 3 files") {
        vector<string> filenames;
        const char* yamls[] = {da_yaml1, da_yaml2, da_yaml3};
        const char* i_filenames[3] = {"/tmp/xxx_1.yml", "/tmp/xxx_2.yml", "/tmp/xxx_3.yml"};
        for (int i=0; i < 3; i++) {
            string fname = string(i_filenames[i]);
            yaml_file_of_discovered_alleles(fname, contigs, yamls[i]);
            filenames.push_back(fname);
        }

        vector<pair<string,size_t>> contigs2;
        discovered_alleles dsals;
        Status s = utils::merge_discovered_allele_files(console, filenames, contigs2, dsals);
        REQUIRE(s.ok());
        REQUIRE(contigs2.size() == contigs.size());
        REQUIRE(dsals.size() == 6);
    }

    SECTION("merging discovered allele files, error, no files provided") {
        vector<string> filenames;
        discovered_alleles dsals;

        Status s = utils::merge_discovered_allele_files(console, filenames, contigs, dsals);
        REQUIRE(s.bad());
    }

    SECTION("merging discovered allele files, error 2") {
        // create two files, with different ranges
        string tmp_file_name1 = "/tmp/xxx_1.yml";
        yaml_file_of_discovered_alleles(tmp_file_name1, contigs, da_yaml1);

        vector<pair<string,size_t>> contigs2;
        contigs2.push_back(make_pair("16",550000));
        contigs2.push_back(make_pair("17",8811991));
        string tmp_file_name2 = "/tmp/xxx_2.yml";
        yaml_file_of_discovered_alleles(tmp_file_name2, contigs2, da_yaml2);

        vector<string> filenames;
        filenames.push_back(tmp_file_name1);
        filenames.push_back(tmp_file_name2);

        vector<pair<string,size_t>> contigs3;
        discovered_alleles dsals;
        Status s = utils::merge_discovered_allele_files(console, filenames, contigs3, dsals);
        REQUIRE(s.bad());
    }

    const char* da_yaml10 = 1 + R"(
- range: {ref: '16', beg: 107, end: 109}
  dna: A
  is_ref: true
  top_AQ: [99]
  zygosity_by_GQ: [[100,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
)";

    const char* da_yaml11 = 1 + R"(
- range: {ref: '16', beg: 107, end: 109}
  dna: G
  is_ref: true
  top_AQ: [99]
  zygosity_by_GQ: [[100,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
- range: {ref: '17', beg: 500, end: 501}
  dna: G
  is_ref: true
  top_AQ: [99]
  zygosity_by_GQ: [[0,0],[10,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,5]]
- range: {ref: '17', beg: 1190, end: 1200}
  dna: G
  is_ref: true
  top_AQ: [99]
  zygosity_by_GQ: [[0,0],[10,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,5]]
)";

    SECTION("merging discovered allele files, error, number of sites doesn't match") {
        vector<string> filenames;
        const char* yamls[] = {da_yaml10, da_yaml11};
        const char* i_filenames[3] = {"/tmp/xxx_1.yml", "/tmp/xxx_2.yml"};
        for (int i=0; i < 2; i++) {
            string fname = string(i_filenames[i]);
            yaml_file_of_discovered_alleles(fname, contigs, yamls[i]);
            filenames.push_back(fname);
        }

        discovered_alleles dsals;
        Status s = utils::merge_discovered_allele_files(console, filenames, contigs, dsals);
    }
}

TEST_CASE("find_containing_range") {
    vector<pair<string,size_t>> contigs;
    contigs.push_back(make_pair("1",100000));
    contigs.push_back(make_pair("2",130001));
    contigs.push_back(make_pair("3",24006));

    SECTION("5 ranges") {
        std::set<range> ranges;
        ranges.insert(range(1,15,20));
        ranges.insert(range(1,108,151));
        ranges.insert(range(2,31,50));
        ranges.insert(range(2,42,48));
        ranges.insert(range(3,1000,1507));

        Status s;
        range ans(-1,-1,-1);

        // beginning
        range pos(1,17,18);
        s = utils::find_containing_range(ranges, pos, ans);
        REQUIRE(s.ok());
        REQUIRE(ans == range(1,15,20));

        // last
        pos = range(3,1200,1218);
        s = utils::find_containing_range(ranges, pos, ans);
        REQUIRE(s.ok());
        REQUIRE(ans == range(3,1000,1507));

        pos = range(3,1000,1000);
        s = utils::find_containing_range(ranges, pos, ans);
        REQUIRE(s.ok());
        REQUIRE(ans == range(3,1000,1507));

        // middle
        pos = range(2,31,35);
        s = utils::find_containing_range(ranges, pos, ans);
        REQUIRE(s.ok());
        REQUIRE(ans == range(2,31,50));

        // missing range
        pos = range(2,301,302);
        s = utils::find_containing_range(ranges, pos, ans);
        REQUIRE(s.bad());

        // missing range
        pos = range(3,2000,2003);
        s = utils::find_containing_range(ranges, pos, ans);
        REQUIRE(s.bad());

        // missing range
        pos = range(1,2,4);
        s = utils::find_containing_range(ranges, pos, ans);
        REQUIRE(s.bad());
    }

    SECTION("1 range") {
        std::set<range> ranges;
        ranges.insert(range(1,15,20));

        Status s;
        range ans(-1,-1,-1);

        // beginning
        range pos(1,17,18);
        s = utils::find_containing_range(ranges, pos, ans);
        REQUIRE(s.ok());
        REQUIRE(ans == range(1,15,20));

        // missing range
        pos = range(3,1200,1218);
        s = utils::find_containing_range(ranges, pos, ans);
        REQUIRE(s.bad());

        // missing range
        pos = range(1,2,5);
        s = utils::find_containing_range(ranges, pos, ans);
        REQUIRE(s.bad());

        // non overlapping
        pos = range(1,15,30);
        s = utils::find_containing_range(ranges, pos, ans);
        REQUIRE(s.bad());
    }

    SECTION("no ranges") {
        Status s;
        range ans(-1,-1,-1);
        std::set<range> ranges;

        range pos(1,17,18);
        s = utils::find_containing_range(ranges, pos, ans);
        REQUIRE(s.bad());
    }

    SECTION("larger example") {
        vector<pair<string,size_t>> contigs;
        contigs.push_back(make_pair("1",260000000));
        contigs.push_back(make_pair("10",140000000));
        contigs.push_back(make_pair("11",100000000));

        set<range> ranges;
        ranges.insert(range(1,30366,30503));
        ranges.insert(range(1,35720,35737));
        ranges.insert(range(1,69090,70009));
        ranges.insert(range(1,249152026,249152059));
        ranges.insert(range(1,249208061,249208079));
        ranges.insert(range(1,249210800,249214146));
        ranges.insert(range(2,92996,94055));
        ranges.insert(range(2,95121,95179));
        ranges.insert(range(2,255828,255989));
        ranges.insert(range(2,135480471,135481678));
        ranges.insert(range(2,135491589,135491842));
        ranges.insert(range(3,168957,169053));
        ranges.insert(range(3,180208,180405));
        ranges.insert(range(3,5685264,5685404));
        ranges.insert(range(3,5700990,5701408));
        ranges.insert(range(3,5717462,5717886));

        range pos(1, 249211350, 249211350);
        range ans(-1,-1,-1);
        Status s = utils::find_containing_range(ranges,pos,ans);
        REQUIRE(s.ok());
        REQUIRE(ans == range(1,249210800,249214146));
    }
}

#include <iostream>
#include <fstream>
#include <cstdio>
#include <memory>
#include "cli_utils.h"
#include "catch.hpp"

using namespace std;
using namespace GLnexus;
using namespace GLnexus::cli;

static auto console = spdlog::stderr_logger_mt("cli_utils_test");

// create a yaml file in the right format
static void capnp_file_of_discovered_alleles(const string &filename,
                                             unsigned int sample_count,
                                             const vector<pair<string,size_t>> &contigs,
                                             const char* buf) {
    YAML::Node node = YAML::Load(buf);
    discovered_alleles dals;
    Status s = discovered_alleles_of_yaml(node, contigs, dals);
    REQUIRE(s.ok());

    // A sanity check
    REQUIRE(GLnexus::capnp_discover_alleles_verify(sample_count, contigs, dals, filename).ok());
    std::remove(filename.c_str());

    s = capnp_of_discovered_alleles(sample_count, contigs, dals, filename);
    REQUIRE(s.ok());
}

TEST_CASE("cli_utils") {
    unsigned N;
    vector<pair<string,size_t>> contigs;
    contigs.push_back(make_pair("16",12345));
    contigs.push_back(make_pair("17",23456));

    const char* da_yaml1 = 1 + R"(
- range: {ref: '16', beg: 100, end: 100}
  dna: A
  is_ref: true
  all_filtered: false
  top_AQ: [99]
  zygosity_by_GQ: [[100,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
- range: {ref: '16', beg: 113, end: 120}
  dna: G
  is_ref: false
  all_filtered: false
  top_AQ: [99]
  zygosity_by_GQ: [[0,0],[10,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,5]]
)";

    const char* da_yaml2 = 1 + R"(
- range: {ref: '17', beg: 100, end: 100}
  dna: A
  is_ref: true
  all_filtered: false
  top_AQ: [99]
  zygosity_by_GQ: [[100,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
- range: {ref: '17', beg: 200, end: 310}
  dna: G
  is_ref: false
  all_filtered: false
  top_AQ: [99]
  zygosity_by_GQ: [[0,0],[10,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,5]]
)";

    const char* da_yaml3 = 1 + R"(
- range: {ref: '16', beg: 107, end: 109}
  dna: A
  is_ref: true
  all_filtered: false
  top_AQ: [99]
  zygosity_by_GQ: [[100,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
- range: {ref: '17', beg: 220, end: 330}
  dna: G
  is_ref: true
  all_filtered: false
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
            Status s = utils::yaml_stream_of_discovered_alleles(1, contigs, dsals, ss);
            REQUIRE(s.ok());
        }
        discovered_alleles dal_empty;
        vector<pair<string,size_t>> contigs2;
        Status s = utils::discovered_alleles_of_yaml_stream(ss, N, contigs2, dal_empty);
        REQUIRE(s.ok());
        REQUIRE(N == 1);
        REQUIRE(contigs == contigs2);
    }

    SECTION("yaml_discovered_alleles, an empty contigs") {
        Status s;

        discovered_alleles dal1;
        vector<pair<string,size_t>> contigs_empty;
        std::stringstream ss;
        s = utils::yaml_stream_of_discovered_alleles(0, contigs_empty, dal1, ss);
        REQUIRE(s.ok());

        discovered_alleles dal_empty;
        vector<pair<string,size_t>> contigs2;
        s = utils::discovered_alleles_of_yaml_stream(ss, N, contigs2, dal_empty);
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
        Status s = utils::yaml_stream_of_discovered_alleles(1, contigs, dsals, ss);
        REQUIRE(s.ok());

        std::vector<std::pair<std::string,size_t> > contigs2;
        discovered_alleles dsals2;
        s = utils::discovered_alleles_of_yaml_stream(ss, N, contigs2, dsals2);
        REQUIRE(N == 1);
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
            Status s = utils::discovered_alleles_of_yaml_stream(ss, N, contigs, dsals);
            REQUIRE(s.bad());
        }

        {
            stringstream ss(bad_yaml_1);
            Status s = utils::discovered_alleles_of_yaml_stream(ss, N, contigs, dsals);
            REQUIRE(s.bad());
        }

        {
            stringstream ss(bad_yaml_2);
            Status s = utils::discovered_alleles_of_yaml_stream(ss, N, contigs, dsals);
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
            Status s = utils::discovered_alleles_of_yaml_stream(ss, N, contigs2, dal_empty);
            REQUIRE(s.bad());
        }

        // Good contigs, but bad discovered alleles
        {
            std::stringstream ss;
            {
                GLnexus::discovered_alleles dsals;
                Status s = utils::yaml_stream_of_discovered_alleles(1, contigs, dsals, ss);
                REQUIRE(s.ok());
            }

            // Extra characters after EOF, that should not be there
            ss << "ssss" << endl;

            discovered_alleles dal_empty;
            vector<pair<string,size_t>> contigs2;
            Status s = utils::discovered_alleles_of_yaml_stream(ss, N, contigs2, dal_empty);
            REQUIRE(s.bad());
        }
    }


    const char* snp = 1 + R"(
range: {ref: '17', beg: 100, end: 100}
alleles: [A, G]
allele_frequencies: [.nan, 0.51]
quality: 317
unification:
  - range: {ref: '17', beg: 100, end: 100}
    dna: A
    to: 0
  - range: {ref: '17', beg: 100, end: 100}
    dna: G
    to: 1
)";


    const char* del = 1 + R"(
range: {ref: '17', beg: 1000, end: 1001}
alleles: [AG, AC, C]
allele_frequencies: [.nan, 0.50, 0.1]
quality: 317
unification:
  - range: {ref: '17', beg: 1000, end: 1001}
    dna: AG
    to: 0
  - range: {ref: '17', beg: 1000, end: 1001}
    dna: AC
    to: 1
  - range: {ref: '17', beg: 1000, end: 1001}
    dna: C
    to: 2
  - range: {ref: '17', beg: 1001, end: 1001}
    dna: C
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
        stringstream ss;
        {
            Status s = utils::yaml_stream_of_unified_sites(sites, contigs, ss);
            REQUIRE(s.ok());
        }

        // convert back and compare
        {
            vector<unified_site> sites2;
            Status s = utils::unified_sites_of_yaml_stream(ss, contigs, sites2);
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
        capnp_file_of_discovered_alleles(tmp_file_name1, 1, contigs, da_yaml1);

        vector<string> filenames;
        filenames.push_back(tmp_file_name1);

        vector<pair<string,size_t>> contigs2;
        vector<discovered_alleles> dsals2;
        Status s = utils::merge_discovered_allele_files(console, 0, filenames, N, contigs2, dsals2);
        REQUIRE(s.ok());
        REQUIRE(N == 1);
        REQUIRE(contigs2.size() == contigs.size());
        REQUIRE(dsals2.size() == contigs.size());
    }

    SECTION("merging discovered allele files, 2 files") {
        // create two files, with different ranges
        string tmp_file_name1 = "/tmp/xxx_1.yml";
        capnp_file_of_discovered_alleles(tmp_file_name1, 1, contigs, da_yaml1);

        string tmp_file_name2 = "/tmp/xxx_2.yml";
        capnp_file_of_discovered_alleles(tmp_file_name2, 2, contigs, da_yaml2);

        vector<string> filenames;
        filenames.push_back(tmp_file_name1);
        filenames.push_back(tmp_file_name2);

        vector<pair<string,size_t>> contigs2;
        vector<discovered_alleles> dsals;
        Status s = utils::merge_discovered_allele_files(console, 0, filenames, N, contigs2, dsals);
        REQUIRE(s.ok());
        REQUIRE(N == 3);
        REQUIRE(contigs2.size() == contigs.size());
        REQUIRE(dsals.size() == contigs.size());
        REQUIRE(dsals[0].size() == 2);
        REQUIRE(dsals[1].size() == 2);
    }

    SECTION("merging discovered allele files, 3 files") {
        vector<string> filenames;
        const char* yamls[] = {da_yaml1, da_yaml2, da_yaml3};
        const char* i_filenames[3] = {"/tmp/xxx_1.yml", "/tmp/xxx_2.yml", "/tmp/xxx_3.yml"};
        for (int i=0; i < 3; i++) {
            string fname = string(i_filenames[i]);
            capnp_file_of_discovered_alleles(fname, 1, contigs, yamls[i]);
            filenames.push_back(fname);
        }

        vector<pair<string,size_t>> contigs2;
        vector<discovered_alleles> dsals;
        Status s = utils::merge_discovered_allele_files(console, 0, filenames, N, contigs2, dsals);
        REQUIRE(s.ok());
        REQUIRE(N == 3);
        REQUIRE(contigs2.size() == contigs.size());
        REQUIRE(dsals.size() == contigs.size());
        REQUIRE(dsals[0].size() == 3);
        REQUIRE(dsals[1].size() == 3);

        discovered_alleles da1, da2, da3;
        YAML::Node node = YAML::Load(da_yaml1);
        s = discovered_alleles_of_yaml(node, contigs, da1); REQUIRE(s.ok());
        node = YAML::Load(da_yaml2);
        s = discovered_alleles_of_yaml(node, contigs, da2); REQUIRE(s.ok());
        node = YAML::Load(da_yaml3);
        s = discovered_alleles_of_yaml(node, contigs, da3); REQUIRE(s.ok());

        vector<pair<allele,discovered_allele_info>> vdsals(dsals[0].begin(), dsals[0].end());
        vector<pair<allele,discovered_allele_info>> vda(da1.begin(), da1.end());
        REQUIRE(vdsals[0] == vda[0]);
        REQUIRE(vdsals[2] == vda[1]);
        vda.assign(da3.begin(), da3.end());
        REQUIRE(vdsals[1] == vda[0]);

        vdsals.assign(dsals[1].begin(), dsals[1].end());
        REQUIRE(vdsals[2] == vda[1]);
        vda.assign(da2.begin(), da2.end());
        REQUIRE(vdsals[0] == vda[0]);
        REQUIRE(vdsals[1] == vda[1]);
    }

    SECTION("merging discovered allele files, error, no files provided") {
        vector<string> filenames;
        vector<discovered_alleles> dsals;

        Status s = utils::merge_discovered_allele_files(console, 0, filenames, N, contigs, dsals);
        REQUIRE(s.bad());
    }

    SECTION("merging discovered allele files, error 2") {
        // create two files, with different ranges
        string tmp_file_name1 = "/tmp/xxx_1.yml";
        capnp_file_of_discovered_alleles(tmp_file_name1, 1, contigs, da_yaml1);

        vector<pair<string,size_t>> contigs2;
        contigs2.push_back(make_pair("16",550000));
        contigs2.push_back(make_pair("17",8811991));
        string tmp_file_name2 = "/tmp/xxx_2.yml";
        capnp_file_of_discovered_alleles(tmp_file_name2, 1, contigs2, da_yaml2);

        vector<string> filenames;
        filenames.push_back(tmp_file_name1);
        filenames.push_back(tmp_file_name2);

        vector<pair<string,size_t>> contigs3;
        vector<discovered_alleles> dsals;
        Status s = utils::merge_discovered_allele_files(console, 0, filenames, N, contigs3, dsals);
        REQUIRE(s.bad());
    }

    const char* da_yaml10 = 1 + R"(
- range: {ref: '16', beg: 107, end: 109}
  dna: A
  is_ref: true
  all_filtered: false
  top_AQ: [99]
  zygosity_by_GQ: [[100,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
)";

    const char* da_yaml11 = 1 + R"(
- range: {ref: '16', beg: 107, end: 109}
  dna: G
  is_ref: true
  all_filtered: false
  top_AQ: [99]
  zygosity_by_GQ: [[100,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
- range: {ref: '17', beg: 500, end: 501}
  dna: G
  is_ref: true
  all_filtered: false
  top_AQ: [99]
  zygosity_by_GQ: [[0,0],[10,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,5]]
- range: {ref: '17', beg: 1190, end: 1200}
  dna: G
  is_ref: true
  all_filtered: false
  top_AQ: [99]
  zygosity_by_GQ: [[0,0],[10,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,5]]
)";

    SECTION("merging discovered allele files, error, number of sites doesn't match") {
        vector<string> filenames;
        const char* yamls[] = {da_yaml10, da_yaml11};
        const char* i_filenames[3] = {"/tmp/xxx_1.yml", "/tmp/xxx_2.yml"};
        for (int i=0; i < 2; i++) {
            string fname = string(i_filenames[i]);
            capnp_file_of_discovered_alleles(fname, 1, contigs, yamls[i]);
            filenames.push_back(fname);
        }

        vector<discovered_alleles> dsals;
        Status s = utils::merge_discovered_allele_files(console, 0, filenames, N, contigs, dsals);
    }
}

TEST_CASE("find_target_range") {
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
        s = find_target_range(ranges, pos, ans);
        REQUIRE(s.ok());
        REQUIRE(ans == range(1,15,20));

        // last
        pos = range(3,1200,1218);
        s = find_target_range(ranges, pos, ans);
        REQUIRE(s.ok());
        REQUIRE(ans == range(3,1000,1507));

        pos = range(3,1000,1001);
        s = find_target_range(ranges, pos, ans);
        REQUIRE(s.ok());
        REQUIRE(ans == range(3,1000,1507));

        // middle
        pos = range(2,31,35);
        s = find_target_range(ranges, pos, ans);
        REQUIRE(s.ok());
        REQUIRE(ans == range(2,31,50));

        // missing range
        pos = range(2,301,302);
        s = find_target_range(ranges, pos, ans);
        REQUIRE(s.bad());

        // missing range
        pos = range(3,2000,2003);
        s = find_target_range(ranges, pos, ans);
        REQUIRE(s.bad());

        // missing range
        pos = range(1,2,4);
        s = find_target_range(ranges, pos, ans);
        REQUIRE(s.bad());

        // dangling not within (from the left)
        pos = range(1,105,110);
        s = find_target_range(ranges, pos, ans);
        REQUIRE(s.ok());
        REQUIRE(ans == range(1,108,151));

        // dangling to the right
        pos = range(1,150,155);
        s = find_target_range(ranges, pos, ans);
        REQUIRE(s.ok());
        REQUIRE(ans == range(1,108,151));
    }

    SECTION("1 range") {
        std::set<range> ranges;
        ranges.insert(range(1,15,20));

        Status s;
        range ans(-1,-1,-1);

        // beginning
        range pos(1,17,18);
        s = find_target_range(ranges, pos, ans);
        REQUIRE(s.ok());
        REQUIRE(ans == range(1,15,20));

        // missing range
        pos = range(3,1200,1218);
        s = find_target_range(ranges, pos, ans);
        REQUIRE(s.bad());

        // missing range
        pos = range(1,2,5);
        s = find_target_range(ranges, pos, ans);
        REQUIRE(s.bad());

        // dangling
        pos = range(1,15,30);
        s = find_target_range(ranges, pos, ans);
        REQUIRE(s.ok());
    }

    SECTION("no ranges") {
        Status s;
        range ans(-1,-1,-1);
        std::set<range> ranges;

        range pos(1,17,18);
        s = find_target_range(ranges, pos, ans);
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
        Status s = find_target_range(ranges,pos,ans);
        REQUIRE(s.ok());
        REQUIRE(ans == range(1,249210800,249214146));
    }
}

static void make_files_in_dir(const string &basedir) {
    auto names = {"a", "b", "c"};
    for (auto cursor : names) {
        string path = basedir + "/" + cursor;
        int retval = system(("touch " + path).c_str());
        REQUIRE(retval == 0);
        REQUIRE(utils::check_file_exists(path));
    }
}

TEST_CASE("file_ops") {
    string basedir = "/tmp/cli_utils";
    int retval = system(("rm -rf " + basedir).c_str());
    REQUIRE(retval == 0);
    retval = system(("mkdir -p " + basedir).c_str());
    REQUIRE(retval == 0);
    REQUIRE(utils::check_dir_exists(basedir));

    string file_path = basedir + "/foo.txt";
    retval = system(("touch " + file_path).c_str());
    REQUIRE(retval == 0);
    REQUIRE(utils::check_file_exists(file_path));

    make_files_in_dir(basedir);

    // create some subdirectories
    auto subdirs = {"X", "Y", "Z"};
    for (auto cursor : subdirs) {
        string path = basedir + "/" + cursor;
        retval = system(("mkdir " + path).c_str());
        REQUIRE(retval == 0);

        make_files_in_dir(path);
    }

    retval = system(("rm -rf " + basedir).c_str());
    REQUIRE(retval == 0);
    REQUIRE(!utils::check_dir_exists(basedir));
}

TEST_CASE("parse_bed_file") {
    Status s;

    vector<pair<string,size_t>> contigs;
    contigs.push_back(make_pair("1",10000000));
    contigs.push_back(make_pair("2",31000000));

    string basedir = "test/data/cli";
    string bedfilename = basedir + "/vcr_test.bed";
    vector<range> ranges;
    s = utils::parse_bed_file(console, bedfilename, contigs, ranges);
    REQUIRE(s.ok());
    REQUIRE(ranges.size() == 10);
}

TEST_CASE("iter_compare") {
    Status s;

    // setup database directory
    string dbdir = "/tmp/iter_compare";
    REQUIRE(system(("rm -rf " + dbdir).c_str()) == 0);
    REQUIRE(system(("mkdir -p " + dbdir).c_str()) == 0);
    string dbpath = dbdir + "/DB";

    string basedir = "test/data/cli";
    string exemplar_gvcf = basedir + "/" + "F1.gvcf.gz";
    vector<pair<string,size_t>> contigs;
    s = cli::utils::db_init(console, dbpath, exemplar_gvcf, contigs);
    REQUIRE(s.ok());
    REQUIRE(contigs.size() >= 1);

    vector<string> gvcfs;
    for (auto fname : {"F1.gvcf.gz", "F2.gvcf.gz"}) {
         gvcfs.push_back(basedir + "/" + fname);
    }
    vector<range> ranges;
    s = cli::utils::db_bulk_load(console, 8, gvcfs, dbpath, ranges, contigs);
    REQUIRE(s.ok());
    REQUIRE(contigs.size() >= 1);

    int n_iter = 50;
    s = cli::utils::compare_db_itertion_algorithms(console, dbpath, n_iter);
    console->info() << "Passed " << n_iter << " iterator comparison tests";
}

#include <iostream>
#include <fstream>
#include "service.h"
#include "types.h"
#include "utils.cc"
#include "catch.hpp"
#include "sys/stat.h"
#include "yaml-cpp/yaml.h"
using namespace std;
using namespace GLnexus;

#define test_contigs {{"A", 0},}

/// GVCFTestCase is the overarching class for unit-testing
/// of a curated "case" for the glnexus workflow.
/// It handles parsing of a yaml file with information on
/// gvcf input. It supports comparing output of the GLnexus workflow
/// with the expected unified_sites and joint genotype output.
class GVCFTestCase {
    #define V(pred,msg) if (!(pred)) return Status::Invalid("GVCFTestCase::load_yml: ", msg);

    // Root of folder in which GVCFTestCases input files
    const string INPUT_ROOT_DIR = "test/data/gvcf_test_cases/";

    // Name of temp dir to store output
    const string TEMP_DIR = "/tmp/GLnexus/";

    // Whether the following input information is expected
    bool test_dsals;
    bool test_usites;
    bool test_genotypes;

    // Limit verbosity, prevent reprinting of the case header
    bool has_printed_header = false;

    // Format fields to validate
    vector<string> validated_formats;

    // Info fields to validate
    vector<string> validated_infos;

    // Name of curated case
    string name;

    // unifier configuration
    unifier_config unifier_cfg;

    // genotyper configuration
    genotyper_config genotyper_cfg;

    // Filenames of input gvcfs
    set<string> input_gvcfs;

    // Path to temp directory where files are written to
    string temp_dir_path;

    // Filename of truth vcf
    string truth_vcf_path;

    // Header text (for each input gvcf & output vcf)
    string header;

    // Parsed "truth" unfiied sites
    vector<unified_site> truth_sites;

    // Parsed "truth" discovered sites
    discovered_alleles truth_dsals;

    // Unified sites based on unify_sites algorithm
    vector<unified_site> sites;

    // Storage for loaded input gvcfs
    unique_ptr<VCFData> data;

    // Service for algorithmic procedures
    unique_ptr<Service> svc;

    // Contigs
    vector<pair<string,size_t>> contigs;

    // Creates the folder based on given path, do not raise
    // error if folder already exists
    Status create_folder(const string path) {
        int status = mkdir(path.c_str(), S_IRWXU);

        // Failed to create folder and it doesn't already exist
        if (status && errno != EEXIST) {
            return Status::Invalid("Failed to create folder", path);
        }

        return Status::OK();
    }

    // Write out a gvcf based on information parsed from yaml file
    Status write_gvcfs(const YAML::Node& n_gvcfs, bool is_input) {

        Status s;
        V(n_gvcfs.IsSequence(), "gvcf YAML node is not a sequence at top level.");

        // Iterate through sequence of input_gvcf
        for(YAML::const_iterator gvcf = n_gvcfs.begin(); gvcf != n_gvcfs.end(); ++gvcf) {
            V(gvcf->IsMap(), "input_gvcf is not a map");

            // Each input_gvcf is a map from filename to vcf text
            for(YAML::const_iterator it=gvcf->begin(); it != gvcf->end(); ++it) {
                string fn = it->first.as<string>();
                string vcf_body = it->second.as<string>();

                // Write each gvcf to file
                string gvcf_out_path = temp_dir_path + fn;
                ofstream fout;
                fout.open(gvcf_out_path);

                // header is shared amongst all input gvcfs
                if (is_input)
                    fout << header << "\t" << vcf_body;
                else
                    fout << vcf_body;

                fout.close();

                if (is_input){
                    // Track gvcf filename fpr input_gvcfs
                    input_gvcfs.insert(fn);
                } else {
                    // Specify the full filepath
                    truth_vcf_path = gvcf_out_path;
                }

            } // Close iterator over gvcf map
        } // Close iterator over input sequence
        return Status::OK();
    }

    // Parses unified_sites from yaml file to memory
    Status parse_unified_sites(const YAML::Node& n_unified_sites, const vector<pair<string,size_t> >& contigs) {
        Status s;

        for (YAML::const_iterator site_p = n_unified_sites.begin(); site_p != n_unified_sites.end(); ++site_p) {

            V(site_p->IsMap(), "unified_site is not a map at top level");

            unified_site site(range(-1, -1, -1));
            s = unified_site::of_yaml(*site_p, contigs, site);
            REQUIRE(s.ok());
            truth_sites.push_back(site);
        }
        return Status::OK();
    }

    // Print headers for cases where unexpected results were obtained
    void print_header() {
        if (!has_printed_header){
            has_printed_header = true;
            cout << "====================================" << endl;
            cout << "Differences for gvcf_test_case: " << name << endl;
            cout << "====================================" << endl;
        }
    }

public:

    // Constructor takes in the name of the curated case, used
    // for locating input yaml file; by convention, the yaml file is found
    // at INPUT_ROOT_DIR/<name>.yml
    GVCFTestCase(const string _name, bool _test_dsals=false,
                 bool _test_usites=true, bool _test_genotypes=true): test_dsals(_test_dsals), test_usites(_test_usites), test_genotypes(_test_genotypes), name(_name) {
        validated_formats = {"GT", "RNC"};
        validated_infos = {};
    }

    // Constructor with specification for format and info fields
    // for validation
    GVCFTestCase(const string _name, vector<string> _validated_formats,
                 vector<string> _validated_infos, bool _test_dsals=false,
                 bool _test_usites=true, bool _test_genotypes=true): test_dsals(_test_dsals), test_usites(_test_usites), test_genotypes(_test_genotypes), validated_formats(_validated_formats), validated_infos(_validated_infos), name(_name) {}

    // Reads in a yml file; writes gvcf records contained within yaml
    // file to disk; load expected unified_sites to memory
    Status load_yml() {
        Status s;

        string path = INPUT_ROOT_DIR + name + ".yml";
        YAML::Node yaml = YAML::LoadFile(path);
        V(yaml.IsMap(), "not a map at top level");

        const auto n_unifier_config = yaml["unifier_config"];
        if (n_unifier_config) {
            S(unifier_config::of_yaml(n_unifier_config, unifier_cfg));
        }

        const auto n_genotyper_config = yaml["genotyper_config"];
        if (n_genotyper_config) {
            S(genotyper_config::of_yaml(n_genotyper_config, genotyper_cfg));
        }

        const auto n_input = yaml["input"];
        const auto n_header = n_input["header"];
        V(n_header.IsScalar(), "header string invalid");
        header = n_header.Scalar();

        // Stage temp folders for output
        S(create_folder(TEMP_DIR));
        temp_dir_path = TEMP_DIR + name + "/";
        S(create_folder(temp_dir_path));

        // write gvcf entries into flat file
        // all test cases should have input gvcf
        const auto n_gvcfs = n_input["body"];
        S(write_gvcfs(n_gvcfs, true));

        s = VCFData::Open(input_gvcfs, data, temp_dir_path);
        REQUIRE(s.ok());
        s = data->contigs(contigs);

        s = Service::Start(service_config(), *data, *data, svc);
        REQUIRE(s.ok());

        // parse discovered_site if test_dsals is set to true
        if (test_dsals) {
            const auto n_dsals = yaml["truth_discovered_alleles"];
            REQUIRE(n_dsals);
            s = discovered_alleles_of_yaml(n_dsals, contigs, truth_dsals);
            REQUIRE(s.ok());
            string tmpfile = "/tmp/discover_alleles_capnp_check";
            REQUIRE(GLnexus::capnp_discover_alleles_verify(1, contigs, truth_dsals, tmpfile).ok());
        }

        // parse truth unified_sites if test_usites is set to true
        if (test_usites) {
            const auto n_unified_sites = yaml["truth_unified_sites"];
            REQUIRE(n_unified_sites);
            S(parse_unified_sites(n_unified_sites, contigs));
            REQUIRE(GLnexus::capnp_unified_sites_verify(truth_sites, "/tmp/unified_sites_capnp_check").ok());
        }

        // write truth gvcf if test_genotypes is set to true
        if (test_genotypes) {
            const auto n_truth_gvcf = yaml["truth_output_vcf"];
            REQUIRE(n_truth_gvcf);
            S(write_gvcfs(n_truth_gvcf, false))
        }
        return Status::OK();
    }

    Status execute_discover_alleles(discovered_alleles& als, const string sampleset, range pos) {
        unsigned N;
        Status s = svc->discover_alleles(sampleset, pos, N, als);

        // verify that the discovered-alleles structure can be serialized correctly with cap'n proto
        string dbg_filename("/tmp/allele_verify.cflat");
        REQUIRE(capnp_discover_alleles_verify(N, contigs, als, dbg_filename).ok());

        return s;
    }

    Status check_discovered_alleles(discovered_alleles& als, const string sampleset="<ALL>", range pos=range(0, 0, 1000000000)) {
        Status s = execute_discover_alleles(als, sampleset, pos);
        REQUIRE(s.ok());

        bool eq = true;
        auto lhs = als.begin(), rhs = truth_dsals.begin();
        for (; lhs != als.end() && rhs != truth_dsals.end(); ++lhs, ++rhs) {
            if (lhs->first != rhs->first || lhs->second != rhs->second) {
                if (eq) print_header();
                cout << lhs->first.str() << ":" << lhs->second.str() << endl;
                cout << rhs->first.str() << ":" << rhs->second.str() << endl;
                eq = false;
            }
        }
        if (eq && (lhs != als.end() || rhs != truth_dsals.end())) {
            print_header();
            cout << "Expected discovered alleles \n";
            for (auto& dsal : truth_dsals) {
                cout << dsal.first.str() << ":" << dsal.second.str() << endl;
            }

            cout << "Generated discovered alleles \n";
            for (auto &dsal : als) {
                cout << dsal.first.str() << ":" << dsal.second.str() << endl;
            }
            eq = false;
        }
        REQUIRE(eq);

        REQUIRE(als == truth_dsals);
        return Status::OK();
    }

    // Verify that unified_site generated by GLnexus service matches
    // the expected "truth"
    Status check_unify_sites(range pos=range(0, 0, 1000000000)) {
        discovered_alleles als;
        unsigned N;
        Status s = svc->discover_alleles("<ALL>", pos, N, als);
        REQUIRE(s.ok());

        unifier_stats stats;
        s = unified_sites(unifier_cfg, N, als, sites, stats);
        REQUIRE(s.ok());

        REQUIRE(is_sorted(sites.begin(), sites.end()));
        sort(truth_sites.begin(), truth_sites.end());

        // Print debug comparison before failing
        if (sites != truth_sites) {
            print_header();
            cout << "Non-matching truth sites \n";
            for (auto& site : truth_sites) {
                bool match = false;
                for (auto& rhs : sites) {
                    if (site == rhs) {
                        match = true;
                    }
                }
                if (!match) {
                    YAML::Emitter yaml;
                    REQUIRE(site.yaml(contigs, yaml).ok());
                    cout << yaml.c_str() << endl;
                }
            }

            cout << "Non-matching generated sites \n";
            for (auto& site: sites) {
                bool match = false;
                for (auto& rhs : truth_sites) {
                    if (site == rhs) {
                        match = true;
                    }
                }
                if (!match) {
                    YAML::Emitter yaml;
                    REQUIRE(site.yaml(contigs, yaml).ok());
                    cout << yaml.c_str() << endl;
                }
            }
        }
        REQUIRE(sites == truth_sites);
        return Status::OK();
    }

    // Verify that the genotypes generated by GLnexus service matches
    // the "truth" vcf output. Assumes that check_unified_sites has
    // already been called, and sites have been populated.
    Status check_genotypes() {
        Status s;

        // There should be at least 1 unified_site, this will fail
        // if check_unified_sites() has not been called beforehand
        REQUIRE(!sites.empty());

        genotyper_cfg.output_format = GLnexusOutputFormat::VCF;
        string out_vcf_path = temp_dir_path + "output.vcf";
        s = svc->genotype_sites(genotyper_cfg, string("<ALL>"), sites, out_vcf_path);
        if (! s.ok()) {
            cout << s.str() << endl;
        }
        REQUIRE(s.ok());

        string diff_out_path = temp_dir_path +  "output.diff";

        string diff_cmd = "python testOutputVcf.py --input " + out_vcf_path + " --truth " + truth_vcf_path;
        diff_cmd += " --quiet";
        if (!validated_formats.empty()) {
            diff_cmd += " --formats ";
            for (auto& validated_format : validated_formats) {
                diff_cmd += validated_format + " ";
            }
        }
        if (!validated_infos.empty()) {
            diff_cmd += " --infos ";
            for (auto& validated_info : validated_infos) {
                diff_cmd +=  validated_info + " ";
            }
        }
        diff_cmd += " > " + diff_out_path;
        int retval = system(diff_cmd.c_str());
        REQUIRE (WIFEXITED(retval));

        int exit_status = WEXITSTATUS(retval);
        if (exit_status != 0) {
            print_header();
            cout << "^^^^^ Output vcf and truth differs. Diff file: " << endl;
            retval = system(("cat " + diff_out_path).c_str());
            REQUIRE(retval == 0);
        }

        REQUIRE(exit_status == 0);

        return Status::OK();
    }

    Status cleanup() {
        int retval = system(("rm -rf " + temp_dir_path).c_str());
        if (retval == 0)
            return Status::OK();
        else
            return Status::Invalid("Failed to remove directory", temp_dir_path);
    }

    Status perform_gvcf_test(range pos=range(0, 0, 1000000000)) {
        Status s;
        s = this->load_yml();
        if (!s.ok()) cout << s.str() << endl;
        REQUIRE(s.ok());

        if (test_dsals) {
            discovered_alleles als;
            s = this->check_discovered_alleles(als, "<ALL>", pos);
            REQUIRE(s.ok());
        }
        if (test_usites) {
            s = this->check_unify_sites(pos);
            REQUIRE(s.ok());
        }
        if (test_genotypes) {
            s = this->check_genotypes();
            REQUIRE(s.ok());
        }
        //s = this->cleanup();
        REQUIRE(s.ok());
        return Status::OK();
    }
};

TEST_CASE("curated case: trim_input") {
    GVCFTestCase trim_input_case("trim_input");
    trim_input_case.perform_gvcf_test();
}

TEST_CASE("curated case: join_gvcfs") {
    GVCFTestCase join_gvcfs_case("join_gvcfs");
    join_gvcfs_case.perform_gvcf_test();
}
TEST_CASE("synthetic case: join_3_gvcfs") {
    GVCFTestCase join_3_gvcfs_case("join_3_gvcfs");
    join_3_gvcfs_case.perform_gvcf_test();
}

TEST_CASE("curated case: join_vcf_gvcf") {
    GVCFTestCase join_vcf_gvcf_case("join_vcf_gvcf");
    join_vcf_gvcf_case.perform_gvcf_test();
}

TEST_CASE("synthetic case: join_gvcf_vcf") {
    GVCFTestCase join_gvcf_vcf_case("join_gvcf_vcf");
    join_gvcf_vcf_case.perform_gvcf_test();
}

TEST_CASE("synthetic case: join_gvcf_vcf_gvcf") {
    GVCFTestCase join_gvcf_vcf_gvcf_case("join_gvcf_vcf_gvcf");
    join_gvcf_vcf_gvcf_case.perform_gvcf_test();
}

TEST_CASE("curated case: inconsistent_trim") {
    GVCFTestCase inconsistent_trim_case("inconsistent_trim");
    inconsistent_trim_case.perform_gvcf_test();
}

TEST_CASE("services test: discover_alleles") {
    GVCFTestCase discover_allele_case("discover_alleles_2_trios", true, false, false);
    discover_allele_case.load_yml();
    discovered_alleles als;

    Status s;

    SECTION("2 trios") {
        s = discover_allele_case.check_discovered_alleles(als, "<ALL>", range(0, 0, 1099));
        REQUIRE(s.ok());
    }

    SECTION("nonexistent sampleset") {
        s = discover_allele_case.execute_discover_alleles(als, "bogus", range(0, 0, 1000000));
        REQUIRE(s == StatusCode::NOT_FOUND);
    }

    SECTION("trio1") {
        s = discover_allele_case.execute_discover_alleles(als, "discover_alleles_trio1", range(0, 0, 1099));
        REQUIRE(s.ok());
        REQUIRE(als.size() == 7);
        auto p = als.find(allele(range(0, 1000, 1001), "G"));
        REQUIRE(p != als.end());
        REQUIRE(p->second.zGQ.copy_number() == 4);
        p = als.find(allele(range(0, 1010, 1012), "CC"));
        REQUIRE(p != als.end());
        REQUIRE(p->second.zGQ.copy_number() == 3);
    }

    SECTION("trio1 partial") {
        s = discover_allele_case.execute_discover_alleles(als, "discover_alleles_trio1", range(0, 1009, 1011));
        REQUIRE(s.ok());

        REQUIRE(als.size() == 2);

        s = discover_allele_case.execute_discover_alleles(als, "discover_alleles_trio1", range(0, 1009, 1012));
        REQUIRE(s.ok());

        REQUIRE(als.size() == 2);
        REQUIRE(als.find(allele(range(0, 1010, 1012), "CC"))->second.zGQ.copy_number() == 3);
    }

    SECTION("spanning allele") {
        s = discover_allele_case.execute_discover_alleles(als, "<ALL>", range(1, 1010, 1012));
        REQUIRE(s.ok());
        REQUIRE(als.size() == 4);
        REQUIRE(als.find(allele(range(1, 1001, 1016), "AAAAAAAAAAAAAAA")) != als.end());

        s = discover_allele_case.execute_discover_alleles(als, "<ALL>", range(1, 1001, 1016));
        REQUIRE(s.ok());
        REQUIRE(als.size() == 6);
        REQUIRE(als.find(allele(range(1, 1001, 1016), "AAAAAAAAAAAAAAA"))->second.zGQ.copy_number() == 3);
    }

    SECTION("detect inconsistent reference alleles") {
        s = discover_allele_case.execute_discover_alleles(als, "<ALL>", range(2, 1000, 1004));
        REQUIRE_FALSE(s.ok());
        REQUIRE(s == StatusCode::INVALID);
        REQUIRE(s.str().find("allele appears as both REF and ALT") != string::npos);
    }

    SECTION("calculate copy numbers from desired samples only") {
        // looking at all three samples in the dataset, we should receive all
        // three alleles at this position
        s = discover_allele_case.execute_discover_alleles(als, "discover_alleles_trio1", range(0, 1001, 1002));
        REQUIRE(s.ok());
        REQUIRE(als.size() == 3);
        REQUIRE(als.find(allele(range(0, 1001, 1002), "C"))->second.zGQ.copy_number() == 2);
        REQUIRE(als.find(allele(range(0, 1001, 1002), "G"))->second.zGQ.copy_number() == 2);
        REQUIRE(als.find(allele(range(0, 1001, 1002), "T"))->second.zGQ.copy_number() == 2);

        // looking only at one homozygous ref individual
        s = discover_allele_case.execute_discover_alleles(als, "trio1.fa", range(0, 1001, 1002));
        REQUIRE(s.ok());
        REQUIRE(als.size() == 0);


        s = discover_allele_case.execute_discover_alleles(als, "discover_alleles_trio1", range(0, 1000, 1001));
        REQUIRE(s.ok());
        REQUIRE(als.size() == 2);
        REQUIRE(als.find(allele(range(0, 1000, 1001), "A"))->second.zGQ.copy_number() == 2);
        REQUIRE(als.find(allele(range(0, 1000, 1001), "G"))->second.zGQ.copy_number() == 4);

        s = discover_allele_case.execute_discover_alleles(als, "trio1.fa", range(0, 1000, 1001));
        REQUIRE(s.ok());
        REQUIRE(als.size() == 2);
        REQUIRE(als.find(allele(range(0, 1000, 1001), "A"))->second.zGQ.copy_number() == 1);
        REQUIRE(als.find(allele(range(0, 1000, 1001), "G"))->second.zGQ.copy_number() == 1);
    }

    discover_allele_case.cleanup();
}

TEST_CASE("services test: discover_alleles_gVCF") {
    GVCFTestCase discover_allele_case("discover_alleles_gvcf", true, false, false);
    discover_allele_case.load_yml();
    discovered_alleles als;

    Status s = discover_allele_case.check_discovered_alleles(als, "NA12878", range(0, 10009460, 10009470));

    REQUIRE(s.ok());

    discover_allele_case.cleanup();
}

TEST_CASE("services test: discover_alleles_gVCF_bogus") {
    GVCFTestCase discover_allele_case("discover_alleles_gvcf_bogus", true, false, false);
    discover_allele_case.load_yml();
    discovered_alleles als;

    Status s = discover_allele_case.check_discovered_alleles(als, "<ALL>", range(0, 10009460, 10009465));
    REQUIRE(s.ok());

    s = discover_allele_case.execute_discover_alleles(als, "<ALL>", range(0, 10009465, 10009466));
    REQUIRE(s == StatusCode::INVALID);

    discover_allele_case.cleanup();
}

TEST_CASE("services test: VCF file services test") {
    GVCFTestCase vcf_services_case("vcf_services");
    vcf_services_case.perform_gvcf_test();
}

TEST_CASE("services test: gVCF file services test") {
    GVCFTestCase gvcf_services_case("gvcf_services");
    gvcf_services_case.perform_gvcf_test();
}

TEST_CASE("services test: gVCF file services test depth12") {
    GVCFTestCase gvcf_services_case("gvcf_services_depth12");
    gvcf_services_case.perform_gvcf_test();
}

TEST_CASE("services test: gVCF file services test depth9") {
    GVCFTestCase gvcf_services_case("gvcf_services_depth9");
    gvcf_services_case.perform_gvcf_test();
}

TEST_CASE("join record logic test: basic") {
    GVCFTestCase join_record_case("join_records_basic");
    join_record_case.perform_gvcf_test();
}

TEST_CASE("join record logic test: incomplete span") {
    GVCFTestCase join_record_case("join_records_incomplete_span");
    join_record_case.perform_gvcf_test();
}

TEST_CASE("join record logic test: unjoinable (multiple variant records)") {
    GVCFTestCase join_record_case("join_records_unjoinable");
    join_record_case.perform_gvcf_test();
}

TEST_CASE("join record logic test: 0 ref depth requirement") {
    GVCFTestCase join_record_case("join_records_basic");
    join_record_case.perform_gvcf_test();
}

TEST_CASE("join record logic test: insufficient ref depth") {

    GVCFTestCase join_record_case("join_records_insufficient_ref_depth");
    join_record_case.perform_gvcf_test();
}

TEST_CASE("liftover field: combi") {
    vector<string> v_formats = {"GT", "RNC", "AD", "DP", "SB"};
    vector<string> v_infos = {};
    GVCFTestCase format_fields_combi("format_fields_combi", v_formats, v_infos);
    format_fields_combi.perform_gvcf_test();
}

TEST_CASE("liftover field: AD_fall_back") {
    vector<string> v_formats = {"GT", "RNC", "AD"};
    vector<string> v_infos = {};
    GVCFTestCase format_fields_AD("format_fields_AD_fallback", v_formats, v_infos);
    format_fields_AD.perform_gvcf_test();
}

TEST_CASE("liftover field: default_type=missing") {
    vector<string> v_formats = {"GT", "RNC", "AD", "SB", "SBF"};
    vector<string> v_infos = {};
    GVCFTestCase format_field_missing("format_fields_default_missing", v_formats, v_infos);
    format_field_missing.perform_gvcf_test();
}

TEST_CASE("liftover field: default_type=none") {
    vector<string> v_formats = {"GT", "RNC", "AD", "SB"};
    vector<string> v_infos = {};
    GVCFTestCase format_field_none("format_fields_default_none", v_formats, v_infos);
    format_field_none.perform_gvcf_test();
}

TEST_CASE("liftover field: ignore_non_variants") {

    vector<string> v_formats = {"GT", "RNC", "DP", "SB", "AD", "GQ"};
    vector<string> v_infos = {};
    GVCFTestCase format_field_ignore_non_variants("format_fields_ignore_non_variants", v_formats, v_infos);
    format_field_ignore_non_variants.perform_gvcf_test();
}

TEST_CASE("liftover field: integration") {

    vector<string> v_formats = {"GT", "RNC", "DP", "SB", "AD", "GQ"};
    vector<string> v_infos = {};
    GVCFTestCase join_record_case("format_fields_integrated", v_formats, v_infos);
    join_record_case.perform_gvcf_test();
}

TEST_CASE("join record logic test: overlapping records") {
    GVCFTestCase join_record_case("join_records_overlapping");
    join_record_case.perform_gvcf_test();
}

TEST_CASE("lost deletion") {
    GVCFTestCase("lost_deletion").perform_gvcf_test();
}

TEST_CASE("lost deletion_monoallelic_site") {
    vector<string> v_formats = {"GT", "RNC", "DP", "AD"};
    GVCFTestCase("lost_deletion_monoallelic_site", v_formats, {}).perform_gvcf_test();
}

TEST_CASE("join records with unifier preference for small alleles") {
    GVCFTestCase("join_records_prefer_small").perform_gvcf_test();
}

TEST_CASE("DP0_noAD") {
    vector<string> v_formats = {"GT", "RNC", "DP", "SB", "AD", "GQ"};
    vector<string> v_infos = {};
    GVCFTestCase DP0_case("DP0_noAD", v_formats, v_infos);
    DP0_case.perform_gvcf_test();
}

TEST_CASE("min_allele_copy_number") {
    GVCFTestCase("min_allele_copy_number").perform_gvcf_test();
}

TEST_CASE("min_quality") {
    GVCFTestCase("min_quality",true,true,true).perform_gvcf_test();
}

TEST_CASE("rs141305015") {
    GVCFTestCase("rs141305015",true).perform_gvcf_test();
}

TEST_CASE("rs11429009 allele unifier test") {
    GVCFTestCase("rs11429009").perform_gvcf_test();
}

TEST_CASE("revise_overlapping") {
    vector<string> v_formats = {"GT", "RNC", "GQ"};
    vector<string> v_infos = {};
    GVCFTestCase revise_overlapping_case("revise_overlapping", v_formats, v_infos);
    revise_overlapping_case.perform_gvcf_test();
}

TEST_CASE("censor_rnc_format_fields") {
    vector<string> v_formats = {"DP", "GT", "RNC", "GQ", "AD"};
    vector<string> v_infos = {};
    GVCFTestCase gtc("censor_rnc_format_fields", v_formats, v_infos);
    gtc.perform_gvcf_test();
}

TEST_CASE("xAtlas") {
    vector<string> v_formats = {"DP", "GT", "GQ", "PL", "RR", "VR", "FT"};
    vector<string> v_infos = {};
    GVCFTestCase("xAtlas", v_formats, v_infos, true).perform_gvcf_test();
}

TEST_CASE("weCall") {
    vector<string> v_formats = {"DP", "GT", "GQ", "PL", "AD", "FT"};
    vector<string> v_infos = {};
    GVCFTestCase("weCall", v_formats, v_infos, false).perform_gvcf_test();
}

TEST_CASE("edge_spanning_deletion") {
    vector<string> v_formats = {"GT", "RNC", "AD", "DP", "SB"};
    vector<string> v_infos = {};
    GVCFTestCase edge_spanning_deletion("edge_spanning_deletion", v_formats, v_infos);
    edge_spanning_deletion.perform_gvcf_test(range(0, 999, 1001));
    // Repeat test with slightly different range, such that the deletion dangles from
    // upstream of the target region into it.
    GVCFTestCase edge_spanning_deletion_1000("edge_spanning_deletion", v_formats, v_infos);
    edge_spanning_deletion_1000.perform_gvcf_test(range(0, 1000, 1001));
}

TEST_CASE("edge_spanning_deletion2") {
    vector<string> v_formats = {"GT", "RNC", "AD", "DP", "SB"};
    vector<string> v_infos = {};
    GVCFTestCase edge_spanning_deletion("edge_spanning_deletion2", v_formats, v_infos);
    edge_spanning_deletion.perform_gvcf_test(range(0, 999, 1000));
}

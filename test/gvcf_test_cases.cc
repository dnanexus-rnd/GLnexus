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
            if (!s.ok()) cout<<s.str()<<endl;
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
                 bool _test_usites=true, bool _test_genotypes=true): name(_name), test_dsals(_test_dsals), test_usites(_test_usites), test_genotypes(_test_genotypes) {
        validated_formats = {"GT", "RNC"};
        validated_infos = {};
    }

    // Constructor with specification for format and info fields
    // for validation
    GVCFTestCase(const string _name, vector<string> _validated_formats,
                 vector<string> _validated_infos, bool _test_dsals=false,
                 bool _test_usites=true, bool _test_genotypes=true): name(_name), test_dsals(_test_dsals), test_usites(_test_usites), test_genotypes(_test_genotypes),validated_formats(_validated_formats), validated_infos(_validated_infos) {}

    // Reads in a yml file; writes gvcf records contained within yaml
    // file to disk; load expected unified_sites to memory
    Status load_yml() {
        Status s;

        string path = INPUT_ROOT_DIR + name + ".yml";
        YAML::Node yaml = YAML::LoadFile(path);
        V(yaml.IsMap(), "not a map at top level");

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
        }

        // parse truth unified_sites if test_usites is set to true
        if (test_usites) {
            const auto n_unified_sites = yaml["truth_unified_sites"];
            REQUIRE(n_unified_sites);
            S(parse_unified_sites(n_unified_sites, contigs));
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
        return svc->discover_alleles(sampleset, pos, als);
    }

    Status check_discovered_alleles(discovered_alleles& als, const string sampleset="<ALL>", range pos=range(0, 0, 1000000000)) {
        Status s = execute_discover_alleles(als, sampleset, pos);
        REQUIRE(s.ok());

        // Print comparison pf dsals before failing
        if (als != truth_dsals) {
            print_header();
            cout << "Expected discovered alleles \n";
            for (auto& dsal : truth_dsals) {
                cout << dsal.first.str() << ":" << dsal.second.str() << endl;
            }

            cout << "Generated discovered alleles \n";
            for (auto &dsal : als) {
                cout << dsal.first.str() << ":" << dsal.second.str() << endl;
            }
        }

        REQUIRE(als == truth_dsals);
        return Status::OK();
    }

    // Verify that unified_site generated by GLnexus service matches
    // the expected "truth"
    Status check_unify_sites() {
        discovered_alleles als;
        Status s = svc->discover_alleles("<ALL>", range(0, 0, 1000000000), als);
        REQUIRE(s.ok());

        s =unified_sites(unifier_config(), als, sites);
        REQUIRE(s.ok());

        REQUIRE(is_sorted(sites.begin(), sites.end()));
        sort(truth_sites.begin(), truth_sites.end());

        // Print debug comparison before failing
        if (sites != truth_sites) {
            print_header();
            cout << "Expected truth sites \n";
            for (auto& site : truth_sites) {
                YAML::Emitter yaml;
                REQUIRE(site.yaml(contigs, yaml).ok());
                cout << yaml.c_str() << endl;
            }

            cout << "Generated unified sites \n";
            for (auto& site: sites) {
                YAML::Emitter yaml;
                REQUIRE(site.yaml(contigs, yaml).ok());
                cout << yaml.c_str() << endl;
            }
        }
        REQUIRE(sites == truth_sites);
        return Status::OK();
    }

    // Verify that the genotypes generated by GLnexus service matches
    // the "truth" vcf output. Assumes that check_unified_sites has
    // already been called, and sites have been populated.
    Status check_genotypes(const genotyper_config& config) {
        Status s;

        // There should be at least 1 unified_site, this will fail
        // if check_unified_sites() has not been called beforehand
        REQUIRE(!sites.empty());

        consolidated_loss losses;
        string out_vcf_path = temp_dir_path + "output.vcf";
        s = svc->genotype_sites(config, string("<ALL>"), sites, out_vcf_path, losses);

        cout << s.str() << endl;
        if (!s.ok()) {
            cout << s.str() << endl;
        }
        REQUIRE(s.ok());

        string diff_out_path = temp_dir_path +  "output.diff";

        string diff_cmd = "python testOutputVcf.py --input " + out_vcf_path + " --truth " + truth_vcf_path;
        // string diff_cmd += " --quiet";
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

    Status perform_gvcf_test(const genotyper_config& config=genotyper_config(GLnexusOutputFormat::VCF)) {
        Status s;
        s = this->load_yml();
        REQUIRE(s.ok());

        if (test_dsals) {
            discovered_alleles als;
            s = this->check_discovered_alleles(als);
            REQUIRE(s.ok());
        }
        if (test_usites) {
            s = this->check_unify_sites();
            REQUIRE(s.ok());
        }
        if (test_genotypes) {
            s = this->check_genotypes(config);
            REQUIRE(s.ok());
        }
        s = this->cleanup();
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
        REQUIRE(p->second.observation_count == 4);
        p = als.find(allele(range(0, 1010, 1012), "CC"));
        REQUIRE(p != als.end());
        REQUIRE(p->second.observation_count == 3);
    }

    SECTION("trio1 partial") {
        s = discover_allele_case.execute_discover_alleles(als, "discover_alleles_trio1", range(0, 1009, 1011));
        REQUIRE(s.ok());

        REQUIRE(als.size() == 0);

        s = discover_allele_case.execute_discover_alleles(als, "discover_alleles_trio1", range(0, 1009, 1012));
        REQUIRE(s.ok());

        REQUIRE(als.size() == 2);
        REQUIRE(als.find(allele(range(0, 1010, 1012), "CC"))->second.observation_count == 3);
    }

    SECTION("spanning allele") {
        s = discover_allele_case.execute_discover_alleles(als, "<ALL>", range(1, 1010, 1012));
        REQUIRE(s.ok());
        REQUIRE(als.size() == 2);
        REQUIRE(als.find(allele(range(1, 1001, 1016), "AAAAAAAAAAAAAAA")) == als.end());

        s = discover_allele_case.execute_discover_alleles(als, "<ALL>", range(1, 1001, 1016));
        REQUIRE(s.ok());
        REQUIRE(als.size() == 6);
        REQUIRE(als.find(allele(range(1, 1001, 1016), "AAAAAAAAAAAAAAA"))->second.observation_count == 3);
    }

    SECTION("detect inconsistent reference alleles") {
        s = discover_allele_case.execute_discover_alleles(als, "<ALL>", range(2, 1000, 1004));
        REQUIRE_FALSE(s.ok());
        REQUIRE(s == StatusCode::INVALID);
        REQUIRE(s.str().find("allele appears as both REF and ALT") != string::npos);
    }

    SECTION("only return alleles observed in desired samples") {
        // looking at all three samples in the dataset, we should receive all
        // three alleles at this position
        s = discover_allele_case.execute_discover_alleles(als, "discover_alleles_trio1", range(0, 1001, 1002));
        REQUIRE(s.ok());
        REQUIRE(als.size() == 3);

        // looking only at one homozygous ref individual, we should get
        // nothing
        s = discover_allele_case.execute_discover_alleles(als, "trio1.fa", range(0, 1001, 1002));
        REQUIRE(s.ok());
        REQUIRE(als.size() == 0);

        // throw in a previous position as positive control
        s = discover_allele_case.execute_discover_alleles(als, "trio1.fa", range(0, 1000, 1002));
        REQUIRE(s.ok());
        REQUIRE(als.size() == 2);
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

    genotyper_config cfg(GLnexusOutputFormat::VCF);
    cfg.required_dp = 13;

    gvcf_services_case.perform_gvcf_test(cfg);
}

TEST_CASE("services test: gVCF file services test depth9") {
    GVCFTestCase gvcf_services_case("gvcf_services_depth9");

    genotyper_config cfg(GLnexusOutputFormat::VCF);
    cfg.required_dp = 10;

    gvcf_services_case.perform_gvcf_test(cfg);
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

    genotyper_config cfg(GLnexusOutputFormat::VCF);
    cfg.required_dp = 1;

    join_record_case.perform_gvcf_test(cfg);
}

TEST_CASE("join record logic test: insufficient ref depth") {

    GVCFTestCase join_record_case("join_records_insufficient_ref_depth");

    genotyper_config cfg(GLnexusOutputFormat::VCF);
    cfg.required_dp = 13;

    join_record_case.perform_gvcf_test(cfg);
}

TEST_CASE("join record logic test: test liftover fields") {

    vector<string> v_formats = {"GT", "RNC", "DP", "SB", "AD", "GQ"};
    vector<string> v_infos = {};
    GVCFTestCase join_record_case("join_records_formats_basic", v_formats, v_infos);

    genotyper_config cfg(GLnexusOutputFormat::VCF);

    vector<string> orig_names;

    retained_format_field dp_field = retained_format_field({"DP", "MIN_DP"}, "DP", RetainedFieldType::INT, FieldCombinationMethod::MIN, RetainedFieldNumber::BASIC, 1);
    dp_field.description = "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">";

    retained_format_field gq_field = retained_format_field({"GQ"}, "GQ", RetainedFieldType::INT, FieldCombinationMethod::MAX, RetainedFieldNumber::BASIC, 1);
    gq_field.description = "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">";

    retained_format_field sb_field = retained_format_field({"SB"}, "SB", RetainedFieldType::INT, FieldCombinationMethod::MAX, RetainedFieldNumber::BASIC, 4);
    sb_field.description = "##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.\">";

    retained_format_field ad_field = retained_format_field({"AD"}, "AD", RetainedFieldType::INT, FieldCombinationMethod::MIN,
        RetainedFieldNumber::ALLELES);
    ad_field.description = "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">";

    cfg.liftover_fields.push_back(ad_field);
    cfg.liftover_fields.push_back(dp_field);
    cfg.liftover_fields.push_back(gq_field);
    cfg.liftover_fields.push_back(sb_field);
    join_record_case.perform_gvcf_test(cfg);
}

TEST_CASE("join record logic test: overlapping records") {
    GVCFTestCase join_record_case("join_records_overlapping");
    join_record_case.perform_gvcf_test();
}

TEST_CASE("lost deletion") {
    GVCFTestCase("lost_deletion").perform_gvcf_test();
}

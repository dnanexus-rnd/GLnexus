#include <iostream>
#include <fstream>
#include "service.h"
#include "types.h"
#include "utils.cc"
#include "catch.hpp"
#include "yaml-cpp/yaml.h"
using namespace std;
using namespace GLnexus;

#define test_contigs {{"A", 0},}

/// CuratedCase is the overarching class for unit-testing
/// of a curated "case" for the glnexus workflow.
/// It handles parsing of a yaml file with information on
/// gvcf input. It supports comparing output of the GLnexus workflow
/// with the expected unified_sites and joint genotype output.

class CuratedCase {
    #define V(pred,msg) if (!(pred)) return Status::Invalid("CuratedCase::load_yml: " msg)

    // Name of curated case
    string name;

    // Filenames of input gvcfs
    set<string> input_gvcfs;

    // Filename of truth vcf
    string truth_vcf_fn;

    // Parsed "truth" unfiied sites
    vector<unified_site> truth_sites;

    // Unified sites based on unify_sites algorithm
    vector<unified_site> sites;

    // Storage for loaded input gvcfs
    unique_ptr<VCFData> data;

    // Service for algorithmic procedures
    unique_ptr<Service> svc;

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
                string gvcf_out_path = "test/data/special_cases/" + name + "/" + fn;
                ofstream fout;
                fout.open(gvcf_out_path);
                fout << vcf_body;
                fout.close();

                // Track gvcf filename fpr input_gvcfs
                if (is_input){
                    // the structure test/data is embedded in VCFData::Open
                    input_gvcfs.insert("special_cases/" + name + "/" + fn);
                } else {
                    // Specify the full filepath
                    truth_vcf_fn = "test/data/special_cases/" + name + "/" + fn;
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
            s = unified_site_of_yaml(*site_p, contigs, site);
            if (!s.ok()) cout<<s.str()<<endl;
            REQUIRE(s.ok());
            truth_sites.push_back(site);
        }
        return Status::OK();
    }

public:

    // Constructor takes in the name of the curated case, used
    // for locating input yaml file; by convention, the yaml file is found
    // at test/data/special_cases/<name>/<name>.yml
    CuratedCase(const string _name): name(_name) {
        cout << "====================================" << endl;
        cout << "Start of test for curated case: " << name << endl;
        cout << "====================================" << endl;
    }

    // Reads in a yml file; writes gvcf records contained within yaml
    // file to disk; load expected unified_sites to memory
    Status load_yml() {
        Status s;

        string path = "test/data/special_cases/" + name + "/" + name + ".yml";
        YAML::Node yaml = YAML::LoadFile(path);
        V(yaml.IsMap(), "not a map at top level");

        // write gvcf entries into flat file
        const auto n_gvcfs = yaml["input"];
        S(write_gvcfs(n_gvcfs, true));

        // parse truth unified_sites
        const auto n_unified_sites = yaml["truth_unified_sites"];
        S(parse_unified_sites(n_unified_sites, test_contigs));

        // write truth gvcf
        const auto n_truth_gvcf = yaml["truth_output_vcf"];
        S(write_gvcfs(n_truth_gvcf, false))
        return Status::OK();
    }

    // Verify that unified_site generated by GLnexus service matches
    // the expected "truth"
    Status check_unify_sites() {
        Status s;
        S(VCFData::Open(input_gvcfs, data));
        s = Service::Start(*data, *data, svc);
        REQUIRE(s.ok());

        discovered_alleles als;
        s = svc->discover_alleles("<ALL>", range(0, 0, 1000000), als);
        REQUIRE(s.ok());

        s =unified_sites(als, sites);
        REQUIRE(s.ok());

        sort(sites.begin(), sites.end());
        sort(truth_sites.begin(), truth_sites.end());

        if (sites != truth_sites) {
            cout << "^^^^^ Unified sites differ from truth:";
            cout << "Expected truth sites \n";
            for (auto& site : truth_sites) {
                cout << site.str() << endl;
            }

            cout << "Generated unified sites \n";
            for (auto& site: sites) {
                cout << site.str() << endl;
            }
        } else {
            cout << "@@@@@ Unified sites agree!" << endl;
        }
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

        consolidated_loss losses;
        string out_fn = "test/data/special_cases/" + name + "/output.vcf";
        s = svc->genotype_sites(genotyper_config(), string("<ALL>"), sites, out_fn, losses, true);

        REQUIRE(s.ok());

        string diff_fn = "test/data/special_cases/" + name + "/output.diff";

        string diff_cmd = "diff " + out_fn + " " + truth_vcf_fn + " > " + diff_fn;

        int retval = system(diff_cmd.c_str());
        REQUIRE (WIFEXITED(retval));

        if (WEXITSTATUS(retval) == 0){
            // no difference found in diff command
            cout << "@@@@@ Output vcf and truth agrees!" << endl;
        } else {
            cout << "^^^^^ Output vcf and truth differs. Diff file: " << endl;
            retval = system(("cat " + diff_fn).c_str());
        }

        cout << "++++++++++++++++++++++++++++++++++++" << endl;
        cout << "End of test for Curated case: " << name << endl;
        cout << "++++++++++++++++++++++++++++++++++++" << endl;
        return Status::OK();
    }

    Status cleanup() {
        // TO-DO: delete intermediate gvcf/vcf files written to disk
        return Status::OK();
    }
};

TEST_CASE("trim_input") {
    CuratedCase trim_input_case("trim_input");
    trim_input_case.load_yml();
    trim_input_case.check_unify_sites();
    trim_input_case.check_genotypes();
}
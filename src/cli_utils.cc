#include <iostream>
#include <exception>
#include <fstream>
#include <sstream>
#include <regex>
#include "cli_utils.h"

// This file has utilities employed by the glnexus applet.
using namespace std;

namespace GLnexus {
namespace cli {
namespace utils {

// Parse a range like chr1:1000-2000. The item can also just be the name of a
// contig, in which case it gets mapped to the contig's full length.
bool parse_range(const vector<pair<string,size_t> >& contigs,
                 const string& range_txt, range& ans) {
    static regex re_range("([^:]+):([0-9,]+)-([0-9,]+)");

    string contig = range_txt;
    size_t beg=1, end=0;

    smatch sm;
    if (regex_match(range_txt, sm, re_range)) {
        assert(sm.size() == 4);
        contig = sm[1];
        string sbeg = sm[2], send = sm[3];
        sbeg.erase(std::remove(sbeg.begin(), sbeg.end(), ','), sbeg.end());
        send.erase(std::remove(send.begin(), send.end(), ','), send.end());
        beg = strtoul(sbeg.c_str(), nullptr, 10);
        end = strtoul(send.c_str(), nullptr, 10);
    }

    int rid=0;
    while (rid < contigs.size() && contig != contigs[rid].first) rid++;
    if (rid >= contigs.size()) {
        return false;
    }
    end = end>0 ? end : contigs[rid].second;
    if (beg < 1 || beg >= end) {
        return false;
    }
    ans = range(rid, beg-1, end);
    return true;
}

// parse a comma-separated list of ranges
bool parse_ranges(const vector<pair<string,size_t> >& contigs,
                  const string& ranges, vector<range>& ans) {
    ans.clear();

    string item;
    stringstream ss(ranges);
    while (std::getline(ss, item, ',')) {
        range range(-1,-1,-1);
        if (!parse_range(contigs, item, range)) {
            return false;
        }
        ans.push_back(range);
    }

    return true;
}

Status yaml_of_contigs_alleles_ranges(const vector<pair<string,size_t> > &contigs,
                                      const vector<range> &ranges,
                                      const vector<discovered_alleles> &valleles,
                                      YAML::Emitter &yaml) {
    Status s;

    if (ranges.size() != valleles.size())
        return Status::Invalid("number of ranges must equal number of sites");
    size_t nm_elem = ranges.size();

    yaml << YAML::BeginMap;

    {
        // write contigs
        yaml << YAML::Key << "contigs";
        yaml << YAML::Value;
        {
            yaml << YAML::BeginSeq;
            for (const std::pair<string,size_t> &pr : contigs) {
                yaml << YAML::BeginMap;
                yaml << YAML::Key << "name";
                yaml << YAML::Value << pr.first;
                yaml << YAML::Key << "size";
                yaml << YAML::Value << pr.second;
                yaml << YAML::EndMap;
            }
            yaml << YAML::EndSeq;
        }
    }

    {
        // Write ranges and alleles
        yaml << YAML::Key << "ranges_alleles";
        yaml << YAML::Value;
        {
            // Note: we need to be careful here, to write out only non empty sites.
            yaml << YAML::BeginSeq;
            for (int i=0; i < nm_elem; ++i) {
                if (valleles[i].size() == 0) continue;

                yaml << YAML::BeginMap;
                yaml << YAML::Key << "containing_range";
                yaml << YAML::Value;
                s = range_yaml(contigs, ranges[i], yaml);
                if (s.bad()) return s;

                yaml << YAML::Key << "discovered_alleles";
                yaml << YAML::Value;
                s = yaml_of_discovered_alleles(valleles[i], contigs, yaml);
                if (s.bad()) return s;
                yaml << YAML::EndMap;
            }
            yaml << YAML::EndSeq;
        }
    }
    yaml << YAML::EndMap;

    return Status::OK();
}

Status contigs_alleles_ranges_of_yaml(const YAML::Node& yaml,
                                      std::vector<std::pair<std::string,size_t> > &contigs,
                                      vector<range> &ranges,
                                      vector<discovered_alleles> &valleles) {
    Status s;
    if (!yaml.IsMap()) {
        return Status::Invalid("not a map at top level");
    }

    contigs.clear();
    ranges.clear();
    valleles.clear();

    // read contigs
    auto n_contigs = yaml["contigs"];
    if (!n_contigs.IsSequence()) {
        return Status::Invalid("contigs should be a yaml sequence");
    }
    for (auto p = n_contigs.begin(); p != n_contigs.end(); ++p) {
        const std::string name = (*p)["name"].as<std::string>();
        size_t size = (*p)["size"].as<size_t>();
        contigs.push_back(make_pair(name, size));
    }

    // read ranges and alleles
    auto n_ranges_alleles = yaml["ranges_alleles"];
    if (!n_ranges_alleles.IsSequence()) {
        return Status::Invalid("ranges and alleles should be a yaml sequence");
    }
    for (YAML::const_iterator p = n_ranges_alleles.begin(); p != n_ranges_alleles.end(); ++p) {
        const auto n_range = (*p)["containing_range"];
        range rng(-1,-1,-1);
        s = range_of_yaml(n_range, contigs, rng);
        if (s.bad()) return s;
        ranges.push_back(rng);

        const auto n_dsals = (*p)["discovered_alleles"];
        discovered_alleles dsals;
        s = discovered_alleles_of_yaml(n_dsals, contigs, dsals);
        if (s.bad()) return s;
        valleles.push_back(dsals);
    }

    return Status::OK();
}

// Serialize the unified sites to yaml format.
//
Status yaml_of_unified_sites(const vector<unified_site> &sites,
                             const vector<pair<string,size_t> > &contigs,
                             YAML::Emitter &yaml) {
    Status s;

    yaml << YAML::BeginSeq;
    for (auto& u_site : sites) {
        s = u_site.yaml(contigs, yaml);
        if (s.bad()) return s;
    }
    yaml << YAML::EndSeq;

    return Status::OK();
}

// Load the unified-sites from a file in yaml format.
//
Status unified_sites_of_yaml(const YAML::Node& yaml,
                             const std::vector<std::pair<std::string,size_t> > &contigs,
                             vector<unified_site> &sites) {
    Status s;
    if (!yaml.IsSequence()) {
        return Status::Invalid("not a sequence at top level");
    }

    sites.clear();
    for (YAML::const_iterator p = yaml.begin(); p != yaml.end(); ++p) {
        unified_site u_site(range(-1, -1, -1));
        s = unified_site::of_yaml(*p, contigs, u_site);
        if (s.bad()) return s;
        sites.push_back(u_site);
    }

    return Status::OK();
}

}}}

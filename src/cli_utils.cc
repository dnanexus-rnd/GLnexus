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

Status LoadYAMLFile(const string& filename, YAML::Node &node) {
    if (filename.size() == 0)
        return Status::Invalid("The YAML file must be specified");

    try {
        node = YAML::LoadFile(filename);
    } catch (exception &e) {
        return Status::Invalid("Error loading yaml file ", filename);
    }

    if (node.IsNull())
        return Status::NotFound("bad YAML file", filename);
    return Status::OK();
}

static void log_num_alleles(std::shared_ptr<spdlog::logger> logger,
                            const vector<discovered_alleles> &valleles,
                            const string& filename) {
    if (logger == nullptr)
        return;
    unsigned ct=0;
    for (const auto& dsals : valleles) {
        ct += dsals.size();
    }
    logger->info() << "loaded " << ct << " alleles from " << filename;
}

Status yaml_of_contigs_alleles(const vector<pair<string,size_t> > &contigs,
                               const vector<discovered_alleles> &valleles,
                               YAML::Emitter &yaml) {
    Status s;

    yaml << YAML::BeginMap;

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

    // Write alleles
    yaml << YAML::Key << "alleles";
    yaml << YAML::Value;
    {
        // Note: we might write out empty sites too
        yaml << YAML::BeginSeq;
        for (auto &dsal : valleles) {
            s = yaml_of_discovered_alleles(dsal, contigs, yaml);
            if (s.bad()) return s;
        }
        yaml << YAML::EndSeq;
    }

    yaml << YAML::EndMap;

    return Status::OK();
}

Status contigs_alleles_of_yaml(const YAML::Node& yaml,
                               std::vector<std::pair<std::string,size_t> > &contigs,
                               vector<discovered_alleles> &valleles) {
    Status s;
    if (!yaml.IsMap()) {
        return Status::Invalid("not a map at top level");
    }

    contigs.clear();
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
    auto n_ranges = yaml["alleles"];
    if (!n_ranges.IsSequence()) {
        return Status::Invalid("alleles should be a yaml sequence");
    }
    for (YAML::const_iterator p = n_ranges.begin(); p != n_ranges.end(); ++p) {
        discovered_alleles dsals;
        s = discovered_alleles_of_yaml((*p), contigs, dsals);
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


// Load first file
Status merge_discovered_allele_files(std::shared_ptr<spdlog::logger> logger,
                                     const vector<string> &filenames,
                                     vector<pair<string,size_t>> &contigs,
                                     vector<discovered_alleles> &valleles) {
    Status s;
    if (filenames.size() == 0)
        return Status::Invalid("no discovered allele files provided");
    contigs.clear();
    valleles.clear();

    // Load the first file
    YAML::Node node;
    const string& first_file = filenames[0];
    S(LoadYAMLFile(first_file, node));
    S(contigs_alleles_of_yaml(node, contigs, valleles));
    log_num_alleles(logger, valleles, first_file);

    if (filenames.size() == 1)
        return Status::OK();

    // Load the rest of the files, and merge
    for (auto it = filenames.begin() + 1; it != filenames.end(); ++it) {
        vector<pair<string,size_t>> contigs2;
        vector<discovered_alleles> valleles2;
        YAML::Node node;
        const string &crnt_file = *it;

        S(LoadYAMLFile(crnt_file, node));
        S(contigs_alleles_of_yaml(node, contigs2, valleles2));
        log_num_alleles(logger, valleles2, crnt_file);

        // sanity: verify that the contigs are the same, and that the number of sites
        // is the same.
        if (valleles.size() != valleles2.size())
            return Status::Invalid("The number of sites is different bewteen",
                                   first_file + " " + crnt_file);
        if (contigs.size() != contigs2.size())
            return Status::Invalid("The number of contigs is different bewteen",
                                   first_file + " " + crnt_file);
        for (int i=0; i < contigs.size(); ++i) {
            if (contigs[i] != contigs2[i])
                return Status::Invalid("The contigs are different bewteen",
                                       first_file + " " + crnt_file + " index=" + to_string(i));
        }

        for (int i=0; i < valleles2.size(); ++i) {
            S(merge_discovered_alleles(valleles2[i], valleles[i]));
        }
    }

    return Status::OK();
}

Status find_containing_range(const std::set<range> &ranges,
                             const range &pos,
                             range &ans) {
    // deal with special cases first
    int nm_elem = ranges.size();
    if (nm_elem == 0) {
        return Status::NotFound();
    }
    if (nm_elem == 1) {
        auto it = ranges.begin();
        if (it->contains(pos)) {
            ans = *it;
            return Status::OK();
        }
        return Status::NotFound();
    }

    // The returned value here is the first element that is
    // greater or equal to [pos].
    auto it = ranges.lower_bound(pos);
    if (it == ranges.end()) {
        it = std::prev(ranges.end());
    }

    if (it->contains(pos)) {
        // we got the right range
        ans = *it;
        return Status::OK();
    }

    // we landed one range after the one we need
    it = std::prev(it);
    if (it->contains(pos)) {
        ans = *it;
        return Status::OK();
    }
    return Status::NotFound();
}

}}}

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

Status yaml_of_contigs(const std::vector<std::pair<std::string,size_t> > &contigs,
                       YAML::Emitter &yaml) {
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

    return Status::OK();
}

Status contigs_of_yaml(const YAML::Node& yaml,
                       std::vector<std::pair<std::string,size_t> > &contigs) {
    contigs.clear();

    // read contigs
    if (!yaml.IsSequence()) {
        return Status::Invalid("contigs should be a yaml sequence");
    }
    for (auto p = yaml.begin(); p != yaml.end(); ++p) {
        const std::string name = (*p)["name"].as<std::string>();
        size_t size = (*p)["size"].as<size_t>();
        contigs.push_back(make_pair(name, size));
    }

    return Status::OK();
}

// We write the stream as follows:
//
// ---
// contigs
// ---
// allele 1
// ---
// allele 2
// ---
// etc.
// allele N
// ...
//
// Each element is transformed to YAML, and then written to the output stream.
// This generates the document in pieces, while keeping it valid YAML.
Status yaml_stream_of_discovered_alleles(const std::vector<std::pair<std::string,size_t> > &contigs,
                                         const discovered_alleles &dsals,
                                         std::ostream &os) {
    Status s;

    // write the contigs
    {
        os << "---" << endl;
        YAML::Emitter yaml;
        S(yaml_of_contigs(contigs, yaml));
        os << yaml.c_str() << endl;
    }

    // Write alleles
    for (auto &pr : dsals) {
        os << "---" << endl;
        YAML::Emitter yaml;
        S(yaml_of_one_discovered_allele(pr.first, pr.second, contigs, yaml));
        os << yaml.c_str() << endl;
    }

    // special notation for end-of-file
    os << "..." << endl;

    return Status::OK();
}


// In a YAML document built out of top-level sub-documents, get the next document.
//
// Notes:
// -  Catch any exceptions, and return an empty node.
static YAML::Node yaml_get_next_document(std::istream &is, bool skip_first_line = false) {
    try {
        char buf[200];
        stringstream ss;
        if (skip_first_line) {
            // the document starts with a "---", skip it
            is.getline(buf, 200);
        }

        while (is.good() && !is.eof()) {
            is.getline(buf, 200);
            size_t buf_len = strlen(buf);

            if (buf_len == 3) {
                // check if we reached the end of this top level document
                string marker(buf);
                if (marker == "---" || marker == "...")
                    break;
            }
            ss.write(buf, buf_len);
            ss.write("\n", 1);
        }

        // We have the entire document in memory. Convert to
        // a YAML node and return.
        return std::move(YAML::Load(ss.str()));
    } catch (exception e) {
        return YAML::Node();
    }
}


Status discovered_alleles_of_yaml_stream(std::istream &is,
                                         std::vector<std::pair<std::string,size_t> > &contigs,
                                         discovered_alleles &dsals) {
    Status s;
    contigs.clear();
    dsals.clear();

    // The first top-level document is the contigs
    YAML::Node doc = yaml_get_next_document(is, true);
    S(contigs_of_yaml(doc, contigs));

    // All other documents are discovered-alleles
    for (doc = yaml_get_next_document(is);
         !doc.IsNull();
         doc = yaml_get_next_document(is)) {
        allele allele(range(-1,-1,-1), "A");
        discovered_allele_info ainfo;

        S(one_discovered_allele_of_yaml(doc, contigs, allele, ainfo));
        dsals[allele] = ainfo;
    }
    if (contigs.size() == 0)
        return Status::Invalid("Empty contigs");
    if (dsals.size() == 0)
        return Status::Invalid("empty discovered alleles");

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
        S(u_site.yaml(contigs, yaml));
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
        S(unified_site::of_yaml(*p, contigs, u_site));
        sites.push_back(u_site);
    }

    return Status::OK();
}

// Load first file
Status merge_discovered_allele_files(std::shared_ptr<spdlog::logger> logger,
                                     const vector<string> &filenames,
                                     vector<pair<string,size_t>> &contigs,
                                     discovered_alleles &dsals) {
    Status s;
    if (filenames.size() == 0)
        return Status::Invalid("no discovered allele files provided");
    contigs.clear();
    dsals.clear();

    // Load the first file
    const string& first_file = filenames[0];
    {
        std::ifstream ifs(first_file.c_str());
        S(discovered_alleles_of_yaml_stream(ifs, contigs, dsals));
        logger->info() << "loaded " << dsals.size() << " alleles from " << first_file;
        ifs.close();
    }

    if (filenames.size() == 1)
        return Status::OK();

    // Load the rest of the files, and merge
    for (auto it = filenames.begin() + 1; it != filenames.end(); ++it) {
        vector<pair<string,size_t>> contigs2;
        discovered_alleles dsals2;
        const string &crnt_file = *it;
        std::ifstream ifs(crnt_file.c_str());

        S(discovered_alleles_of_yaml_stream(ifs, contigs2, dsals2));
        logger->info() << "loaded " << dsals2.size() << " alleles from " << crnt_file;
        ifs.close();

        // verify that the contigs are the same
        if (contigs != contigs2) {
            return Status::Invalid("The contigs are different between", first_file + " " + crnt_file);
        }

        S(merge_discovered_alleles(dsals2, dsals));
    }

    return Status::OK();
}

Status find_containing_range(const std::set<range> &ranges,
                             const range &pos,
                             range &ans) {
    if (ranges.size() == 0) {
        return Status::NotFound();
    }

    // The returned value here is the first element that is
    // greater or equal to [pos].
    auto it = ranges.lower_bound(pos);
    if (it == ranges.end() ||  !it->contains(pos)) {
        // we landed one range after the one we need
        if (it != ranges.begin()) {
            it = std::prev(it);
        }
    }

    if (it->contains(pos)) {
        // we got the right range
        ans = *it;
        return Status::OK();
    }
    return Status::NotFound();
}

}}}

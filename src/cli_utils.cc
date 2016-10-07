#include <iostream>
#include <exception>
#include <fstream>
#include <sstream>
#include <regex>
#include "cli_utils.h"
#include "ctpl_stl.h"

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
// N, contigs
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
Status yaml_stream_of_discovered_alleles(unsigned N, const std::vector<std::pair<std::string,size_t> > &contigs,
                                         const discovered_alleles &dsals,
                                         std::ostream &os) {
    Status s;

    // write the contigs
    {
        os << "---" << endl;
        YAML::Emitter yaml;
        yaml << YAML::BeginMap;
        yaml << YAML::Key << "N" << YAML::Value << N;
        yaml << YAML::Key << "contigs" << YAML::Value;
        S(yaml_of_contigs(contigs, yaml));
        yaml << YAML::EndMap;
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


static string yaml_begin_doc = "---";
static string yaml_end_doc_list = "...";

// Verify the document begins with '---' line. Move the stream position to the next line.
static Status yaml_verify_begin_doc_list(std::istream &is) {
    try {
        string marker;
        std::getline(is, marker);
        if (marker != yaml_begin_doc)
            return Status::Invalid("document does not start with `---`");
        return Status::OK();
    } catch (exception e) {
        string err = e.what();
        return Status::Failure("exception caught in yaml_verify_begin_doc_list: ", err);
    }
}

// Verify that we have reached the end of the file
static Status yaml_verify_eof(std::istream &is) {
    try {
        // The end-of-file flag is lit up only after reading past the
        // end of file
        char c;
        is.get(c);
        if (!is.eof()) {
            return Status::Invalid("Found YAML end-of-document marker, but there is unread data");
        }
        return Status::OK();
    } catch (exception e) {
        string err = e.what();
        return Status::Failure("exception caught in yaml_verify_eof: ", err);
    }
}

// In a YAML document built as a of document list, get the next document.
//
// Notes:
// -  Catch any exceptions, and convert to bad status
static Status yaml_get_next_document(std::istream &is,
                                     YAML::Node &ans,
                                     string& next_marker) {
    try {
        stringstream ss;

        while (is.good() && !is.eof()) {
            string line;
            std::getline(is, line);
            if (line.size() == 3) {
                // check if we reached the end of this top level document
                if (line == yaml_begin_doc ||
                    line == yaml_end_doc_list) {
                    // We have the entire document in memory. Convert to
                    // a YAML node and return.
                    next_marker = line;
                    ans = std::move(YAML::Load(ss.str()));
                    return Status::OK();
                }
            }
            ss.write(line.c_str(), line.size());
            ss.write("\n", 1);
        }

        if (!is.good())
            return Status::IOError("reading yaml stream");
        return Status::Invalid("premature end of document, did not find end-of-document marker");
    } catch (exception e) {
        string err = e.what();
        return Status::Failure("exception caught in yaml_get_next_document: ", err);
    }
}


Status discovered_alleles_of_yaml_stream(std::istream &is,
                                         unsigned &N, std::vector<std::pair<std::string,size_t> > &contigs,
                                         discovered_alleles &dsals) {
    Status s;
    string next_marker;

    contigs.clear();
    dsals.clear();
    S(yaml_verify_begin_doc_list(is));

    // The first top-level document has N & contigs
    {
        YAML::Node doc;
        S(yaml_get_next_document(is, doc, next_marker));

        if (!doc.IsMap()) return Status::Invalid("discovered alleles header missing");
        if (!doc["N"] || !doc["N"].IsScalar()) return Status::Invalid("discovered alleles header missing N");
        N = doc["N"].as<unsigned>();

        if (!doc["contigs"] || !doc["contigs"].IsSequence()) return Status::Invalid("discovered alleles header missing contigs");
        S(contigs_of_yaml(doc["contigs"], contigs));
    }

    // All other documents are discovered-alleles
    while (next_marker != yaml_end_doc_list) {
        YAML::Node doc;
        S(yaml_get_next_document(is, doc, next_marker));

        allele allele(range(-1,-1,-1), "A");
        discovered_allele_info ainfo;

        S(one_discovered_allele_of_yaml(doc, contigs, allele, ainfo));
        dsals[allele] = ainfo;
    }
    S(yaml_verify_eof(is));

    if (contigs.size() == 0)
        return Status::Invalid("Empty contigs");
    return Status::OK();
}

// Serialize the unified sites to yaml format.
//
Status yaml_stream_of_unified_sites(const std::vector<unified_site> &sites,
                                    const std::vector<std::pair<std::string,size_t> > &contigs,
                                    std::ostream &os) {
    Status s;
    for (auto& u_site : sites) {
        os << "---" << endl;
        YAML::Emitter yaml;
        S(u_site.yaml(contigs, yaml));
        os << yaml.c_str() << endl;
    }

    // special notation for end-of-file
    os << "..." << endl;
    return Status::OK();
}

// Load the unified-sites from a file in yaml format.
//
Status unified_sites_of_yaml_stream(std::istream &is,
                                    const std::vector<std::pair<std::string,size_t> > &contigs,
                                    std::vector<unified_site> &sites) {
    Status s;
    string next_marker;
    sites.clear();

    S(yaml_verify_begin_doc_list(is));

    // Read a list of documents representing unified-sites
    do {
        YAML::Node doc;
        S(yaml_get_next_document(is, doc, next_marker));
        unified_site u_site(range(-1, -1, -1));
        S(unified_site::of_yaml(doc, contigs, u_site));
        sites.push_back(u_site);
    } while (next_marker != yaml_end_doc_list);

    S(yaml_verify_eof(is));
    return Status::OK();
}

// Load first file
Status merge_discovered_allele_files(std::shared_ptr<spdlog::logger> logger,
                                     size_t nr_threads,
                                     const vector<string> &filenames,
                                     unsigned &N, vector<pair<string,size_t>> &contigs,
                                     discovered_alleles &dsals) {
    Status s;
    mutex mu;

    if (nr_threads == 0) {
        // by default, we use the number of cores for the thread pool size.
        nr_threads = thread::hardware_concurrency();
    }
    ctpl::thread_pool threadpool(nr_threads);

    if (filenames.size() == 0)
        return Status::Invalid("no discovered allele files provided");
    contigs.clear();
    dsals.clear();
    N = 0;

    // Load the files in parallel, and merge
    vector<future<GLnexus::Status>> statuses;
    for (const auto& dsal_file: filenames) {
        auto fut = threadpool.push([&, dsal_file](int tid){
            discovered_alleles dsals2;
            vector<pair<string,size_t>> contigs2;
            unsigned int N2;

            //std::ifstream ifs(dsal_file.c_str());
            //Status s = discovered_alleles_of_yaml_stream(ifs, N2, contigs2, dsals2);
            //ifs.close();
            Status s = discovered_alleles_of_capnp(dsal_file, N2, contigs2, dsals2);
            if (!s.ok()) {
                logger->info() << "Error loading alleles from " << dsal_file;
                return s;
            }
            logger->info() << "loaded " << dsals2.size() << " alleles from " << dsal_file << ", N = " << N2;

            lock_guard<mutex> lock(mu);
            if (contigs.empty()) {
                // This is the first file, initialize the result data-structures
                contigs = move(contigs2);
                dsals = move(dsals2);
                N = N2;
                return Status::OK();
            }

            // We have already loaded and merged some files.
            // Verify that the contigs are the same
            if (contigs != contigs2) {
                return Status::Invalid("The contigs do not match");
            }

            // Merge the discovered alleles
            N += N2;
            return merge_discovered_alleles(dsals2, dsals);
        });
        statuses.push_back(move(fut));
    }

    // collect results, wait for all threads to complete operation
    vector<GLnexus::Status> failures;
    for (size_t i = 0; i < filenames.size(); i++) {
        GLnexus::Status s_i(move(statuses[i].get()));
        if (!s_i.ok()) {
            failures.push_back(move(s_i));
        }
    }
    if (!failures.empty())
        return move(failures[0]);

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

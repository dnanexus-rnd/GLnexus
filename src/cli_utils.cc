#include "cli_utils.h"
#include "ctpl_stl.h"
#include <exception>
#include <fts.h>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "crc32c.h"
#include "service.h"
#include "compare_queries.h"

#include "BCFKeyValueData.h"

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

Status parse_bed_file(std::shared_ptr<spdlog::logger> logger,
                      const string &bedfilename,
                      const vector<pair<string,size_t> > &contigs,
                      vector<range> &ranges) {
    if (bedfilename.empty()) {
        return Status::Invalid("Empty bed file");
    }
    if (!check_file_exists(bedfilename)) {
        return Status::IOError("bed file does not exist", bedfilename);
    }

    // read BED file
    string rname, beg_txt, end_txt;
    ifstream bedfile(bedfilename);
    while (bedfile >> rname >> beg_txt >> end_txt) {
        int rid = 0;
        for(; rid<contigs.size(); rid++)
            if (contigs[rid].first == rname)
                break;
        if (rid == contigs.size()) {
            return Status::Invalid("Unknown contig ", rname);
        }
        ranges.push_back(range(rid,
                               strtol(beg_txt.c_str(), nullptr, 10),
                               strtol(end_txt.c_str(), nullptr, 10)));
    }
    if (bedfile.bad() || !bedfile.eof()) {
        return Status::IOError( "Error reading ", bedfilename);
    }

    sort(ranges.begin(), ranges.end());
    for (auto& query : ranges) {
        if (query.beg < 0 || query.end < 1 || query.end <= query.beg) {
            return Status::Invalid("query range ", query.str(contigs));
        }
        if (query.end > contigs[query.rid].second) {
            query.end = contigs[query.rid].second;
            logger->warn() << "Truncated query range at end of contig: " << query.str(contigs);
        }
    }

    // make an orderly pass on the ranges, and verify that there are no overlaps.
    if (ranges.size() > 1) {
        const range &prev = *(ranges.begin());
        for (auto it = ranges.begin() + 1; it != ranges.end(); ++it) {
            if (prev.overlaps(*it))
                return Status::Invalid("overlapping ranges ", prev.str(contigs) + " " + it->str(contigs));
        }
    }

    return Status::OK();
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

// Write the discovered alleles to a file
Status yaml_write_discovered_alleles_to_file(const discovered_alleles &dsals,
                                             const vector<pair<string,size_t>> &contigs,
                                             unsigned int sample_count,
                                             const string &filename) {
    Status s;

    ofstream ofs(filename, std::ofstream::out | std::ofstream::trunc);
    if (ofs.bad())
        return Status::IOError("could not open file for writing", filename);
    S(yaml_stream_of_discovered_alleles(sample_count, contigs, dsals, ofs));
    ofs.close();

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

// Write the unified-sites to a file
Status write_unified_sites_to_file(const vector<unified_site> &sites,
                                   const vector<pair<string,size_t>> &contigs,
                                   const string &filename) {
    Status s;

    ofstream ofs(filename, std::ofstream::out | std::ofstream::trunc);
    if (ofs.bad())
        return Status::IOError("could not open file for writing", filename);

    S(utils::yaml_stream_of_unified_sites(sites, contigs, ofs));
    ofs.close();

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

// Merge a bunch of discovered-allele files, all in capnp format.
Status merge_discovered_allele_files(std::shared_ptr<spdlog::logger> logger,
                                     size_t nr_threads,
                                     const vector<string> &filenames,
                                     unsigned &N, vector<pair<string,size_t>> &contigs,
                                     discovered_alleles& dsals) {
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
    vector<future<Status>> statuses;
    for (const auto& dsal_file: filenames) {
        auto fut = threadpool.push([&, dsal_file](int tid){
            // load one file
            discovered_alleles dsals2;
            vector<pair<string,size_t>> contigs2;
            unsigned int N2;

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
            S(merge_discovered_alleles(dsals2, dsals));
            return Status::OK();
        });
        statuses.push_back(move(fut));
    }

    // collect results, wait for all threads to complete operation
    vector<Status> failures;
    for (size_t i = 0; i < filenames.size(); i++) {
        Status s_i(move(statuses[i].get()));
        if (!s_i.ok()) {
            failures.push_back(move(s_i));
        }
    }
    if (!failures.empty())
        return move(failures[0]);

    return Status::OK();
}

// Check if a file exists
bool check_file_exists(const string &filename) {
    ifstream ifs(filename);
    if (!ifs) {
        return false;
    }
    return true;
}

bool check_dir_exists(const string &path) {
    struct stat info;

    if (stat(path.c_str(), &info) != 0) {
        // Could not access the directory
        return false;
    }
    if (info.st_mode & S_IFDIR) {
        return true;
    }
    return false;
}


// hard-coded configuration presets for unifier & genotyper. TODO: these
// should reside in some user-modifiable yml file
static const char* config_presets_yml = R"eof(
gatk:
    description: Merge and joint-call GATK-style gVCFs
    unifier_config:
        min_AQ1: 70
        min_AQ2: 40
        min_GQ: 40
        monoallelic_sites_for_lost_alleles: true
    genotyper_config:
        required_dp: 1
        revise_genotypes: true
        liftover_fields:
            - orig_names: [GQ]
              name: GQ
              description: '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">'
              type: int
              number: basic
              combi_method: min
              count: 1
              ignore_non_variants: true
            - orig_names: [DP, MIN_DP]
              name: DP
              description: '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">'
              type: int
              combi_method: min
              number: basic
              count: 1
            - orig_names: [AD]
              name: AD
              description: '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">'
              type: int
              number: alleles
              combi_method: min
              default_type: zero
              count: 0
            - orig_names: [SB]
              name: SB
              description: '##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fishers Exact Test to detect strand bias.">'
              type: int
              combi_method: missing
              number: basic
              count: 4
gatk_unfiltered:
    description: Merge GATK-style gVCFs with no filtering or genotype revision. Without any allele quality thresholds, large-cohort pVCFs may contain excessive false positive sites.
    unifier_config:
        min_AQ1: 0
        min_AQ2: 0
        min_GQ: 0
        monoallelic_sites_for_lost_alleles: true
    genotyper_config:
        required_dp: 1
        revise_genotypes: false
        liftover_fields:
            - orig_names: [GQ]
              name: GQ
              description: '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">'
              type: int
              number: basic
              combi_method: min
              count: 1
              ignore_non_variants: true
            - orig_names: [PL]
              name: PL
              description: '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype Likelihoods">'
              type: int
              number: genotype
              combi_method: missing
              count: 0
              ignore_non_variants: true
            - orig_names: [DP, MIN_DP]
              name: DP
              description: '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">'
              type: int
              combi_method: min
              number: basic
              count: 1
            - orig_names: [AD]
              name: AD
              description: '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">'
              type: int
              number: alleles
              combi_method: min
              default_type: zero
              count: 0
            - orig_names: [SB]
              name: SB
              description: '##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fishers Exact Test to detect strand bias.">'
              type: int
              combi_method: missing
              number: basic
              count: 4
xAtlas:
    unifier_config:
        drop_filtered: true
        monoallelic_sites_for_lost_alleles: true
    genotyper_config:
        required_dp: 0
        allow_partial_data: true
        revise_genotypes: false
        # TODO: ref_dp_format=DPX[0] would be more precise
        ref_dp_format: DP
        liftover_fields:
            - orig_names: [GQ]
              name: GQ
              description: '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">'
              type: int
              number: basic
              combi_method: min
              count: 1
              ignore_non_variants: true
            - orig_names: [PL]
              name: PL
              description: '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype Likelihoods">'
              type: int
              number: genotype
              combi_method: missing
              count: 0
              ignore_non_variants: true
            - orig_names: [DP]
              name: DP
              description: '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">'
              type: int
              combi_method: min
              number: basic
              count: 1
            - orig_names: [RR]
              name: RR
              description: '##FORMAT=<ID=RR,Number=1,Type=Integer,Description="Reference Read Depth">'
              type: int
              combi_method: min
              number: basic
              count: 1
            - orig_names: [VR]
              name: VR
              description: '##FORMAT=<ID=VR,Number=1,Type=Integer,Description="Major Variant Read Depth">'
              type: int
              combi_method: min
              number: basic
              count: 1
              ignore_non_variants: true
            - orig_names: [FILTER]
              name: FT
              description: '##FORMAT=<ID=FT,Number=1,Type=String,Description="FILTER field from sample gVCF (other than PASS)">'
              type: string
              combi_method: missing
              number: basic
              count: 1
              ignore_non_variants: true
weCall:
    unifier_config:
        min_AQ1: 60
        min_AQ2: 30
        min_GQ: 30
        monoallelic_sites_for_lost_alleles: true
    genotyper_config:
        required_dp: 0
        revise_genotypes: false
        liftover_fields:
            - orig_names: [GQ]
              name: GQ
              description: '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">'
              type: int
              number: basic
              combi_method: min
              count: 1
              ignore_non_variants: true
            - orig_names: [MIN_DP, DP]
              name: DP
              description: '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">'
              type: int
              combi_method: min
              number: basic
              count: 1
            - orig_names: [AD]
              name: AD
              description: '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">'
              type: int
              number: alleles
              combi_method: min
              default_type: zero
              count: 0
            - orig_names: [FILTER]
              name: FT
              description: '##FORMAT=<ID=FT,Number=1,Type=String,Description="FILTER field from sample gVCF">'
              type: string
              combi_method: missing
              number: basic
              count: 1
              ignore_non_variants: true
DeepVariant:
    unifier_config:
        min_AQ1: 0
        min_AQ2: 0
        min_GQ: 0
        monoallelic_sites_for_lost_alleles: true
    genotyper_config:
        required_dp: 0
        revise_genotypes: false
        liftover_fields:
            - orig_names: [GQ]
              name: GQ
              description: '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">'
              type: int
              number: basic
              combi_method: min
              count: 1
              ignore_non_variants: true
            - orig_names: [PL]
              name: PL
              description: '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype Likelihoods">'
              type: int
              number: genotype
              combi_method: missing
              count: 0
              ignore_non_variants: true
            - orig_names: [DP, MIN_DP]
              name: DP
              description: '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">'
              type: int
              combi_method: min
              number: basic
              count: 1
            - orig_names: [AD]
              name: AD
              description: '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">'
              type: int
              number: alleles
              combi_method: min
              default_type: zero
              count: 0
)eof";

Status load_config(std::shared_ptr<spdlog::logger> logger,
                   const YAML::Node& config,
                   unifier_config& unifier_cfg,
                   genotyper_config& genotyper_cfg,
                   std::string& config_crc32c) {
    Status s;
    if (config["unifier_config"]) {
        S(unifier_config::of_yaml(config["unifier_config"], unifier_cfg));
    }
    if (config["genotyper_config"]) {
        S(genotyper_config::of_yaml(config["genotyper_config"], genotyper_cfg));
    }
    YAML::Emitter em;
    em << YAML::BeginMap
       << YAML::Key << "unifier_config" << YAML::Value;
    S(unifier_cfg.yaml(em));
    em << YAML::Key << "genotyper_config" << YAML::Value;
    S(genotyper_cfg.yaml(em));
    em << YAML::EndMap;
    std::string config_text = em.c_str();
    logger->info() << "config:\n" << config_text;

    config_crc32c = std::to_string(rocksdb::crc32c::Value(config_text.c_str(), config_text.size()));
    logger->info() << "config CRC32C = " << config_crc32c;
    return Status::OK();
}

Status load_config(std::shared_ptr<spdlog::logger> logger,
                   const std::string& name,
                   unifier_config& unifier_cfg,
                   genotyper_config& genotyper_cfg,
                   std::string& config_crc32c) {
    Status s;
    if (name.size() > 4 && name.substr(name.size() - 4) == ".yml") {
        try {
            logger->info() << "Loading config YAML file " << name;
            YAML::Node config = YAML::LoadFile(name);
            if (!config || !config.IsMap()) {
                return Status::IOError("loading configuration YAML file", name);
            }
            return load_config(logger, config, unifier_cfg, genotyper_cfg, config_crc32c);
        } catch (YAML::Exception& exn) {
            return Status::IOError("loading configuration YAML file", name);
        }
    } else {
        logger->info() << "Loading config preset " << name;
        YAML::Node presets = YAML::Load(config_presets_yml);
        if (!presets || !presets.IsMap() || !presets[name] || !presets[name].IsMap()) {
            return Status::NotFound("unknown configuration preset", name);
        }
        return load_config(logger, presets[name], unifier_cfg, genotyper_cfg, config_crc32c);
    }
}

RocksKeyValue::prefix_spec* GLnexus_prefix_spec() {
    static unique_ptr<RocksKeyValue::prefix_spec> p;
    if (!p) {
        p = make_unique<RocksKeyValue::prefix_spec>("bcf", BCFKeyValueDataPrefixLength());
    }
    return p.get();
}


// Initialize a database
Status db_init(std::shared_ptr<spdlog::logger> logger,
               const string &dbpath,
               const string &exemplar_gvcf,
               vector<pair<string,size_t>> &contigs,
               size_t bucket_size) {
    Status s;
    logger->info() << "init database, exemplar_vcf=" << exemplar_gvcf;
    if (check_dir_exists(dbpath)) {
        return Status::IOError("Database directory already exists", dbpath);
    }

    // load exemplar contigs
    unique_ptr<vcfFile, void(*)(vcfFile*)> vcf(bcf_open(exemplar_gvcf.c_str(), "r"),
                                               [](vcfFile* f) { bcf_close(f); });
    if (!vcf) {
        return Status::IOError("Failed to open exemplar gVCF file at ", exemplar_gvcf);
    }
    unique_ptr<bcf_hdr_t, void(*)(bcf_hdr_t*)> hdr(bcf_hdr_read(vcf.get()), &bcf_hdr_destroy);
    if (!hdr) {
        return Status::IOError("Failed to read gVCF file header from", exemplar_gvcf);
    }
    int ncontigs = 0;
    const char **contignames = bcf_hdr_seqnames(hdr.get(), &ncontigs);
    for (int i = 0; i < ncontigs; i++) {
        if (hdr->id[BCF_DT_CTG][i].val == nullptr) {
            return Status::Invalid("Invalid gVCF header in ", exemplar_gvcf);
        }
        contigs.push_back(make_pair(string(contignames[i]),
                                    hdr->id[BCF_DT_CTG][i].val->info[0]));
    }
    free(contignames);

    // create and initialize the database
    unique_ptr<KeyValue::DB> db;
    S(RocksKeyValue::Initialize(dbpath, db, GLnexus_prefix_spec()));
    S(BCFKeyValueData::InitializeDB(db.get(), contigs, bucket_size));

    // report success
    logger->info() << "Initialized GLnexus database in " << dbpath;
    logger->info() << "bucket size: " << bucket_size;

    stringstream ss;
    ss << "contigs:";
    for (const auto& contig : contigs) {
        ss << " " << contig.first;
    }
    logger->info() << ss.str();
    S(db->flush());
    db.reset();

    return Status::OK();
}


Status db_get_contigs(std::shared_ptr<spdlog::logger> logger,
                      const string &dbpath,
                      std::vector<std::pair<std::string,size_t> > &contigs) {
    Status s;
    logger->info() << "db_get_contigs " << dbpath;

    unique_ptr<KeyValue::DB> db;

    S(RocksKeyValue::Open(dbpath, db, GLnexus_prefix_spec(),
                          RocksKeyValue::OpenMode::READ_ONLY));
    {
        unique_ptr<BCFKeyValueData> data;
        S(BCFKeyValueData::Open(db.get(), data));
        S(data->contigs(contigs));
    }

    return Status::OK();
}

Status db_bulk_load(std::shared_ptr<spdlog::logger> logger,
                    size_t nr_threads,
                    const vector<string> &gvcfs,
                    const string &dbpath,
                    const vector<range> &ranges_i,
                    std::vector<std::pair<std::string,size_t> > &contigs, // output param
                    bool delete_gvcf_after_load) {
    Status s;

    set<range> ranges;
    for (auto &r : ranges_i)
        ranges.insert(r);

    // open the database
    unique_ptr<KeyValue::DB> db;
    S(RocksKeyValue::Open(dbpath, db, GLnexus_prefix_spec(),
                          RocksKeyValue::OpenMode::BULK_LOAD));
    unique_ptr<BCFKeyValueData> data;
    S(BCFKeyValueData::Open(db.get(), data));

    unique_ptr<MetadataCache> metadata;
    S(MetadataCache::Start(*data, metadata));
    contigs = metadata->contigs();

     if (ranges.size()) {
        ostringstream ss;
        for (const auto& rng : ranges) {
            ss << " " << rng.str(contigs);
        }
        logger->info() << "Beginning bulk load of records overlapping:" << ss.str();
    } else {
        logger->info() << "Beginning bulk load with no range filter.";
    }

    ctpl::thread_pool threadpool(nr_threads);
    vector<future<Status>> statuses;
    set<string> datasets_loaded;
    BCFKeyValueData::import_result stats;
    mutex mu;
    string dataset;

    // load the gVCFs on the thread pool
    for (const auto& gvcf : gvcfs) {
        // infer dataset name as the gVCF filename minus path and extension
        size_t p = gvcf.find_last_of('/');
        if (p != string::npos && p < gvcf.size()-1) {
            dataset = gvcf.substr(p+1);
        } else {
            dataset = gvcf;
        }
        for (const string& ext : {".bgzip",".gz",".gvcf",".g.vcf",".vcf"}) {
            if (dataset.size() > ext.size()) {
                p = dataset.size() - ext.size();
                if (dataset.rfind(ext) == p) {
                    dataset.erase(p);
                }
            }
        }

        auto fut = threadpool.push([&, gvcf, dataset](int tid) {
                BCFKeyValueData::import_result rslt;
                Status ls = data->import_gvcf(*metadata, dataset, gvcf, ranges, rslt);
                if (ls.ok()) {
                    if (delete_gvcf_after_load && unlink(gvcf.c_str())) {
                        logger->warn() << "Loaded " << gvcf << " successfully, but failed deleting it afterwards.";
                    }
                    lock_guard<mutex> lock(mu);
                    if (rslt.records == 0) {
                        logger->warn() << "No data loaded from " << gvcf << " after range and other filters.";
                    }
                    stats += rslt;
                    datasets_loaded.insert(dataset);
                    size_t n = datasets_loaded.size();
                    if (n % 100 == 0) {
                        logger->info() << n << " (" << dataset << ")...";
                    }
                }
                return ls;
            });
        statuses.push_back(move(fut));
        dataset.clear();
    }

    // collect results
    vector<pair<string,Status>> failures;
    for (size_t i = 0; i < gvcfs.size(); i++) {
        Status s_i(move(statuses[i].get()));
        if (!s_i.ok()) {
            failures.push_back(make_pair(gvcfs[i],move(s_i)));
        }
    }

    // report results
    logger->info() << "Loaded " << datasets_loaded.size() << " datasets with "
                    << stats.samples.size() << " samples; "
                    << stats.bytes << " bytes in "
                    << stats.records << " BCF records ("
                    << stats.duplicate_records << " duplicate) in "
                    << stats.buckets << " buckets. "
                    << "Bucket max " << stats.max_bytes << " bytes, max "
                    << stats.max_records << " records. "
                    << stats.skipped_records << " BCF records skipped due to caller-specific exceptions.";

    // call all_samples_sampleset to create the sample set including
    // the newly loaded ones. By doing this now we make it possible
    // for other CLI functions to open the database in purely read-
    // only mode (since the sample set has to get written into the
    // database to be used)
    string sampleset;
    S(data->all_samples_sampleset(sampleset));
    logger->info() << "Created sample set " << sampleset;

    if (failures.size()) {
        for (const auto& p : failures) {
            logger->error() << p.first << " " << p.second.str();
        }
        return Status::Failure("FAILED to load ", failures.size() + " datasets:");
    }

    logger->info() << "Flushing and compacting database...";
    S(db->flush());
    db.reset();
    logger->info() << "Bulk load complete!";
    return Status::OK();
}

Status discover_alleles(std::shared_ptr<spdlog::logger> logger,
                        size_t nr_threads,
                        const string &dbpath,
                        const vector<range> &ranges,
                        const std::vector<std::pair<std::string,size_t> > &contigs,
                        discovered_alleles &dsals,
                        unsigned &sample_count) {
    Status s;
    unique_ptr<KeyValue::DB> db;
    unique_ptr<BCFKeyValueData> data;
    dsals.clear();

    // open the database in read-only mode
    S(RocksKeyValue::Open(dbpath, db, GLnexus_prefix_spec(),
                          RocksKeyValue::OpenMode::READ_ONLY));
    S(BCFKeyValueData::Open(db.get(), data));

    // start service, discover alleles
    service_config svccfg;
    svccfg.threads = nr_threads;
    unique_ptr<Service> svc;
    S(Service::Start(svccfg, *data, *data, svc));

    string sampleset;
    S(data->all_samples_sampleset(sampleset));
    logger->info() << "found sample set " << sampleset;

    logger->info() << "discovering alleles in " << ranges.size() << " range(s)";
    vector<discovered_alleles> valleles;
    S(svc->discover_alleles(sampleset, ranges, sample_count, valleles));

    for (auto it = valleles.begin(); it != valleles.end(); ++it) {
        S(merge_discovered_alleles(*it, dsals));
        it->clear(); // free some memory
    }
    logger->info() << "discovered " << dsals.size() << " alleles";
    return Status::OK();
}

Status unify_sites(std::shared_ptr<spdlog::logger> logger,
                   const unifier_config &unifier_cfg,
                   const vector<pair<string,size_t> > &contigs,
                   discovered_alleles &dsals,
                   unsigned sample_count,
                   vector<unified_site> &sites,
                   unifier_stats& stats) {
    Status s;
    S(unified_sites(unifier_cfg, sample_count, dsals, sites, stats));

    // sanity check, sites are in-order
    if (sites.size() > 1) {
        auto p = sites.begin();
        for (auto q = p+1; q != sites.end(); ++p, ++q) {
            if (q->pos < p->pos) {
                return Status::Failure(
                    "BUG: unified sites failed sanity check -- sites are out of order",
                    p->pos.str(contigs)  + " " + q->pos.str(contigs));
            }
        }
    }

    return Status::OK();
}


Status genotype(std::shared_ptr<spdlog::logger> logger,
                size_t nr_threads,
                const string &dbpath,
                const genotyper_config &genotyper_cfg,
                const vector<unified_site> &sites,
                const vector<string>& extra_header_lines,
                const string &output_filename) {
    Status s;
    logger->info() << "Lifting over " << genotyper_cfg.liftover_fields.size() << " fields.";

    // open the database in read-only mode
    unique_ptr<KeyValue::DB> db;
    S(RocksKeyValue::Open(dbpath, db, GLnexus_prefix_spec(),
                                   RocksKeyValue::OpenMode::READ_ONLY));
    unique_ptr<BCFKeyValueData> data;
    S(BCFKeyValueData::Open(db.get(), data));

    std::vector<std::pair<std::string,size_t> > contigs;
    S(data->contigs(contigs));

    // start service, discover alleles, unify sites, genotype sites
    service_config svccfg;
    svccfg.threads = nr_threads;
    svccfg.extra_header_lines = extra_header_lines;
    unique_ptr<Service> svc;
    S(Service::Start(svccfg, *data, *data, svc));

    string sampleset;
    S(data->all_samples_sampleset(sampleset));
    logger->info() << "found sample set " << sampleset;

    S(svc->genotype_sites(genotyper_cfg, sampleset, sites, output_filename));
    logger->info() << "genotyping complete!";

    auto stalls_ms = svc->threads_stalled_ms();
    if (stalls_ms) {
        logger->info() << "worker threads were cumulatively stalled for " << stalls_ms << "ms";
    }

    std::shared_ptr<StatsRangeQuery> statsRq = data->getRangeStats();
    logger->info() << statsRq->str();

    return Status::OK();
}

Status compare_db_itertion_algorithms(std::shared_ptr<spdlog::logger> logger,
                                      const std::string &dbpath,
                                      int n_iter) {
    Status s;
    unique_ptr<KeyValue::DB> db;
    unique_ptr<BCFKeyValueData> data;
    string sampleset;
    S(RocksKeyValue::Open(dbpath, db, GLnexus_prefix_spec(),
                          RocksKeyValue::OpenMode::READ_ONLY));
    S(BCFKeyValueData::Open(db.get(), data));

    unique_ptr<MetadataCache> metadata;
    S(MetadataCache::Start(*data, metadata));

    S(data->all_samples_sampleset(sampleset));
    logger->info() << "using sample set " << sampleset;

    // get samples and datasets
    shared_ptr<const set<string>> samples, datasets;
    S(metadata->sampleset_datasets(sampleset, samples, datasets));

    // compare queries
    S(compare_queries::compare_n_queries(n_iter, *data, *metadata, sampleset));
    return Status::OK();
}

}}}

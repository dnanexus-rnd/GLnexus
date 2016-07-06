#ifndef GLNEXUS_APPLET_UTILS_H
#define GLNEXUS_APPLET_UTILS_H

#include <string>
#include <vector>
#include <map>
#include "types.h"
#include "spdlog/spdlog.h"

namespace GLnexus {
namespace cli {
namespace utils {

// contig, in which case it gets mapped to the contig's full length.
bool parse_range(const std::vector<std::pair<std::string,size_t> >& contigs,
                 const std::string& range_txt, range& ans);

// parse a comma-separated list of ranges
bool parse_ranges(const std::vector<std::pair<std::string,size_t> >& contigs,
                  const std::string& ranges, std::vector<range>& ans);

// Read from disk and parse a YAML file
Status LoadYAMLFile(const std::string& filename, YAML::Node &node);

// Serialize the contigs, and discovered alleles to YAML
Status yaml_of_contigs_alleles(const std::vector<std::pair<std::string,size_t> > &contigs,
                               const discovered_alleles &dsals,
                               YAML::Emitter &yaml);

// Load a YAML file, previously created with the above function
Status contigs_alleles_of_yaml(const YAML::Node& yaml,
                               std::vector<std::pair<std::string,size_t> > &contigs,
                               discovered_alleles &dsals);

// Serialize a vector of unified-alleles to YAML.
Status yaml_of_unified_sites(const std::vector<unified_site> &sites,
                             const std::vector<std::pair<std::string,size_t> > &contigs,
                             YAML::Emitter &yaml);

// Load from a file, data previously serialized with the above function
Status unified_sites_of_yaml(const YAML::Node& yaml,
                             const std::vector<std::pair<std::string,size_t> > &contigs,
                             std::vector<unified_site> &sites);

Status merge_all(const std::vector<discovered_alleles> &valleles,
                 discovered_alleles &dsals);

// Merge a bunch of discovered-allele files, all in YAML format.
//
Status merge_discovered_allele_files(std::shared_ptr<spdlog::logger> logger,
                                     const std::vector<std::string> &filenames,
                                     std::vector<std::pair<std::string,size_t>> &contigs,
                                     discovered_alleles &dsals);


// Find which range contains [pos].
Status find_containing_range(const std::set<range> &ranges,
                             const range &pos,
                             range &ans);

}}}

#endif

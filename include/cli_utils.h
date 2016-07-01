#ifndef GLNEXUS_APPLET_UTILS_H
#define GLNEXUS_APPLET_UTILS_H

#include <string>
#include <vector>
#include <map>
#include "types.h"

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

// Serialize the contigs, ranges, and discovered alleles to YAML
Status yaml_of_contigs_alleles_ranges(const std::vector<std::pair<std::string,size_t> > &contigs,
                                      const std::vector<range> &ranges,
                                      const std::vector<discovered_alleles> &valleles,
                                      YAML::Emitter &yaml);

// Load a YAML file, previously created with the above function
Status contigs_alleles_ranges_of_yaml(const YAML::Node& yaml,
                                      std::vector<std::pair<std::string,size_t> > &contigs,
                                      std::vector<range> &ranges,
                                      std::vector<discovered_alleles> &valleles);

// Serialize a vector of unified-alleles to YAML.
Status yaml_of_unified_sites(const std::vector<unified_site> &sites,
                             const std::vector<std::pair<std::string,size_t> > &contigs,
                             YAML::Emitter &yaml);

// Load from a file, data previously serialized with the above function
Status unified_sites_of_yaml(const YAML::Node& yaml,
                             const std::vector<std::pair<std::string,size_t> > &contigs,
                             std::vector<unified_site> &sites);

Status merge_discovered_allele_files(const std::vector<std::string> &filenames,
                                     std::vector<std::pair<std::string,size_t>> &contigs,
                                     std::vector<range> &ranges,
                                     std::vector<discovered_alleles> &valleles);
}}}

#endif

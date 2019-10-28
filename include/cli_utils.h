#ifndef GLNEXUS_APPLET_UTILS_H
#define GLNEXUS_APPLET_UTILS_H

#include <string>
#include <vector>
#include <map>
#include "types.h"
#include "spdlog/spdlog.h"
#include "RocksKeyValue.h"
#include "BCFKeyValueData.h"
#include "unifier.h"

namespace GLnexus {
namespace cli {
namespace utils {

bool detect_jemalloc(std::shared_ptr<spdlog::logger> logger);

// contig, in which case it gets mapped to the contig's full length.
bool parse_range(const std::vector<std::pair<std::string,size_t> >& contigs,
                 const std::string& range_txt, range& ans);

// parse a comma-separated list of ranges
bool parse_ranges(const std::vector<std::pair<std::string,size_t> >& contigs,
                  const std::string& ranges, std::vector<range>& ans);

Status parse_bed_file(std::shared_ptr<spdlog::logger> logger,
                      const std::string &bedfilename,
                      const std::vector<std::pair<std::string,size_t> > &contigs,
                      std::vector<range> &ranges);

// Note: YAML serialization generates nice human readable
// files. However, it tends to be slow. An alternate implementation
// that uses cap'n proto (https://capnproto.org/index.html) has
// headers in file capnp_serialize.h. This is done only for the cases
// where efficiency matters.

// Read from disk and parse a YAML file
Status LoadYAMLFile(const std::string& filename, YAML::Node &node);


// Serialize the contigs to YAML
Status yaml_of_contigs(const std::vector<std::pair<std::string,size_t> > &contigs,
                       YAML::Emitter &yaml);
Status contigs_of_yaml(const YAML::Node& yaml,
                       std::vector<std::pair<std::string,size_t> > &contigs);

// Serialize N (sample count), contigs, and discovered alleles to YAML, in a
// streaming fashion.
Status yaml_stream_of_discovered_alleles(unsigned N, const std::vector<std::pair<std::string,size_t> > &contigs,
                                         const discovered_alleles &dsals,
                                         std::ostream &os);

// Load a YAML file, previously created with the above function
Status discovered_alleles_of_yaml_stream(std::istream &is,
                                         unsigned &N, std::vector<std::pair<std::string,size_t> > &contigs,
                                         discovered_alleles &dsals);

// Write the discovered alleles to a file
Status yaml_write_discovered_alleles_to_file(const discovered_alleles &dsals,
                                             const std::vector<std::pair<std::string,size_t>> &contigs,
                                             unsigned int sample_count,
                                             const std::string &filename);

// Serialize a vector of unified-alleles to YAML.
Status yaml_stream_of_unified_sites(const std::vector<unified_site> &sites,
                                    const std::vector<std::pair<std::string,size_t> > &contigs,
                                    std::ostream &os);

// Write the unified-sites to a file
Status write_unified_sites_to_file(const std::vector<unified_site> &sites,
                                   const std::vector<std::pair<std::string,size_t>> &contigs,
                                   const std::string &filename);

// Load from a file, data previously serialized with the above function
Status unified_sites_of_yaml_stream(std::istream &is,
                                    const std::vector<std::pair<std::string,size_t> > &contigs,
                                    std::vector<unified_site> &sites);

// Check if a file exists
bool check_file_exists(const std::string &path);

// Check if a directory exists
bool check_dir_exists(const std::string &path);

// Load configuration for the unifier and genotyper
// If name ends in ".yml" then the YAML file is parsed
// Otherwise name selects a hardcoded configuration preset.
Status load_config(std::shared_ptr<spdlog::logger> logger,
                   const std::string& name,
                   unifier_config& unifier_cfg,
                   genotyper_config& genotyper_cfg,
                   std::string& config_txt,
                   std::string& config_crc32c,
                   bool more_PL = false,
                   bool squeeze = false,
                   bool trim_uncalled_alleles = false);

std::string describe_config_presets();

RocksKeyValue::prefix_spec* GLnexus_prefix_spec();

// Initialize a database. Fills in the contigs.
Status db_init(std::shared_ptr<spdlog::logger> logger,
               const std::string &dbpath,
               const std::string &exemplar_gvcf,
               std::vector<std::pair<std::string,size_t>> &contigs, // output parameter
               size_t bucket_size = BCFKeyValueData::default_bucket_size);

// Read the contigs from a database
Status db_get_contigs(std::shared_ptr<spdlog::logger> logger,
                      const std::string &dbpath,
                      std::vector<std::pair<std::string,size_t> > &contigs);

// Load gvcf files into a database in parallel
Status db_bulk_load(std::shared_ptr<spdlog::logger> logger,
                    size_t mem_budget, size_t nr_threads,
                    const std::vector<std::string> &gvcfs,
                    const std::string &dbpath,
                    const std::vector<range> &ranges,   // limit the bulk load to these ranges
                    std::vector<std::pair<std::string,size_t>> &contigs, // output param
                    bool delete_gvcf_after_load = false);

// Discover alleles in the database. Return discovered alleles, and the sample count.
Status discover_alleles(std::shared_ptr<spdlog::logger> logger,
                        size_t mem_budget, size_t nr_threads,
                        const std::string &dbpath,
                        const std::vector<range> &ranges,
                        const std::vector<std::pair<std::string,size_t> > &contigs,
                        discovered_alleles &dsals,
                        unsigned &sample_count);

// Run unifier on given discovered alleles.
// input dsals is cleared by side-effect to save memory
// output sites is appended to (not cleared!)
Status unify_sites(std::shared_ptr<spdlog::logger> logger,
                   const unifier_config &unifier_cfg,
                   const std::vector<std::pair<std::string,size_t> > &contigs,
                   discovered_alleles &dsals,
                   unsigned sample_count,
                   std::vector<unified_site> &sites,
                   GLnexus::unifier_stats& stats);

// if the file name is "-", then output is written to stdout.
Status genotype(std::shared_ptr<spdlog::logger> logger,
                size_t mem_budget, size_t nr_threads,
                const std::string &dbpath,
                const GLnexus::genotyper_config &genotyper_cfg,
                const std::vector<unified_site> &sites,
                const std::vector<std::string> &extra_header_lines,
                const std::string &output_filename);

// compare different implementations of database iteration methods.
//
// n_iter: how many random queries to try
Status compare_db_itertion_algorithms(std::shared_ptr<spdlog::logger> logger,
                                      const std::string &dbpath,
                                      int n_iter);
}}}

#endif

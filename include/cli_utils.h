#ifndef GLNEXUS_APPLET_UTILS_H
#define GLNEXUS_APPLET_UTILS_H

#include <string>
#include <vector>
#include <map>
#include "types.h"
#include "spdlog/spdlog.h"
#include "RocksKeyValue.h"
#include "KeyValue.h"

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

// Merge a bunch of discovered-allele files, all in YAML format.
// N is summed across the files, and the contigs must be the same in all.
Status merge_discovered_allele_files(std::shared_ptr<spdlog::logger> logger,
                                     size_t nr_threads,
                                     const std::vector<std::string> &filenames,
                                     unsigned& N, std::vector<std::pair<std::string,size_t>> &contigs,
                                     discovered_alleles &dsals);


// Find which range contains [pos]. The ranges are assumed to be non-overlapping.
Status find_containing_range(const std::set<range> &ranges,
                             const range &pos,
                             range &ans);

// Check if a file exists
bool check_file_exists(const std::string &path);

// Check if a directory exists
bool check_dir_exists(const std::string &path);

// Recursively remove a directory path. Do not return an error
// if the path does not exist.
//
// Note: this could be done with boost::remove_all, but we do not
// want to require that dependency.
Status recursive_delete(const std::string &path);

// Load a named configuration for the unifier and genotyper
Status load_config_preset(std::shared_ptr<spdlog::logger> logger,
                          const std::string& name,
                          unifier_config& unifier_cfg,
                          genotyper_config& genotyper_cfg);

RocksKeyValue::prefix_spec* GLnexus_prefix_spec();

const int default_bucket_size = 30000;

// Initialize a database. Fills in the contigs.
Status db_init(std::shared_ptr<spdlog::logger> logger,
               const std::string &dbpath,
               const std::string &exemplar_gvcf,
               std::vector<std::pair<std::string,size_t>> &contigs, // output parameter
               size_t bucket_size = default_bucket_size);

// Read the contigs from a database
Status db_get_contigs(const std::string &dbpath,
                      std::vector<std::pair<std::string,size_t> > &contigs);

// Load gvcf files into a database in parallel
Status db_bulk_load(std::shared_ptr<spdlog::logger> logger,
                    const std::vector<std::string> &gvcfs,
                    const std::string &dbpath,
                    const std::vector<range> &ranges,   // limit the bulk load to these ranges
                    int nr_threads,
                    std::vector<std::pair<std::string,size_t>> &contigs, // output param
                    bool delete_gvcf_after_load = false);

// Discover alleles in the database. Return discovered alleles, and the sample count.
Status discover_alleles(std::shared_ptr<spdlog::logger> logger,
                        const std::string &dbpath,
                        const std::vector<range> &ranges,
                        const std::vector<std::pair<std::string,size_t> > &contigs,
                        int nr_threads,
                        discovered_alleles &dsals,
                        unsigned &sample_count);

Status unify_sites(std::shared_ptr<spdlog::logger> logger,
                   const unifier_config &unifier_cfg,
                   const std::vector<range> &ranges,
                   const std::vector<std::pair<std::string,size_t> > &contigs,
                   int nr_threads,
                   const discovered_alleles &dsals,
                   unsigned sample_count,
                   std::vector<unified_site> &sites);

// if the file name is "-", then output is written to stdout.
Status genotype(std::shared_ptr<spdlog::logger> logger,
                int nr_threads,
                const std::string &dbpath,
                const GLnexus::genotyper_config &genotyper_cfg,
                const std::vector<unified_site> &sites,
                const std::string &output_filename);

}}}

#endif

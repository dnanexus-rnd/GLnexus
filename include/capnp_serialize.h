#ifndef GLNEXUS_CAPNP_SERIALIZE_H
#define GLNEXUS_CAPNP_SERIALIZE_H

#include <string>
#include <vector>
#include <memory>
#include "types.h"

namespace GLnexus {
namespace capnp {

// Serialization of data structures with cap'n proto (https://capnproto.org/index.html)
//
// The issue this module tries to solve is that YAML serialization is slow for
// big data structures, and discovered-alleles tends to be large. Capnp serialization
// is very fast.

// write discovered_alleles structure to a file, with cap'n proto serialization
Status write_discovered_alleles(unsigned int sample_count,
                                const std::vector<std::pair<std::string,size_t> >& contigs,
                                const discovered_alleles &dsals,
                                const std::string &filename);

// write discovered_alleles structure to a file descriptor.
//
// Note: this will not work with a C++ stream, only a low level file descriptor.
Status write_discovered_alleles_fd(unsigned int sample_count,
                                   const std::vector<std::pair<std::string,size_t> >& contigs,
                                   const discovered_alleles &dsals,
                                   int fd);

// read discovered_alleles structure from a file, as serialized by cap'n proto
Status read_discovered_alleles(const std::string &filename,
                               unsigned int &sample_count,
                               std::vector<std::pair<std::string,size_t> >& contigs,
                               discovered_alleles &dsals);

// Verify that we can serialize and deserialize a discovered-alleles
// structure. Temporary results are written to [filename].
//
// Note: this is a debugging function
Status discover_alleles_verify(unsigned int sample_count,
                               const std::vector<std::pair<std::string,size_t> >& contigs,
                               const discovered_alleles &dsals,
                               const std::string &filename);

}}

#endif

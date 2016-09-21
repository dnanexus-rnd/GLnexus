#ifndef GLNEXUS_CAPNP_SERIALIZE_H
#define GLNEXUS_CAPNP_SERIALIZE_H

#include <string>
#include <vector>
#include <memory>
#include "types.h"

namespace GLnexus {
namespace capnp {

// Serialization of data structures with cap'n proto (https://capnproto.org/index.html)

// write discovered_alleles structure to a file-descriptor, with cap'n proto serialization
Status write_discovered_alleles(int fd,
                                const std::vector<std::pair<std::string,size_t> >& contigs,
                                const discovered_alleles &dsals);

// read discovered_alleles structure from a file-descriptor, as serialized by cap'n proto
Status read_discovered_alleles(int fd,
                               const std::vector<std::pair<std::string,size_t> >& contigs,
                               discovered_alleles &dsals);

// Verify that we can serialize and deserialize a discovered-alleles
// structure. Temporary results are written to [filename].
Status discover_alleles_verify(const std::vector<std::pair<std::string,size_t> >& contigs,
                               const discovered_alleles &dsals,
                               const std::string &filename);

}}

#endif

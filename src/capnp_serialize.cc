#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <capnp/message.h>
#include <capnp/serialize-packed.h>
#include <iostream>

#include <defs.capnp.h>
#include <capnp_serialize.h>

using namespace std;

namespace GLnexus {
namespace capnp {

// write discovered_alleles structure to a file-descriptor, with cap'n proto serialization
Status write_discovered_alleles(int fd,
                                const discovered_alleles&,
                                const std::vector<std::pair<std::string,size_t> >& contigs) {
    return Status::Invalid("Not implemented yet");
}

// read discovered_alleles structure from a file-descriptor, as serialized by cap'n proto
Status read_discovered_alleles(int fd,
                               const std::vector<std::pair<std::string,size_t> >& contigs,
                               discovered_alleles&) {
    return Status::Invalid("Not implemented yet");
}

Status discover_alleles_verify(const discovered_alleles &dal) {
    return Status::Invalid("Not implemented yet");
}

}} // namespace GLnexus

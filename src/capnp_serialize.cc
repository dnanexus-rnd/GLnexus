#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
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
                                const std::vector<std::pair<std::string,size_t> >& contigs,
                                const discovered_alleles &dsals) {
    return Status::Invalid("Not implemented yet");
}

// read discovered_alleles structure from a file-descriptor, as serialized by cap'n proto
Status read_discovered_alleles(int fd,
                               const std::vector<std::pair<std::string,size_t> >& contigs,
                               discovered_alleles&) {
    return Status::Invalid("Not implemented yet");
}

Status discover_alleles_verify(const std::vector<std::pair<std::string,size_t> >& contigs,
                               const discovered_alleles &dsals) {
    Status s;

    // create temporary file
    std::FILE* tmpf = std::tmpfile();
    int fd = fileno(tmpf);

    // write to it
    S(write_discovered_alleles(fd, contigs, dsals));

    // read from it
    discovered_alleles dsals2;
    S(read_discovered_alleles(fd, contigs, dsals2));

    // Close the file-descriptor, this will also delete the file
    // once the program is closed.
    int rc = close(fd);
    if (rc != 0)
        return Status::IOError("Could not close file-descriptor to temporary file");

    // verify we get the same alleles back
    if (dsals == dsals2) {
        return Status::OK();
    }
    return Status::Invalid("capnp serialization/deserialization of discovered alleles does not return original value");
}

}} // namespace GLnexus

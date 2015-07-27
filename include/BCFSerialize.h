#ifndef GLNEXUS_BCF_SERIALIZE_H
#define GLNEXUS_BCF_SERIALIZE_H
#include "vcf.h"
#include "types.h"

namespace GLnexus {

// Classes for serializing BCF structures to/from RAM

// The following three functions, prefixed with bcf_raw, are copied
// and modified from the htslib sources. They are used to read/write
// uncompressed BCF records from/to memory. They are declared for
// testing purposes, they are not to be used without the proper wrappers
// provided by the BCFWriter and BCFReader classes.

// Ccalculate the amount of bytes it would take to pack this bcf1 record.
int bcf_raw_calc_packed_len(bcf1_t *v);

//  Write the BCF record directly to a memory location. Return how much
//  space was used.
void bcf_raw_write_to_mem(bcf1_t *v, int reclen, char *addr);

//  Read a BCF record from memory, return the length of the packed record in RAM.
int bcf_raw_read_from_mem(const char *addr, bcf1_t *v);

class BCFWriter {
 private:
    static int STACK_ALLOC_LIMIT;

    std::ostringstream oss_;
    int valid_bytes_ = 0;

 public:
    BCFWriter();

    /// Serializes BCF records into a memory buffer (without the header)
    static Status Open(std::unique_ptr<BCFWriter>& ans);
    ~BCFWriter();
    Status write(bcf1_t* x);
    Status contents(std::string& ans);

    /// Writes a BCF header into a memory buffer
    static std::string write_header(const bcf_hdr_t *hdr);
};

class BCFReader {
 private:
    BCFReader(const char* buf, size_t bufsz);

 public:
    /// Reads BCF records from a memory buffer (without the header)
    static Status Open(const char* buf, size_t bufsz, std::unique_ptr<BCFReader>& ans);
    ~BCFReader();
    Status read(std::shared_ptr<bcf1_t>& ans);

    /// Reads a BCF header from a memory buffer
    static Status read_header(const char* buf, int len, int& consumed, std::shared_ptr<bcf_hdr_t>& ans);
};

} // namespace GLnexus

#endif

#ifndef GLNEXUS_BCF_SERIALIZE_H
#define GLNEXUS_BCF_SERIALIZE_H
#include "vcf.h"
#include "types.h"

namespace GLnexus {
// Classes for serializing BCF structures to/from RAM

class BCFWriter {
 public:
    static int INIT_SIZE;
    static int SIZE_MULTIPLIER;
    static int MAX_BUF_SIZE;
    static int STACK_ALLOC_LIMIT;

    std::ostringstream oss_;
    int valid_bytes_ = 0;
    char *scratch_pad_ = NULL;
    int scratch_pad_size_ = 0;

    BCFWriter();

    static Status Open(std::unique_ptr<BCFWriter>& ans);
    ~BCFWriter();
    Status write(bcf1_t* x);
    Status contents(std::string& ans);

    static std::string write_header(const bcf_hdr_t *hdr);
    };

class BCFReader {
 public:
    const char* buf_ = nullptr;
    size_t bufsz_;
    int current_ = 0;

    BCFReader(const char* buf, size_t bufsz);

    static Status Open(const char* buf, size_t bufsz,  std::unique_ptr<BCFReader>& ans);
    ~BCFReader();
    Status read(std::shared_ptr<bcf1_t>& ans);

    static bcf_hdr_t* read_header(const char* buf, int len);
};

} // namespace GLnexus

#endif

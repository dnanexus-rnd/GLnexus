#ifndef GLNEXUS_BCF_SERIALIZE_H
#define GLNEXUS_BCF_SERIALIZE_H
#include "vcf.h"
#include "types.h"

namespace GLnexus {
// Classes for serializing BCF structures to/from RAM

// Read BCF header
Status bcf_raw_read_header(const uint8_t* buf, int hdrlen,
                           int& consumed, std::shared_ptr<bcf_hdr_t>& ans);

// Write BCF header
std::string bcf_write_header(const bcf_hdr_t *hdr);

// The following three functions, prefixed with bcf_raw, are copied
// and modified from the htslib sources. They are used to read/write
// uncompressed BCF records from/to memory. They are declared for
// testing purposes, they are not to be used without the proper wrappers
// provided by the BCFWriter and BCFScanner classes.

// Calculate the amount of bytes it would take to pack this bcf1 record.
int bcf_raw_calc_packed_len(bcf1_t *v);

// Read a BCF record from memory, return the length of the packed record in [ans].
// Return Invalid in case of error (out of bounds memory access).
//
// Returns:
// 1) A BCF record
// 2) length of memory extent read
// 3) Status code
Status bcf_raw_read_from_mem(const uint8_t *buf, int start, size_t len, bcf1_t *v,
                             int &ans);

// Quickly read the range from the BCF record (without deserializing it entirely)
Status bcf_raw_range(const uint8_t *buf, int start, size_t len, range& rng);

// convert a BCF record into a VCF string (no newline)
std::shared_ptr<std::string> bcf1_to_string(const bcf_hdr_t *hdr, const bcf1_t *bcf);

/*
  Write the BCF record directly to a memory location. This function is
  unsafe, it assumes that the caller allocated sufficient space, by
  calling [bcf_raw_calc_packed_len] beforehand.
*/
void bcf_raw_write_to_mem(bcf1_t *v, int reclen, uint8_t *addr);

// Return 1 if the records are the same, 0 otherwise.
// This compares most, but not all, fields.
int bcf_shallow_compare(const bcf1_t *x, const bcf1_t *y);

} // namespace GLnexus

#endif

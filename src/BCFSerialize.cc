/*
 Helper classes to de/serialize BCF records to memory buffers The code
 that reads and writes bcf1_t records/headers is adapted from the
 htslib library. We include htslib license below.


vcf.c -- VCF/BCF API functions.

    Copyright (C) 2012, 2013 Broad Institute.
    Copyright (C) 2012-2014 Genome Research Ltd.
    Portions copyright (C) 2014 Intel Corporation.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/

#include <assert.h>
#include <alloca.h>
#include <iostream>
#include "BCFSerialize.h"
using namespace std;

namespace GLnexus {

// Ccalculate the amount of bytes it would take to pack this bcf1 record.
int bcf_raw_calc_packed_len(bcf1_t *v)
{
    return 32 + v->shared.l + v->indiv.l;
}

/*
  Write the BCF record directly to a memory location. Return how much
  space was used.
   Note: the code is adapted from the bcf_write routine in htslib/vcf.c.
  The original prototype is:
      int bcf_write(htsFile *hfp, const bcf_hdr_t *h, bcf1_t *v)
  The original routine writes to a file, not to memory.
*/
void bcf_raw_write_to_mem(bcf1_t *v, int reclen, char *addr) {
    int loc = 0;
    uint32_t x[8];
    assert(sizeof(x) == 32);
    x[0] = v->shared.l + 24; // to include six 32-bit integers
    x[1] = v->indiv.l;
    memcpy(x + 2, v, 16);
    x[6] = (uint32_t)v->n_allele<<16 | v->n_info;
    x[7] = (uint32_t)v->n_fmt<<24 | v->n_sample;
     memcpy(&addr[loc], (char*)x, sizeof(x));
    loc += sizeof(x);
    memcpy(&addr[loc], v->shared.s, v->shared.l);
    loc += v->shared.l;
    memcpy(&addr[loc], v->indiv.s, v->indiv.l);
    loc += v->indiv.l;
     assert(loc == reclen);
}

/*
    Read a BCF record from memory, return the length of the packed record in RAM.

    Note: the code is adapted from the bcf_read1_core routine in htslib/vcf.c.
    The original prototype is:
         int bcf_read1_core(BGZF *fp, bcf1_t *v)
    The original routine reads from a file, not from memory.
*/
int bcf_raw_read_from_mem(const char *addr, bcf1_t *v) {
    int loc = 0;
    uint32_t x[8];
    memcpy(x, &addr[loc], 32);
    loc += 32;

    assert(x[0] > 0);
    x[0] -= 24; // to exclude six 32-bit integers
    ks_resize(&v->shared, x[0]);
    ks_resize(&v->indiv, x[1]);
    memcpy(v, (char*)&x[2], 16);
    v->n_allele = x[6]>>16; v->n_info = x[6]&0xffff;
    v->n_fmt = x[7]>>24; v->n_sample = x[7]&0xffffff;
    v->shared.l = x[0], v->indiv.l = x[1];

    // silent fix of broken BCFs produced by earlier versions of
    // bcf_subset, prior to and including bd6ed8b4
    if ( (!v->indiv.l || !v->n_sample) && v->n_fmt ) v->n_fmt = 0;

    memcpy(v->shared.s, &addr[loc], v->shared.l);
    loc += v->shared.l;
    memcpy(v->indiv.s, &addr[loc], v->indiv.l);
    loc += v->indiv.l;

    return loc;
}


int BCFWriter::STACK_ALLOC_LIMIT = 32 * 1024;

BCFWriter::BCFWriter() {}

Status BCFWriter::Open(unique_ptr<BCFWriter>& ans) {
    ans.reset(new BCFWriter);
    ans->valid_bytes_ = 0;
    //ans->oss_;
    return Status::OK();
}

BCFWriter::~BCFWriter() {
    oss_.clear();
    valid_bytes_ = 0;
}

Status BCFWriter::write(bcf1_t* x) {
    int reclen = bcf_raw_calc_packed_len(x);

    // Note: allocation on the stack for small memory
    // sizes, this should be the normal usage case. The idea
    // is to avoid contention if multiple threads access this
    // method.
    char *scratch_pad;
    bool heap_allocation = false;
    if (reclen <= STACK_ALLOC_LIMIT) {
        scratch_pad = (char*) alloca(reclen);
    } else {
        heap_allocation = true;
        scratch_pad = (char*) malloc(reclen);
    }

    // Separate the C code, from the C++ code
    // Serialize the bcf1_t stuct into a [char*], then
    // append it to the end of the buffer.
    bcf_raw_write_to_mem(x, reclen, scratch_pad);
    oss_.write(scratch_pad, reclen);
    valid_bytes_ += reclen;

    if (heap_allocation) {
        free(scratch_pad);
    }
    return Status::OK();
}

Status BCFWriter::contents(string& ans) {
    ans.clear();
    ans.append(oss_.str(), 0, valid_bytes_);
    return Status::OK();
}

/* Adapted from [htslib::vcf.c::bcf_hdr_write] to write
   to memory instead of disk.
*/
std::string BCFWriter::write_header(const bcf_hdr_t *hdr) {
    int hlen;
    char *htxt = bcf_hdr_fmt_text(hdr, 1, &hlen);
    hlen++; // include the \0 byte

    char *buf = (char*) malloc(5 + 4 + hlen);
    assert(buf != NULL);
    int loc = 0;
    memcpy(&buf[loc], "BCF\2\2", 5);
    loc += 5;
    memcpy(&buf[loc], &hlen, 4);
    loc += 4;
    memcpy(&buf[loc], htxt, hlen);
    loc += hlen;

    // cleanup and return
    string rc(buf, loc);
    free(htxt);
    free(buf);
    return rc;
}



// constructor
BCFReader::BCFReader(const char* buf, size_t bufsz) :
    buf_(buf), bufsz_(bufsz)
{}

Status BCFReader::Open(const char* buf,
                       size_t bufsz,
                       unique_ptr<BCFReader>& ans) {
    ans.reset(new BCFReader(buf, bufsz));
    return Status::OK();
}

BCFReader::~BCFReader() {
    buf_ = NULL;
    bufsz_ = 0;
    current_ = 0;
}

Status BCFReader::read(shared_ptr<bcf1_t>& ans) {
    if (!ans) {
        ans = shared_ptr<bcf1_t>(bcf_init(), &bcf_destroy);
    }
    if ((size_t)current_ >= bufsz_)
        return Status::NotFound();
    int reclen = bcf_raw_read_from_mem(&buf_[current_], ans.get());
    current_ += reclen;
    return Status::OK();
}

/* Adapted from [htslib::vcf.c::bcf_hdr_read] to read
   from memory instead of disk.
*/
bcf_hdr_t* BCFReader::read_header(const char* buf, int hdrlen) {
    if (strncmp(buf, "BCF\2\2", 5) != 0) {
        //fprintf("invalid BCF2 magic string");
        return NULL;
    }

    bcf_hdr_t *hdr = bcf_hdr_init("r");
    int loc = 5;
    int hlen;
    char *htxt;
    memcpy(&hlen, &buf[loc], 4);
    loc += 4;

    // sanity check; make sure we do not  over the record length
    assert(loc + hlen <= hdrlen);

    htxt = (char*)malloc(hlen);
    memcpy(htxt, &buf[loc], hlen);
    bcf_hdr_parse(hdr, htxt);

    // release resources
    free(htxt);
    return hdr;
}

} // namespace GLnexus

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

static void range_check(int x, size_t y) {
    // FIXME: this should -not- be removed in Release builds
    assert(x <= y);
}

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
// FIXME: add range checks
int bcf_raw_read_from_mem(const char *buf, int start, size_t len, bcf1_t *v) {
    int loc = 0;
    uint32_t x[8];
    range_check(start + 32, len);
    memcpy(x, &buf[start + loc], 32);
    loc += 32;

    assert(x[0] > 0);
    x[0] -= 24; // to exclude six 32-bit integers
    range_check(start + loc + x[0] + x[1], len);
    ks_resize(&v->shared, x[0]);
    ks_resize(&v->indiv, x[1]);
    memcpy(v, (char*)&x[2], 16);
    v->n_allele = x[6]>>16; v->n_info = x[6]&0xffff;
    v->n_fmt = x[7]>>24; v->n_sample = x[7]&0xffffff;
    v->shared.l = x[0], v->indiv.l = x[1];

    // silent fix of broken BCFs produced by earlier versions of
    // bcf_subset, prior to and including bd6ed8b4
    if ( (!v->indiv.l || !v->n_sample) && v->n_fmt ) v->n_fmt = 0;

    memcpy(v->shared.s, &buf[start + loc], v->shared.l);
    loc += v->shared.l;
    memcpy(v->indiv.s, &buf[start + loc], v->indiv.l);
    loc += v->indiv.l;

    return loc;
}

// Calculate the total length of a record, without unpacking it.
uint32_t bcf_raw_calc_rec_len(const char *buf, int start, size_t len) {
    range_check(start + 32, len);
    uint32_t *x = (uint32_t*) &buf[start];

    assert(x[0] > 0);
    uint32_t rc = 32 + (x[0] -24) + x[1];
    return rc;
}

// Check if the record that starts at memory address [addr] overlaps
// the range [rng]. The trick is to do this without unpacking the record.
bool bcf_raw_overlap(const char *buf, int start, size_t len, const range &rng) {
    range_check(start + 32, len);
    uint32_t *x = (uint32_t*) &buf[start];

    int32_t rid = x[2];
    int32_t beg = x[3];
    int32_t rlen = x[4];
    int32_t end = beg + rlen;

    return rid == rng.rid && end > rng.beg && beg < rng.end;
}

// Return 1 if the records are the same, 0 otherwise.
// This compares most, but not all, fields.
int bcf_shallow_compare(const bcf1_t *x, const bcf1_t *y) {
    if (x->rid != y->rid) return 0;
    if (x->pos != y->pos) return 0;
    if (x->rlen != y->rlen) return 0;

    // I think nan != nan, so we are skipping this comparison
    //if (x->qual != y->qual) return 0;

    if (x->n_info != y->n_info) return 0;
    if (x->n_allele != y->n_allele) return 0;
    if (x->n_sample != y->n_sample) return 0;
    if (x->shared.l != y->shared.l) return 0;
    if (x->indiv.l != y->indiv.l) return 0;

    return 1;
}

/* Adapted from [htslib::vcf.c::bcf_hdr_read] to read
   from memory instead of disk.
*/
static Status bcf_raw_read_header(const char* buf,
                                  int hdrlen,
                                  int& consumed,
                                  shared_ptr<bcf_hdr_t>& ans) {
    if (strncmp(buf, "BCF\2\2", 5) != 0) {
        return Status::Invalid("BCFSerialize::bcf_raw_read_header invalid BCF2 magic string");
    }

    ans = shared_ptr<bcf_hdr_t>(bcf_hdr_init("r"), &bcf_hdr_destroy);
    int loc = 5;
    int hlen;
    memcpy(&hlen, &buf[loc], 4);
    loc += 4;

    // make sure we do not read beyond the buffer end
    if (loc + hlen > hdrlen) {
        return Status::Invalid("BCFSerialize::bcf_raw_read_header truncated header");
    }

    auto htxt = make_unique<char[]>(hlen);
    memcpy(htxt.get(), &buf[loc], hlen);

    // FIXME bcf_hdr_parse does not seem to check bounds -- potential security issue
    if (bcf_hdr_parse(ans.get(), htxt.get()) != 0) {
        return Status::Invalid("BCFSerialize::bcf_raw_read_header parse error");
    }

    consumed = loc + hlen;
    return Status::OK();
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
    num_entries_ = 0;
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
    num_entries_ ++;

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

int BCFWriter::get_num_entries() const {
    return num_entries_;
}

/* Adapted from [htslib::vcf.c::bcf_hdr_write] to write
   to memory instead of disk.
*/
std::string BCFWriter::write_header(const bcf_hdr_t *hdr) {
    int hlen;
    char *htxt = bcf_hdr_fmt_text(hdr, 1, &hlen);
    hlen++; // include the \0 byte

    char *buf = (char*) malloc(5 + 4 + hlen);
    assert(buf != nullptr);
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
    buf_ = nullptr;
    bufsz_ = 0;
    current_ = 0;
}

Status BCFReader::read(shared_ptr<bcf1_t>& ans) {
    if (!ans) {
        ans = shared_ptr<bcf1_t>(bcf_init(), &bcf_destroy);
    }
    if ((size_t)current_ >= bufsz_)
        return Status::NotFound();
    int reclen = bcf_raw_read_from_mem(buf_, current_, bufsz_, ans.get());
    current_ += reclen;
    return Status::OK();
}

Status BCFReader::read_header(const char* buf, int hdrlen, int& consumed, shared_ptr<bcf_hdr_t>& ans) {
    return bcf_raw_read_header(buf, hdrlen, consumed, ans);
}

// BCFScanner
// constructor
BCFScanner::BCFScanner(const char* buf, size_t bufsz) :
    buf_(buf), bufsz_(bufsz)
{}

Status BCFScanner::Open(const char* buf,
                       size_t bufsz,
                       unique_ptr<BCFScanner>& ans) {
    ans.reset(new BCFScanner(buf, bufsz));
    return Status::OK();
}

BCFScanner::~BCFScanner() {
    buf_ = nullptr;
    bufsz_ = 0;
    current_ = 0;
}

 // move the cursor to the next record
Status BCFScanner::next() {
    if ((size_t)current_ >= bufsz_)
        return Status::NotFound();
    current_ += bcf_raw_calc_rec_len(buf_, current_, bufsz_);
    if ((size_t)current_ >= bufsz_)
        return Status::NotFound();
    return Status::OK();
}

Status BCFScanner::read(shared_ptr<bcf1_t>& ans) {
    if (!ans) {
        ans = shared_ptr<bcf1_t>(bcf_init(), &bcf_destroy);
    }
    bcf_raw_read_from_mem(buf_, current_, bufsz_, ans.get());
    return Status::OK();
}

// check if the current record overlaps a range
bool BCFScanner::overlaps(const range &rng) {
    return bcf_raw_overlap(buf_, current_, bufsz_, rng);
}

/* Adapted from [htslib::vcf.c::bcf_hdr_read] to read
   from memory instead of disk.
*/
Status BCFScanner::read_header(const char* buf, int hdrlen, int& consumed, shared_ptr<bcf_hdr_t>& ans) {
    return bcf_raw_read_header(buf, hdrlen, consumed, ans);
}

} // namespace GLnexus

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

// sugar for memory bounds-checking during BCF deserialization
#define BOUNDS_CHECK(desired, limit, errmsg) \
    {if ((desired) > (limit)) { \
        string msg = string(errmsg) + ". Reading offset "; \
        msg += to_string(desired) + " in buffer of length "; \
        msg += std::to_string(limit); \
        return Status::Invalid("Failed memory bounds check", msg.c_str()); \
    }}

// Calculate the amount of bytes it would take to pack this bcf1 record.
int bcf_raw_calc_packed_len(bcf1_t *v)
{
    return 32 + v->shared.l + v->indiv.l;
}

/*
  Write the BCF record directly to a memory location. This function is
  unsafe, it assumes that the caller allocated sufficient space, by
  calling [bcf_raw_calc_packed_len] beforehand.

   Note: the code is adapted from the bcf_write routine in htslib/vcf.c.
  The original prototype is:
      int bcf_write(htsFile *hfp, const bcf_hdr_t *h, bcf1_t *v)
  The original routine writes to a file, not to memory.
*/
void bcf_raw_write_to_mem(bcf1_t *v, int reclen, uint8_t *addr) {
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
    Read a BCF record from memory, returns:
1) A BCF record
2) length of memory extent read
3) Status code

    Note: the code is adapted from the bcf_read1_core routine in htslib/vcf.c.
    The original prototype is:
         int bcf_read1_core(BGZF *fp, bcf1_t *v)
    The original routine reads from a file, not from memory.

    If v is being reused (rather than a newly allocated bcf1_t) then the
    caller should first sanitize it with bcf_clear1 as in bcf_read1_core.
*/
Status bcf_raw_read_from_mem(const uint8_t *buf, int start, size_t len, bcf1_t *v,
                             int &ans) {
    int loc = 0;
    uint32_t x[8];

    BOUNDS_CHECK(start + 32, len, "reading header of BCF record");
    memcpy(x, &buf[start + loc], 32);
    loc += 32;
    assert(x[0] > 0);
    x[0] -= 24; // to exclude six 32-bit integers

    BOUNDS_CHECK(start + loc + x[0] + x[1], len, "reading BCF record");
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

    ans = loc;
    return Status::OK();
}

// Calculate the total length of a record, without unpacking it.
Status bcf_raw_calc_rec_len(const uint8_t *buf, int start, size_t len, uint32_t &ans) {
    BOUNDS_CHECK(start + 32, len, "calculating BCF record length");

    uint32_t *x = (uint32_t*) &buf[start];
    assert(x[0] > 0);
    ans = 32 + (x[0] -24) + x[1];
    return Status::OK();
}

// Get the range of the serialized BCF record starting at memory address
// [buf+start], without deserializing the whole record.
Status bcf_raw_range(const uint8_t *buf, int start, size_t len, range& rng) {
    BOUNDS_CHECK(start + 32, len, "reading BCF record range");

    uint32_t *x = (uint32_t*) &buf[start];
    uint32_t rid = x[2];
    uint32_t beg = x[3];
    uint32_t rlen = x[4];
    rng = range(rid, beg, beg + rlen);
    return Status::OK();
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
Status bcf_raw_read_header(const uint8_t* buf, int hdrlen,
                           int& consumed, shared_ptr<bcf_hdr_t>& ans) {
    if (hdrlen < 9 || strncmp((const char*)buf, "BCF\2\2", 5) != 0) {
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

std::string bcf_write_header(const bcf_hdr_t *hdr) {
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

// convert a BCF record into a string
shared_ptr<string> bcf1_to_string(const bcf_hdr_t *hdr, const bcf1_t *bcf) {
    kstring_t kstr;
    memset((void*)&kstr, 0, sizeof(kstring_t));
    vcf_format(hdr, bcf, &kstr);
    auto retval = make_shared<string>(kstr.s, kstr.l);
    if (retval->length() > 0)
        retval->erase(retval->find_last_not_of(" \n\r\t")+1); // get rid of extra spaces at the end

    // cleanup
    if (kstr.s != NULL)
        free(kstr.s);

    return retval;
}

} // namespace GLnexus

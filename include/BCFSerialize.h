#ifndef GLNEXUS_BCF_SERIALIZE_H
#define GLNEXUS_BCF_SERIALIZE_H
#include "vcf.h"

namespace GLnexus {

class BCFWriter {
 private:
    // A string that that keeps serialized records. We add to its end addtional records
    // as these are written.
    string buf_;
    int valid_bytes_ = 0;
    char *scratch_pad_ = NULL;
    int scratch_pad_size_ = 0;
    
public:

    static Status Open(unique_ptr<BCFWriter>& ans) {
        ans.reset(new BCFWriter);
        ans.valid_bytes_ = 0;
        ans.buf_.resize(INIT_SIZE);
        ans.scratch_pad_ = (char*) malloc(INIT_SIZE);
        assert(ans.scratch_pad_ != NULL);
        ans.scratch_pad_size_ = INIT_SIZE;
        
        return Status::OK();
    }

    virtual ~BCFWriter() {
        buf_.resize(0);
        valid_bytes_ = 0;
        free(scratch_pad_);
        scratch_pad_ = NULL;
        scratch_pad_size_ = 0;
    }

    Status write(bcf1_t* x) {
        int reclen = bcf_calc_packed_len(x);

        // Make sure the buffer has enough space for the
        // new record. Multiply the size by a constant factor
        // until we have sufficient space.
        int remaining_len = buf_.capacity() - valid_bytes_;
        while (remaining_len < reclen) {
            buf_.resize(buf_.capacity() * SIZE_MULTIPLIER);
            remaining_len = buf_.capacity() - valid_bytes_;
            assert(buf_.capacity() < MAX_REASONABLE_SIZE);
        }

        // make sure we have enough scratch space
        if (scratch_pad_size_ < reclen) {
            free(scratch_pad_);
            scratch_pad_ = (char*) malloc(reclen);
            assert(scratch_pad_ != NULL);
            scratch_pad_size_ = reclen;
        }

        // Separate the C code, from the C++ code
        // Serialize the bcf1_t stuct into a [char*], then
        // append it to the end of the buffer.
        bcf_write_to_mem(x, reclen, scratch_pad_);
        buf_.append(scratch_pad_, reclen);
        return Status::OK();
    }

    Status contents(string& ans) {
        ans.clear();
        ans.replace(0, valid_bytes_, buf_);
        return Status::OK();
    }
};

class BCFReader {
    const char* buf_ = nullptr;
    size_t bufsz_;
    int current_ = 0;
    
    BCFReader(const char* buf, size_t bufsz) : buf_(buf), bufsz_(bufsz) {}

    /*
      Read a BCF record from memory, return the length of the packed record in RAM.

      Note: the code is adapted from the bcf_read1_core routine in htslib/vcf.c.
      The original prototype is:
          int bcf_read1_core(BGZF *fp, bcf1_t *v)
      The original routine reads from a file, not from memory.
    */
    static int bcf_read_from_mem(char *addr, bcf1_t *v) {
        int loc = 0;
        uint32_t x[8];
        memcpy(x, &addr[loc], 32);
        loc += 32;

        x[0] -= 24; // to exclude six 32-bit integers
        ks_resize(&v->shared, x[0]);
        ks_resize(&v->indiv, x[1]);
        memcpy(v, x + 2, 16);
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

public:
    static Status Open(const char* buf, size_t bufsz, unique_ptr<BCFReader>& ans) {
        ans.reset(new BCFReader(buf, bufsz));
        return Status::OK();
    }

    virtual ~BCFReader() {
        buf_ = NULL;
        bufsz_ = 0;
        current_ = 0;
    }

    Status read(shared_ptr<bcf1_t>& ans) {
        if (!ans) {
            ans = shared_ptr<bcf1_t>(bcf_init(), &bcf_destroy);
        }

        if (current_ >= bufsz)
            return Status::Not_Found;
        int reclen = bcf_read_from_mem(&buf_[current_], ans.get());
        current_ += reclen;
        return Status::OK();
    }
};

} // namespace GLnexus

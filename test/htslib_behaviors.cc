// Test certain non-obvious behaviors of htslib 

#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <memory>

#include <vcf.h>
#include <hfile.h>

#include "BCFSerialize.h"
#include "catch.hpp"

using namespace std;

// sugar for declaring a unique_ptr with a custom deleter function
#define UPD(T,name,ini,del) std::unique_ptr<T, void(*)(T*)> up_##name((ini), (del)); auto name = up_##name.get();

TEST_CASE("htslib VCF missing data representation") {
    // Verify how htslib data structures represent missing data such as ./.
    // genotypes, . genotype likelihoods, etc.
    // https://github.com/dnanexus/htslib/blob/dnanexus/htslib/vcf.h
    UPD(vcfFile, vcf, bcf_open("test/data/missingdata.vcf", "r"), [](vcfFile* f) { bcf_close(f); });
    UPD(bcf_hdr_t, hdr, bcf_hdr_read(vcf), &bcf_hdr_destroy);
    UPD(bcf1_t, vt, bcf_init(), &bcf_destroy);

    REQUIRE(bcf_hdr_nsamples(hdr) == 3);

    int *gt = nullptr, gtsz = 0, glsz = 0, dpsz = 0;
    float *gl = nullptr;
    int32_t *dp = nullptr;

    REQUIRE(bcf_read(vcf, hdr, vt) == 0);
    REQUIRE(vt->pos == 105);

    REQUIRE(bcf_read(vcf, hdr, vt) == 0);
    REQUIRE(vt->pos == 106);
    REQUIRE(bcf_unpack(vt, BCF_UN_ALL) == 0);

    // alleles
    REQUIRE(strcmp(vt->d.allele[0], "CA") == 0);
    REQUIRE(strcmp(vt->d.allele[1], "GA") == 0);
    REQUIRE(strcmp(vt->d.allele[2], "GT") == 0);

    // genotypes
    REQUIRE(bcf_get_genotypes(hdr, vt, &gt, &gtsz) == 6);
    REQUIRE(bcf_gt_allele(gt[0]) == 0);  // 0/0
    REQUIRE(bcf_gt_allele(gt[1]) == 0);
    REQUIRE(bcf_gt_allele(gt[2]) == 0);  // 0/2
    REQUIRE(bcf_gt_allele(gt[3]) == 2);
    REQUIRE(bcf_gt_allele(gt[4]) == -1); // ./.
    REQUIRE(bcf_gt_allele(gt[5]) == -1);

    // genotype likelihoods
    // 0,-1.80618,-24.7,.,.,.
    REQUIRE(bcf_get_format_float(hdr, vt, "GL", &gl, &glsz) == 18);
    REQUIRE(isfinite(gl[0]));
    REQUIRE(isfinite(gl[1]));
    REQUIRE(isfinite(gl[2]));
    REQUIRE(isnan(gl[3]));
    REQUIRE(isnan(gl[4]));
    REQUIRE(isnan(gl[5]));

    REQUIRE(bcf_read(vcf, hdr, vt) == 0);
    REQUIRE(vt->pos == 107);
    REQUIRE(bcf_unpack(vt, BCF_UN_ALL) == 0);

    // more genotypes
    REQUIRE(bcf_get_genotypes(hdr, vt, &gt, &gtsz) == 6);
    REQUIRE(bcf_gt_is_missing(gt[0]));    // .
    REQUIRE(gt[1] == bcf_int32_vector_end);
    REQUIRE(bcf_gt_is_missing(gt[2]));    // ./.
    REQUIRE(bcf_gt_is_missing(gt[3]));
    REQUIRE_FALSE(bcf_gt_is_missing(gt[4]));
    REQUIRE_FALSE(bcf_gt_is_missing(gt[5]));
    REQUIRE(bcf_gt_allele(gt[4]) == 0);   // 0/1
    REQUIRE(bcf_gt_allele(gt[5]) == 1);

    // genotype likelihoods '.'
    REQUIRE(bcf_get_format_float(hdr, vt, "GL", &gl, &glsz) == 9);
    for (int i = 0; i < 6; i++) {
        REQUIRE(isnan(gl[i]));
    }

    // DP :.: :.: :6:
    REQUIRE(bcf_get_format_int32(hdr, vt, "DP", &dp, &dpsz) == 3);
    REQUIRE(dp[0] == bcf_int32_missing);
    REQUIRE(dp[1] == bcf_int32_missing);
    REQUIRE(dp[2] == 6);

    if (gt) free(gt);
    if (gl) free(gl);
    if (dp) free(dp);
}

TEST_CASE("htslib VCF header chrom injection") {
    // verify our method to inject contig entries into a vcf_hdr_t to preclude
    // BCF_ERR_CTG_UNDEF errors and stderr warnings like
    // "[W::vcf_parse] contig '21' is not defined in the header."
    // when a VCF file's header doesn't enumerate the contigs.
    UPD(vcfFile, vcf, bcf_open("test/data/trio_denovo.vcf", "r"), [](vcfFile* f) { bcf_close(f); });
    UPD(bcf_hdr_t, hdr, bcf_hdr_read(vcf), &bcf_hdr_destroy);
    UPD(bcf1_t, vt, bcf_init(), &bcf_destroy);

    // read a record with no injection, we should see BCF_ERR_CTG_UNDEF.
    REQUIRE(bcf_read(vcf, hdr, vt) == 0);
    REQUIRE(vt->pos == 4760);
    REQUIRE((vt->errcode & BCF_ERR_CTG_UNDEF) != 0);

    // inject
    bcf_hdr_append(hdr,"##contig=<ID=21,length=1000000>");
    bcf_hdr_sync(hdr);

    // read another record - we shouldn't get the error code
    vt->errcode = 0;
    REQUIRE(bcf_read(vcf, hdr, vt) == 0);
    REQUIRE(vt->pos == 70549);
    REQUIRE(vt->errcode == 0);
}

TEST_CASE("htslib VCF header synthesis") {
    shared_ptr<bcf_hdr_t> hdr(bcf_hdr_init("w"), &bcf_hdr_destroy);
    
    REQUIRE(bcf_hdr_append(hdr.get(),"##contig=<ID=A,length=1000000>") == 0);
    REQUIRE(bcf_hdr_append(hdr.get(),"##contig=<ID=B,length=100000>") == 0);
    REQUIRE(bcf_hdr_append(hdr.get(),"##contig=<ID=C,length=10000>") == 0);

    REQUIRE(bcf_hdr_add_sample(hdr.get(),"fa") == 0);
    REQUIRE(bcf_hdr_add_sample(hdr.get(),"mo") == 0);
    REQUIRE(bcf_hdr_add_sample(hdr.get(),"ch") == 0);

    bcf_hdr_sync(hdr.get());

    int ncontigs = 0;
    const char **contignames = bcf_hdr_seqnames(hdr.get(), &ncontigs);
    REQUIRE(ncontigs == 3);
    REQUIRE(hdr->n[BCF_DT_CTG] == 3);

    REQUIRE(string(contignames[0]) == "A");
    REQUIRE(hdr->id[BCF_DT_CTG][0].val != nullptr);
    REQUIRE(hdr->id[BCF_DT_CTG][0].val->info[0] == 1000000);

    REQUIRE(string(contignames[1]) == "B");
    REQUIRE(hdr->id[BCF_DT_CTG][1].val != nullptr);
    REQUIRE(hdr->id[BCF_DT_CTG][1].val->info[0] == 100000);

    REQUIRE(string(contignames[2]) == "C");
    REQUIRE(hdr->id[BCF_DT_CTG][2].val != nullptr);
    REQUIRE(hdr->id[BCF_DT_CTG][2].val->info[0] == 10000);

    free(contignames);

    REQUIRE(bcf_hdr_nsamples(hdr.get()) == 3);
    REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 0)) == "fa");
    REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 1)) == "mo");
    REQUIRE(string(bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, 2)) == "ch");
}


TEST_CASE("DNAnexus VCF/BCF serialization") {
    // load trio_denovo.vcf
    UPD(vcfFile, vcf, bcf_open("test/data/trio_denovo.vcf", "r"), [](vcfFile* f) { bcf_close(f); });
    UPD(bcf_hdr_t, hdr, bcf_hdr_read(vcf), &bcf_hdr_destroy);
    shared_ptr<bcf1_t> vt;
    vector<shared_ptr<bcf1_t>> records;

    bcf_hdr_append(hdr,"##contig=<ID=21,length=1000000>");
    bcf_hdr_sync(hdr);

    do {
        if (vt) {
            records.push_back(vt);
        }
        vt = shared_ptr<bcf1_t>(bcf_init(), &bcf_destroy);
    } while (bcf_read(vcf, hdr, vt.get()) == 0);

    REQUIRE(records.size() == 3);

    // serialize the BCF records into a memory buffer - without the header
    int memlen = 0;
    for (const auto& rec : records) {
        memlen += GLnexus::bcf_raw_calc_packed_len(rec.get());
    }
    std::cout << "memlen=" << memlen << std::endl;
    char *buf = (char*) malloc(memlen);
    REQUIRE(buf != NULL);

    int loc = 0;
    for (const auto& rec : records) {
        int reclen = GLnexus::bcf_raw_calc_packed_len(rec.get());
        GLnexus::bcf_raw_write_to_mem(rec.get(), reclen, &buf[loc]);
        loc += reclen;
    }

    // now read the records back from memory & ensure they match the originals
    loc = 0;
    for (const auto& rec : records) {
        loc += GLnexus::bcf_raw_read_from_mem(&buf[loc], vt.get());
        REQUIRE(vt->rid == rec->rid);
        REQUIRE(vt->pos == rec->pos);
        REQUIRE(vt->n_sample == rec->n_sample);
        REQUIRE(vt->n_allele == rec->n_allele);
        REQUIRE(vt->n_info == rec->n_info);
        REQUIRE(bcf_unpack(vt.get(), BCF_UN_ALL) == 0);
    }

    if (buf) {
        free(buf);
    }
}

/* Ccalculate the amount of bytes it would take to pack this bcf1 record.
 */
static int bcf_calc_packed_len(const std::shared_ptr<bcf1_t> &v)
{
    return 32 + v.get()->shared.l + v.get()->indiv.l;
}


/*
  Write the BCF record directly to a memory location. Return how much
  space was used.

  Note: the code is adapted from the bcf_write routine in htslib/vcf.c.
  The original prototype is:
      int bcf_write(htsFile *hfp, const bcf_hdr_t *h, bcf1_t *v)
  The original routine writes to a file, not to memory.
*/
static int bcf_write_to_mem(bcf1_t *v, char *addr) {
    int loc = 0;

    uint32_t x[8];
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

    return loc;
}

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


TEST_CASE("DNAnexus VCF/BCF serialization") {
    // load trio_denovo.vcf
    UPD(vcfFile, vcf, bcf_open("test/data/trio_denovo.vcf", "r"), [](vcfFile* f) { bcf_close(f); });
    UPD(bcf_hdr_t, hdr, bcf_hdr_read(vcf), &bcf_hdr_destroy);
    shared_ptr<bcf1_t> vt;
    vector<shared_ptr<bcf1_t>> records;

    bcf_hdr_append(hdr,"##contig=<ID=21,length=1000000>");
    bcf_hdr_sync(hdr);

    do {
        if (vt) {
            records.push_back(vt);
        }
        vt = shared_ptr<bcf1_t>(bcf_init(), &bcf_destroy);
    } while (bcf_read(vcf, hdr, vt.get()) == 0);

    REQUIRE(records.size() == 3);

    // serialize the BCF records into a memory buffer - without the header
    int memlen = 0;
    for (const auto& rec : records) {
        memlen += bcf_calc_packed_len(rec);
    }
    std::cout << "memlen=" << memlen << std::endl;
    char *buf = (char*) malloc(memlen);
    REQUIRE(buf != NULL);

    int loc = 0;
    for (const auto& rec : records) {
        loc += bcf_write_to_mem(rec.get(), &buf[loc]);
    }

    // now read the records back from memory & ensure they match the originals
    loc = 0;
    for (const auto& rec : records) {
        loc += bcf_read_from_mem(&buf[loc], vt.get());
        REQUIRE(vt->rid == rec->rid);
        REQUIRE(vt->pos == rec->pos);
        REQUIRE(vt->n_sample == rec->n_sample);
        REQUIRE(vt->n_allele == rec->n_allele);
        REQUIRE(vt->n_info == rec->n_info);
        REQUIRE(bcf_unpack(vt.get(), BCF_UN_ALL) == 0);
    }

    if (buf) {
        free(buf);
    }
}

TEST_CASE("htslib gVCF representation") {
    UPD(vcfFile, vcf, bcf_open("test/data/NA12878D_HiSeqX.21.10009462-10009469.gvcf", "r"), [](vcfFile* f) { bcf_close(f); });
    UPD(bcf_hdr_t, hdr, bcf_hdr_read(vcf), &bcf_hdr_destroy);
    shared_ptr<bcf1_t> vt;
    vector<shared_ptr<bcf1_t>> records;

    do {
        if (vt) {
            REQUIRE(bcf_unpack(vt.get(), BCF_UN_ALL) == 0);
            records.push_back(vt);
        }
        vt = shared_ptr<bcf1_t>(bcf_init(), &bcf_destroy);
    } while (bcf_read(vcf, hdr, vt.get()) == 0);

    REQUIRE(records.size() == 5);

    REQUIRE(records[0]->pos == 10009461);
    REQUIRE(records[0]->n_allele == 2);
    REQUIRE(string(records[0]->d.allele[0]) == "T");
    REQUIRE(string(records[0]->d.allele[1]) == "<NON_REF>");
    bcf_info_t *info = bcf_get_info(hdr, records[0].get(), "END");
    REQUIRE(info->type == BCF_BT_INT32);
    REQUIRE(info->len == 1);
    REQUIRE(info->v1.i == 10009463); // nb END stays 1-based!

    REQUIRE(records[1]->pos == 10009463);
    REQUIRE(records[1]->n_allele == 3);
    REQUIRE(string(records[1]->d.allele[0]) == "TA");
    REQUIRE(string(records[1]->d.allele[1]) == "T");
    REQUIRE(string(records[1]->d.allele[2]) == "<NON_REF>");
    REQUIRE(bcf_get_info(hdr, records[1].get(), "END") == NULL);

    REQUIRE(records[4]->pos == 10009468);
    REQUIRE(records[4]->n_allele == 2);
    REQUIRE(string(records[4]->d.allele[0]) == "A");
    REQUIRE(string(records[4]->d.allele[1]) == "<NON_REF>");
    REQUIRE(bcf_get_info(hdr, records[4].get(), "END")->v1.i == 10009471); // nb END stays 1-based!
}


/*
Ensure the code we've torn out remains functionally equivalent going
forward -- i.e. the test should break in the unlikely event a future
change in htslib causes its BCF [de]serialization to diverge from
ours.

Step 1: Read a VCF file, write it to a new file in uncompressed BCF
format. Also, store the results in a memory buffer using our own serialization routines.

Step 2: Read the BCF file, and compare that it is byte-for-byte equal
with the memory buffer.
*/
TEST_CASE("Ensure uncompressed BCF encoding remains consistent") {
    const char *tmp_bcf_file = "test/data/tmp.bcf";

    // Read a VCF file into an array of in-memory records
    vector<shared_ptr<bcf1_t>> records;
    {
        UPD(vcfFile, vcf, bcf_open("test/data/NA12878D_HiSeqX.21.10009462-10009469.gvcf", "r"), [](vcfFile* f) { bcf_close(f); });
        UPD(bcf_hdr_t, hdr, bcf_hdr_read(vcf), &bcf_hdr_destroy);
        shared_ptr<bcf1_t> vt;

        do {
            if (vt) {
                REQUIRE(bcf_unpack(vt.get(), BCF_UN_ALL) == 0);
                records.push_back(vt);
            }
            vt = shared_ptr<bcf1_t>(bcf_init(), &bcf_destroy);
        } while (bcf_read(vcf, hdr, vt.get()) == 0);
    }

    // Write the records to disk in uncompressed BCF format
    {
        std::remove(tmp_bcf_file);
        UPD(htsFile, fp, bcf_open(tmp_bcf_file, "wbu"),
            [](htsFile* f) { bcf_close(f); });
        shared_ptr<bcf_hdr_t> hdr(bcf_hdr_init("w"), &bcf_hdr_destroy);

        // Crafg a dummy VCF header. We can get away
        // with not matching the original file header, because it does not matter
        // (in this case).
        bcf_hdr_append(hdr.get(), "##fileDate=20090805");
        bcf_hdr_append(hdr.get(), "##FORMAT=<ID=UF,Number=1,Type=Integer,Description=\"Unused FORMAT\">");
        bcf_hdr_add_sample(hdr.get(), "NA00001");
        bcf_hdr_add_sample(hdr.get(), NULL);      // to update internal structures
        bcf_hdr_write(fp, hdr.get());

        REQUIRE(hts_get_format(fp)->compression == no_compression);
        //REQUIRE(hts_get_format(fp)->format == bcf);
        // force BCF format; not sure why this doesn't just work. Opening
        // the file in 'wbu' mode is supposed to
        fp->format.format = bcf;

        // write to the file
        for (const auto& rec : records) {
            bcf_write1(fp, hdr.get(), rec.get());
        }

        // note: the file is now closed, we can read it in the next phase.
    }

    // serialize the BCF records into a memory buffer  (without the header)
    int memlen = 0;
    char *buf;
    {
        for (const auto& rec : records) {
            memlen += GLnexus::bcf_raw_calc_packed_len(rec.get());
        }
        buf = (char*) malloc(memlen);
        REQUIRE(buf != NULL);

        int loc = 0;
        for (const auto& rec : records) {
            int reclen = GLnexus::bcf_raw_calc_packed_len(rec.get());
            GLnexus::bcf_raw_write_to_mem(rec.get(), reclen, &buf[loc]);
            loc += reclen;
        }
    }

    // compare the file data, while skipping the header
    //
    // read file data
    std::ifstream ifs(tmp_bcf_file);
    std::string content( (std::istreambuf_iterator<char>(ifs) ),
                         (std::istreambuf_iterator<char>()    ) );
    int filelen = content.length();
    const char *filebuf = content.c_str();
    std::cout << "filelen=" << filelen << " memlen=" << memlen << std::endl;
    assert(filelen >= memlen);

    // compare the last [memlen] bytes
    int rc = memcmp(&filebuf[filelen - memlen], buf, memlen);
    if (rc == 0) {
        std::cout << "BCF formats match" << std::endl;
    }
    else {
        std::cout << "BCF format mismatch" << std::endl;
    }
    REQUIRE(rc == 0);

    // cleanup
    if (buf != NULL) {
        free(buf);
    }
    std::remove(tmp_bcf_file);
}

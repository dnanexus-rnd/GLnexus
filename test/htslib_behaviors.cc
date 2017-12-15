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

TEST_CASE("memory bounds-checking in BCFSerialize reader") {
    // A BCF record that is too short, reading should fail
    char buf[16];
    int reclen = 0;
    shared_ptr<bcf1_t> vt = shared_ptr<bcf1_t>(bcf_init(), &bcf_destroy);
    GLnexus::Status s = GLnexus::bcf_raw_read_from_mem(buf, 16, sizeof(buf), vt.get(), reclen);
    REQUIRE(s == GLnexus::StatusCode::INVALID);
}

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
    //std::cout << "memlen=" << memlen << std::endl;
    char *buf = (char*) calloc(1, memlen);
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
        int reclen = 0;
        REQUIRE(GLnexus::bcf_raw_read_from_mem(buf, loc, memlen, vt.get(), reclen).ok());
        loc += reclen;
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
    REQUIRE(records[0]->rlen == 2); // htslib got this from END! it is not strlen(records[0]->d.allele[0])!

    REQUIRE(records[1]->pos == 10009463);
    REQUIRE(records[1]->n_allele == 3);
    REQUIRE(string(records[1]->d.allele[0]) == "TA");
    REQUIRE(string(records[1]->d.allele[1]) == "T");
    REQUIRE(string(records[1]->d.allele[2]) == "<NON_REF>");
    REQUIRE(bcf_get_info(hdr, records[1].get(), "END") == NULL);
    REQUIRE(records[1]->rlen == 2); // END is absent so this time it is strlen(records[1]->d.allele[0])

    REQUIRE(records[4]->pos == 10009468);
    REQUIRE(records[4]->n_allele == 2);
    REQUIRE(string(records[4]->d.allele[0]) == "A");
    REQUIRE(string(records[4]->d.allele[1]) == "<NON_REF>");
    REQUIRE(bcf_get_info(hdr, records[4].get(), "END")->v1.i == 10009471); // nb END stays 1-based!
    REQUIRE(records[4]->rlen == 3);
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
        buf = (char*) calloc(1, memlen);
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
    //std::cout << "filelen=" << filelen << " memlen=" << memlen << std::endl;
    assert(filelen >= memlen);

    // compare the last [memlen] bytes
    int rc = memcmp(&filebuf[filelen - memlen], buf, memlen);
    if (rc == 0) {
        //std::cout << "BCF formats match" << std::endl;
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

/*
Some bcf1_t accessor functions take a non-const bcf1_t*. We believe this is
because they need to "unpack" the record if it hasn't been already. Conversely
we assume that if the record is already unpacked, the accessor functions don't
in fact mutate the record, and they're safe to use from multiple threads. This
test is supposed to detect obvious violations of this assumption (but doesn't
prove its correctness)
*/
TEST_CASE("Check bcf1_t immutability") {
    const char *tmp_bcf_file = "test/data/tmp.bcf";

    // Read a VCF file into an array of in-memory records
    vector<shared_ptr<bcf1_t>> records;
    {
        UPD(vcfFile, vcf, bcf_open("test/data/NA12878D_HiSeqX.21.10009462-10009469.gvcf", "r"), [](vcfFile* f) { bcf_close(f); });
        UPD(bcf_hdr_t, hdr, bcf_hdr_read(vcf), &bcf_hdr_destroy);
        shared_ptr<bcf1_t> vt;

        do {
            if (vt) {
                records.push_back(vt);
            }
            vt = shared_ptr<bcf1_t>(bcf_init(), &bcf_destroy);
        } while (bcf_read(vcf, hdr, vt.get()) == 0);

        UPD(void, buf, malloc(sizeof(bcf1_t)), free);
        for (const auto& record : records) {
            // unpack and verify (positive control) the record changes in memory
            memcpy(buf, record.get(), sizeof(bcf1_t));
            REQUIRE(bcf_unpack(record.get(), BCF_UN_ALL) == 0);
            REQUIRE(memcmp(buf, record.get(), sizeof(bcf1_t)) != 0);
            memcpy(buf, record.get(), sizeof(bcf1_t));

            // call accessors and verify the record does NOT obviously change in memory
            int *gt_arr = nullptr, ngt_arr = 0;
            REQUIRE(bcf_get_genotypes(hdr, record.get(), &gt_arr, &ngt_arr) > 0);
            if (gt_arr) free(gt_arr);
            int32_t *v = nullptr;
            int vsz = 0;
            REQUIRE(bcf_get_format_int32(hdr, record.get(), "AD", &v, &vsz) != 0);
            if (v) free(v);
            REQUIRE(memcmp(buf, record.get(), sizeof(bcf1_t)) == 0);
        }
    }
}

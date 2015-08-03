// Test certain non-obvious behaviors of htslib 

#include <iostream>
#include <vcf.h>
#include <hfile.h>
#include <string.h>
#include <math.h>
#include <memory>
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

TEST_CASE("htslib hfile_mem VCF/BCF serialization") {
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
    char *buf = nullptr;
    size_t bufsz;
    char fn[4+sizeof(void**)+sizeof(size_t*)];
    char *pfn = fn;

    memcpy(pfn, "mem:", 4);
    *(char***)(pfn+4) = &buf;
    *(size_t**)(pfn+4+sizeof(void**)) = &bufsz;

    // bcf_open mode "u" is supposed to cause htslib to write "uncompressed
    // BCF", per the hts_open documentation in hts.h, but instead it writes
    // text VCF. Indeed the implementation of bcf_write in vcf.c hard-codes
    // bgzf_write, so it's incapable of writing out "uncompressed BCF" without
    // the BGZF container. So too for the BCF reading functions which hard-code
    // bgzf_read.
    //
    // We could use 0 to get BCF wrapped in the BGZF container format (with no
    // compression), but that would be less efficient, and carry some useless
    // block compression baggage.
    //
    // We may need to write a significant htslib patch to make it actually
    // able to read & write "uncompressed BCF"

    up_vcf.reset(bcf_open(pfn, "wu"));
    vcf = up_vcf.get();

    //REQUIRE(bcf_hdr_write(vcf, hdr) == 0);
    for (const auto& rec : records) {
        REQUIRE(bcf_write(vcf, hdr, rec.get()) == 0);
    }
    size_t memlen = htell(vcf->fp.hfile);
    REQUIRE(hflush(vcf->fp.hfile) == 0);

    //cout << string(buf, memlen);

    // now read the records back from memory & ensure they match the originals
    REQUIRE(bufsz >= memlen);
    bufsz = memlen;
    up_vcf.reset(bcf_open(pfn, "r"));
    vcf = up_vcf.get();

    REQUIRE(hts_get_format(vcf)->compression == no_compression);
    // TODO: without the header, htslib auto-detects the file format as SAM
    // instead of VCF. This is okay as long as we're using vcf_read below. But
    // if we use bcf_read, then because the file format is not set to VCF, it
    // tries to read it as bgzipped BCF.
    //REQUIRE(hts_get_format(vcf)->format == htsExactFormat::vcf);

    for (const auto& rec : records) {
        REQUIRE(vcf_read(vcf, hdr, vt.get()) == 0);
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
    REQUIRE(bcf_get_info(hdr, records[0].get(), "END")->v1.i == 10009463); // nb END stays 1-based!

    REQUIRE(records[1]->pos == 10009463);
    REQUIRE(records[1]->n_allele == 3);
    REQUIRE(string(records[1]->d.allele[0]) == "TA");
    REQUIRE(string(records[1]->d.allele[1]) == "T");
    REQUIRE(string(records[1]->d.allele[2]) == "<NON_REF>");

    REQUIRE(records[4]->pos == 10009468);
    REQUIRE(records[4]->n_allele == 2);
    REQUIRE(string(records[4]->d.allele[0]) == "A");
    REQUIRE(string(records[4]->d.allele[1]) == "<NON_REF>");
    REQUIRE(bcf_get_info(hdr, records[4].get(), "END")->v1.i == 10009471); // nb END stays 1-based!
}

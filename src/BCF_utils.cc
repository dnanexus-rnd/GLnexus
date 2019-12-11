#include "BCF_utils.h"

#include <set>

namespace GLnexus {

const uint64_t MAX_RECORD_LEN = 100000;

bool gvcf_compatible(const MetadataCache &metadata, const bcf_hdr_t *hdr) {
  Status s;
  const auto &contigs = metadata.contigs();

  // verify contigs match exactly. even the order matters

  int ncontigs = 0;
  const char **contignames = bcf_hdr_seqnames(hdr, &ncontigs);
  bool ans = true;

  if (((uint)ncontigs) != contigs.size()) {
    ans = false;
  } else {
    for (int i = 0; i < ncontigs; i++) {
      if (std::string(contignames[i]) != contigs[i].first ||
          hdr->id[BCF_DT_CTG][i].val->info[0] != contigs[i].second) {
        ans = false;
        break;
      }
    }
  }

  free(contignames);
  return ans;
}

bool xAtlas_ingestion_exceptions(const bcf_hdr_t *hdr, bcf1_t *bcf) {
  if (bcf_has_filter(hdr, bcf, "VRFromDeletion") == 1) {
    return true;
  }

  htsvecbox<int32_t> vr;
  int nVR = bcf_get_format_int32(hdr, bcf, "VR", &vr.v, &vr.capacity);
  if (nVR == 1 && vr[0] != bcf_int32_missing) {
    htsvecbox<int32_t> rr;
    int nRR = bcf_get_format_int32(hdr, bcf, "RR", &rr.v, &rr.capacity);
    if (nRR == 1 && rr[0] != bcf_int32_missing && vr[0] + rr[0] >= 65536) {
      return true;
    }
  }

  return false;
}

Status validate_bcf(const std::vector<std::pair<std::string, size_t> > &contigs,
                    const std::string &filename, const bcf_hdr_t *hdr,
                    bcf1_t *bcf, int prev_rid, int prev_pos,
                    bool &skip_ingestion) {
  skip_ingestion = false;
  if (bcf_unpack(bcf, BCF_UN_ALL) != 0 || bcf->errcode != 0) {
    return Status::Invalid(
        "invalid VCF record (corrupt format, or undeclared info/format fields)",
        filename + " " + range(bcf).str(contigs));
  }

  // A few hard-coded cases where we, reluctantly, skip ingestion
  // MAX_RECORD_LEN: blows up database (due to repetition across buckets)
  //                 and usually arises from gVCF caller bug anyway
  if (bcf->rlen >= MAX_RECORD_LEN || xAtlas_ingestion_exceptions(hdr, bcf)) {
    skip_ingestion = true;
    return Status::OK();
  }

  // Check that bcf->rlen is calculated correctly based on POS,END if
  // available or POS,strlen(REF) otherwise
  bcf_info_t *info = bcf_get_info(hdr, bcf, "END");
  if (info) {
    if (info->type != BCF_BT_INT8 && info->type != BCF_BT_INT16 &&
        info->type != BCF_BT_INT32) {
      return Status::Invalid("gVCF record's END field has unexpected type",
                             filename + " " + range(bcf).str(contigs));
    }
    if (info->len != 1) {
      return Status::Invalid("gVCF record has multiple END fields",
                             filename + " " + range(bcf).str(contigs));
    }
    if (info->v1.i < bcf->pos) {
      return Status::Invalid("gVCF record has END < POS",
                             filename + " " + range(bcf).str(contigs));
    }
    if (info->v1.i - bcf->pos != bcf->rlen) {
      return Status::Invalid("gVCF record END-POS doesn't match rlen",
                             filename + " " + range(bcf).str(contigs));
    }
  } else {
    if (bcf->d.allele == nullptr) {
      return Status::Invalid("gVCF allele is null",
                             filename + " " + range(bcf).str(contigs));
    }
    if (bcf->rlen != (int)strlen(bcf->d.allele[0])) {
      return Status::Invalid(
          "gVCF rlen doesn't match strlen(REF) (and no END field)",
          filename + " " + range(bcf).str(contigs));
    }
  }

  // verify record ordering is non-decreasing within a contig
  if (prev_rid == bcf->rid && prev_pos > bcf->pos) {
    return Status::Invalid("gVCF records are out-of-order ",
                           filename + " " + std::to_string(prev_pos + 1) +
                               " >= " + range(bcf).str(contigs));
  }

  // verify record does not go over the length of the contig
  const std::string& contig_name = contigs[bcf->rid].first;
  size_t contig_len = contigs[bcf->rid].second;
  if (bcf->pos + bcf->rlen > contig_len) {
    return Status::Invalid("gVCF record is longer than contig ",
                           filename + " " + range(bcf).str(contigs) + " " +
                               std::to_string(contig_len) + " " + contig_name);
  }

  // check that alleles are all distinct, and that all alleles are valid strings
  // of IUPAC nucleotides, except the last ALT allele which may be symbolic.
  std::set<std::string> alleles;
  for (int i = 0; i < bcf->n_allele; i++) {
    const std::string allele_i(bcf->d.allele[i]);
    if (!(is_iupac_nucleotides(allele_i) || allele_i == "*" ||
          (i == bcf->n_allele - 1 && is_symbolic_allele(allele_i.c_str())))) {
      return Status::Invalid(
          "allele is not a DNA sequence ",
          filename + " " + allele_i + " " + range(bcf).str(contigs));
    }
    alleles.insert(allele_i);
  }
  if (bcf->n_allele < 1 || alleles.size() != bcf->n_allele) {
    return Status::Invalid("alleles are not distinct ",
                           filename + " " + range(bcf).str(contigs));
  }

  // validate genotypes (all entries either missing or in [0, n_allele))
  if (bcf->n_sample != bcf_hdr_nsamples(hdr)) {
    return Status::Invalid("gVCF record doesn't have expected # of samples",
                           filename + " " + range(bcf).str(contigs));
  }
  htsvecbox<int> gt;
  int nGT = bcf_get_genotypes(hdr, bcf, &gt.v, &gt.capacity);
  if (nGT != 2 * bcf->n_sample &&
      /* Strelka2 accommodation: */ !(bcf->n_sample == 1 && nGT == 1)) {
    return Status::Invalid("gVCF record doesn't have expected # of GT entries",
                           filename + " " + range(bcf).str(contigs));
  }
  for (int i = 0; i < nGT; i++) {
    if (!bcf_gt_is_missing(gt[i]) &&
        (bcf_gt_allele(gt[i]) < 0 || bcf_gt_allele(gt[i]) >= bcf->n_allele)) {
      return Status::Invalid("invalid GT entry in gVCF record",
                             filename + " " + range(bcf).str(contigs));
    }
  }

  // validate genotype likelihoods (PL: all entries nonnegative)
  htsvecbox<int32_t> pl;
  int nPL = bcf_get_format_int32(hdr, bcf, "PL", &pl.v, &pl.capacity);
  if (nPL >= 0) {
    for (int i = 0; i < nPL && pl[i] != bcf_int32_vector_end; i++) {
      if (pl[i] != bcf_int32_missing && pl[i] < 0) {
        return Status::Invalid("negative PL entry in gVCF record",
                               filename + " " + range(bcf).str(contigs));
      }
    }
  }
  /*
      htsvecbox<float> gl;
      int nGL = bcf_get_format_float(hdr, bcf, "GL", &gl.v, &gl.capacity);
      if (nGL >= 0) {
          if (nGL != bcf->n_sample * diploid::genotypes(bcf->n_allele)) {
              return Status::Invalid("gVCF record doesn't have expected # of GL
     entries", filename + " " + range(bcf).str(contigs));
          }
          for (int i = 0; i < nGL; i++) {
              if (gl[i] > 0.0) {
                  return Status::Invalid("positive GL entry in gVCF record",
     filename + " " + range(bcf).str(contigs));
              }
          }
      }
  */
  return Status::OK();
}

Status vcf_validate_basic_facts(MetadataCache &metadata,
                                const std::string &dataset,
                                const std::string &filename, bcf_hdr_t *hdr,
                                vcfFile *vcf,
                                std::set<std::string> &samples_out) {
  if (!hdr) return Status::IOError("reading gVCF header", filename);
  if (!gvcf_compatible(metadata, hdr)) {
    return Status::Invalid(
        "Incompatible gVCF. The reference contigs must match the database "
        "configuration exactly.",
        filename);
  }

  std::vector<std::string> samples;
  unsigned n = bcf_hdr_nsamples(hdr);
  if (n == 0) {
    return Status::Invalid("gVCF contains no samples",
                           dataset + " (" + filename + ")");
  }
  for (unsigned i = 0; i < n; i++) {
    std::string sample(bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i));
    if (!regex_match(sample, regex_id)) {
      return Status::Invalid("gVCF contains invalid sample name",
                             dataset + " (" + filename + ")" + " " + sample);
    }
    samples.push_back(std::move(sample));
  }
  samples_out.clear();
  samples_out.insert(samples.begin(), samples.end());
  if (samples.size() != samples_out.size()) {
    return Status::Invalid("gVCF sample names are not unique",
                           dataset + " (" + filename + ")");
  }

  return Status::OK();
}

}  // namespace GLnexus

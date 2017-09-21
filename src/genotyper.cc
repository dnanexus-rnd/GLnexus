#include <assert.h>
#include <math.h>
#include <algorithm>
#include "genotyper.h"
#include "diploid.h"

using namespace std;

#include "genotyper_utils.h"

namespace GLnexus {

// Helper: determine if all genotypes in the record are 0/0
static inline bool is_homozygous_ref(const bcf_hdr_t* hdr, bcf1_t* record) {
    htsvecbox<int> gt;
    int nGT = bcf_get_genotypes(hdr, record, &gt.v, &gt.capacity);
    for (int i = 0; i < record->n_sample; i++) {
        assert(i*2+1<nGT);
        if (bcf_gt_is_missing(gt[i*2]) || bcf_gt_allele(gt[i*2]) != 0 ||
            bcf_gt_is_missing(gt[i*2+1]) || bcf_gt_allele(gt[i*2+1]) != 0) {
            return false;
        }
    }
    return true;
}

// Helper: given REF and ALT DNA, determine if the ALT represents a deletion
// with respect to REF. Left-alignment is assumed and reference padding on the
// left is tolerated.
static inline bool is_deletion(const string& ref, const string& alt) {
    if (alt.size() >= ref.size()) {
        return false;
    }
    for (int j = 0; j < alt.size(); j++) {
        if (alt[j] != ref[j]) {
            // some kind of complex edit where the ALT allele is both shorter
            // and with differing basis in the prefix...
            return false;
        }
    }
    return true;
}

// Pre-process a bcf1_t record to cache some useful info that we'll use repeatedly
 Status preprocess_record(const unified_site& site, const bcf_hdr_t* hdr, const shared_ptr<bcf1_t>& record,
                          bcf1_t_plus& ans) {
    range rng(record);
    assert(rng.overlaps(site.pos));

    ans.p = record;

    ans.is_ref = is_gvcf_ref_record(record.get());

    if(bcf_get_genotypes(hdr, record.get(), &ans.gt.v, &ans.gt.capacity) != 2*record->n_sample || !ans.gt.v) {
        return Status::Failure("genotyper::preprocess_record: unexpected result from bcf_get_genotypes");
    }

    ans.allele_mapping.assign(record->n_allele, -1);
    ans.allele_mapping[0] = 0;
    ans.deletion_allele.assign(record->n_allele, false);
    ans.has_lost_allele = false;

    // map the bcf1_t alt alleles according to unification
    // checking for valid dna regex match
    string ref_al(record->d.allele[0]);
    for (int i = 1; i < record->n_allele; i++) {
        string al(record->d.allele[i]);
        if (regex_match(al, regex_dna)) {
            auto p = site.unification.find(allele(rng, al));
            if (p != site.unification.end()) {
                ans.allele_mapping[i] = p->second;
            } else {
                ans.has_lost_allele = true;
            }
        }
        if (al.size() < rng.size() && rng.size() == ref_al.size()) {
            ans.deletion_allele[i] = is_deletion(ref_al, al);
        }
    }

    return Status::OK();
}

///////////////////////////////////////////////////////////////////////////////
// Genotyper core
///////////////////////////////////////////////////////////////////////////////

/// Revise genotypes which are initially called with lost alleles. Frequently these
/// are low quality and "round down" to homozygous ref. Mutates record.
Status revise_genotypes(const genotyper_config& cfg, const unified_site& us, const map<int, int>& sample_mapping,
                        const bcf_hdr_t* hdr, bcf1_t_plus& vr) {
    assert(!vr.is_ref);
    if (!vr.has_lost_allele) {
        // below would be a no-op anyway
        return Status::OK();
    }

    // start by replacing the record with a duplicate, since it may not be safe to
    // mutate the "original"
    auto record = shared_ptr<bcf1_t>(bcf_dup(vr.p.get()), &bcf_destroy);
    vr.p = record;
    if (bcf_unpack(record.get(), BCF_UN_ALL)) return Status::Failure("genotyper::prepare_dataset_records bcf_unpack");
    unsigned nGT = diploid::genotypes(record->n_allele);
    range rng(record);

    // extract input genotype likelihoods and GQ
    vector<double> gll;
    Status s = diploid::bcf_get_genotype_log_likelihoods(hdr, record.get(), gll);
    if (!s.ok()) {
        return Status::Failure("genotyper::revise_genotypes: couldn't find genotype likelihoods in gVCF record", s.str());
    }
    assert(gll.size() == nGT*record->n_sample);
    htsvecbox<int32_t> gq;
    if(bcf_get_format_int32(hdr, record.get(), "GQ", &gq.v, &gq.capacity) != record->n_sample || !gq.v) {
        return Status::Failure("genotyper::revise_genotypes: unexpected result from bcf_get_format_int32 GQ");
    }
    assert(gq.capacity >= record->n_sample);

    // construct "prior" over input ALT alleles which penalizes lost ones (otherwise flat)
    const float lost_log_prior = log(std::max(us.lost_allele_frequency, cfg.min_assumed_allele_frequency));
    vector<double> gt_log_prior(diploid::genotypes(record->n_allele), 0.0);
    for (int i = 0; i < record->n_allele; i++) {
        if (vr.allele_mapping[i] == -1) {
            assert(i > 0);
            gt_log_prior[i] = lost_log_prior;
        }
    }

    // proceed through designated samples
    for (const auto& sample : sample_mapping) {
        assert(sample.first < record->n_sample);
        // add "priors" to genotype likelihoods; keep track of MAP and 2nd (silver)
        double* sample_gll = gll.data() + sample.first*nGT;
        double map_gll = log(0), silver_gll = log(0);
        int map_gt = -1;
        for (int g = 0; g < nGT; g++) {
            const auto alleles = diploid::gt_alleles(g);
            // Use the smaller of the priors on the two alleles and not their product.
            // If we view this as "penalizing" the likelihoods of genotypes which include
            // lost alleles, one such penalty is sufficient.
            auto g_ll = sample_gll[g] + std::min(gt_log_prior[alleles.first], gt_log_prior[alleles.second]);
            if (g_ll > map_gll) {
                silver_gll = map_gll;
                map_gll = g_ll;
                map_gt = g;
            } else if (g_ll > silver_gll) {
                silver_gll = g_ll;
            }
        }
        assert(map_gt >= 0 && map_gt < nGT);
        assert(map_gll >= silver_gll);
        assert(silver_gll > log(0));

        // record MAP genotype and recalculate GQ
        const auto revised_alleles = diploid::gt_alleles(map_gt);
        vr.gt.v[sample.first*2] = bcf_gt_unphased(revised_alleles.first);
        vr.gt.v[sample.first*2+1] = bcf_gt_unphased(revised_alleles.second);
        gq.v[sample.first] = std::min(99, (int) round(10.0*(map_gll - silver_gll)/log(10.0)));
    }

    // write GT and GQ back into record
    if (bcf_update_format_int32(hdr, record.get(), "GQ", gq.v, record->n_sample)) {
        return Status::Failure("genotyper::revise_genotypes: bcf_update_format_int32 GQ failed");
    }
    if (bcf_update_genotypes(hdr, record.get(), vr.gt.v, 2*record->n_sample)) {
        return Status::Failure("genotyper::revise_genotypes: bcf_update_genotypes failed");
    }

    return Status::OK();
}

/// Given a unified site and the set of gVCF records overlapping it in some
/// dataset, check that they span the site, preprocess them, and separate
/// the reference and variant records.
///
/// Ideally and often there's either zero or one variant records -- zero if
/// the dataset exhibits no variation at the site, and one if it does.
/// Unfortunately, there are a number of circumstances under which some
/// variant callers produce multiple overlapping records.
///
/// The variant records, if any, are returned through variant_records. It is
/// then the job of translate_genotypes to figure out what to do with the
/// cluster of variant records.
///
/// This function also modifies min_ref_depth which tracks the minimum depth
/// of ref-records for a given sample. min_ref_depth is expected to be
/// initialized to -1 for all samples, and the entry for a sample will be
/// modified only if one or more ref records are present.
///
/// =====================================================
/// Pre conditions:
/// min_ref_depth should be initialized to -1 for all samples
/// ======================================================
/// Post conditions:
/// No records given
///      variant_records empty, min_ref_depth unchanged, rnc = MissingData
/// Records do not span the entire range of site
///      variant_records empty, min_ref_depth updated accordingly, rnc = PartialData
/// Records span entire range, and consist of all reference confidence records
///      variant_records empty, min_ref_depth updated accordingly, rnc = N_A
/// Records span entire range, and include one or more variant records which
///      variant_records filled in, min_ref_depth updated accordingly, rnc = N_A
///
///
/// FIXME: detect & complain if the reference confidence records actually overlap the
///        variant records
Status prepare_dataset_records(const genotyper_config& cfg, const unified_site& site,
                               const string& dataset, const bcf_hdr_t* hdr, int bcf_nsamples,
                               const map<int, int>& sample_mapping,
                               const vector<shared_ptr<bcf1_t>>& records,
                               AlleleDepthHelper& depth,
                               NoCallReason& rnc,
                               vector<int>& min_ref_depth,
                               vector<shared_ptr<bcf1_t_plus>>& all_records,
                               vector<shared_ptr<bcf1_t_plus>>& variant_records) {
    // initialize outputs
    rnc = NoCallReason::MissingData;
    all_records.clear();
    variant_records.clear();

    Status s;

    // collect the ranges covered by the records
    vector<range> record_rngs;
    record_rngs.reserve(records.size());
    for (const auto& record: records) {
        range record_rng(record.get());
        assert(record_rng.overlaps(site.pos));
        record_rngs.push_back(record_rng);
        assert(bcf_nsamples == record->n_sample);
    }

    // check the records span the site, otherwise we need to produce
    // PartialData no-calls
    if (!site.pos.spanned_by(record_rngs)) {
        if (!record_rngs.empty()) {
            rnc = NoCallReason::PartialData;
        }
        return Status::OK();
    }

    vector<shared_ptr<bcf1_t_plus>> ref_records;
    for (const auto& record : records) {
        auto rp = make_shared<bcf1_t_plus>();
        S(preprocess_record(site, hdr, record, *rp));
        if (rp->is_ref) {
            ref_records.push_back(rp);
        } else {
            if (cfg.revise_genotypes) {
                S(revise_genotypes(cfg, site, sample_mapping, hdr, *rp));
            }
            variant_records.push_back(rp);
        }
        all_records.push_back(rp);
    }

    // compute min_ref_depth across the reference confidence records
    S(update_min_ref_depth(dataset, hdr, bcf_nsamples, sample_mapping,
                           ref_records, depth, min_ref_depth));

    // Success...
    rnc = NoCallReason::N_A;
    return Status::OK();
}

/// Based on the cluster of variant records and min_ref_depth produced by
/// prepare_dataset_records, fill genotypes for this dataset's samples with
/// appropriate calls (currently by translation of the input hard-calls).
/// Updates genotypes and may modify min_ref_depth.
///
/// FIXME: not coded to deal with multi-sample gVCFs properly.
static Status translate_genotypes(const genotyper_config& cfg, const unified_site& site,
                                  const string& dataset, const bcf_hdr_t* dataset_header,
                                  int bcf_nsamples, const map<int,int>& sample_mapping,
                                  const vector<shared_ptr<bcf1_t_plus>>& variant_records,
                                  AlleleDepthHelper& depth,
                                  vector<int>& min_ref_depth,
                                  vector<one_call>& genotypes) {
    assert(genotypes.size() == 2*min_ref_depth.size());
    assert(!site.monoallelic);
    Status s;

    // Scan the variant records to pull out those with 0/0 genotype calls
    // from those actually presenting variation
    vector<shared_ptr<bcf1_t_plus>> records_00, records_non00;
    for (const auto& a_record : variant_records) {
        assert(!a_record->is_ref);
        bool non00 = false;
        for (int i = 0; i < 2*bcf_nsamples; i++) {
            assert(a_record->gt.capacity > i);
            if (bcf_gt_is_missing(a_record->gt[i]) || bcf_gt_allele(a_record->gt[i]) != 0) {
                non00 = true;
            }
        }
        if (non00) {
            records_non00.push_back(a_record);
        } else {
            records_00.push_back(a_record);
        }
    }

    // update min_ref_depth with found 0/0 records
    S(update_min_ref_depth(dataset, dataset_header, bcf_nsamples, sample_mapping,
                           records_00, depth, min_ref_depth));

    bcf1_t_plus* record = nullptr;
    bool half_call = false;

    if (records_non00.size() == 0) {
        // no variation represented in this dataset; make homozygous ref calls
        // if min_ref_depth indicates sufficient coverage for this sample
        for (const auto& ij : sample_mapping) {
            assert(ij.second < min_ref_depth.size());
            int rd = min_ref_depth[ij.second];

            if (rd >= cfg.required_dp) {
                genotypes[2*ij.second] =
                    genotypes[2*ij.second+1] =
                        one_call(bcf_gt_unphased(0), NoCallReason::N_A);
            } else {
                genotypes[2*ij.second].RNC =
                    genotypes[2*ij.second+1].RNC = NoCallReason::InsufficientDepth;
            }
        }
        return Status::OK();
    } else if (records_non00.size() == 1)  {
        // simple common case: one variant record overlapping the unified site
        record = records_non00[0].get();
    } else {
        // complex situation: multiple non-0/0 records overlapping the unified site.
        //
        // If the records don't all share at least one position in common (i.e. their
        // ALT alleles aren't mutually exclusive on one chromosome), punt with
        // UnphasedVariants. We'll improve this in the future.
        //
        // If at least one record is a heterozygous 0/X call where X is a known
        // allele in the unified site, and none of the records call >1 ALT allele,
        // generate a half-call from the highest-quality such record.
        // This at least recovers some of the information when the GVCF has two
        // overlapping 0/X records for one sample (we'd rather it present one record
        // heterozygous for two ALTs)
        //
        // Otherwise: punt with OverlappingVariants

        half_call = true;
        range intersection(records_non00[0]->p);
        for (auto& a_record : records_non00) {
            range record_rng(a_record->p);
            assert(record_rng.rid == intersection.rid);
            intersection.beg = max(record_rng.beg, intersection.beg);
            intersection.end = min(record_rng.end, intersection.end);

            for (int i = 0; half_call && i < bcf_nsamples; i++) {
                assert(a_record->gt.capacity > 2*i);
                if (!bcf_gt_is_missing(a_record->gt[2*i]) && bcf_gt_allele(a_record->gt[2*i]) != 0) {
                    half_call = false;
                } else if (!bcf_gt_is_missing(a_record->gt[2*i+1])) {
                    auto al = bcf_gt_allele(a_record->gt[2*i+1]);
                    if (al > 0 && a_record->allele_mapping[al] > 0 &&
                        (!record || record->p->qual < a_record->p->qual)) {
                        record = a_record.get();
                    }
                }
            }
        }

        if (intersection.beg >= intersection.end) {
            for (int i = 0; i < bcf_nsamples; i++) {
                genotypes[sample_mapping.at(i)*2].RNC =
                    genotypes[sample_mapping.at(i)*2+1].RNC =
                        NoCallReason::UnphasedVariants;
            }
            return Status::OK();
        }

        if (!record || !half_call) {
            for (int i = 0; i < bcf_nsamples; i++) {
                genotypes[sample_mapping.at(i)*2].RNC =
                    genotypes[sample_mapping.at(i)*2+1].RNC =
                        NoCallReason::OverlappingVariants;
            }
            return Status::OK();
        }
    }

    // Now, translating genotypes from one variant BCF record.
    assert(record != nullptr);

    // get the genotype calls
    htsvecbox<int> gt;
    int nGT = bcf_get_genotypes(dataset_header, record->p.get(), &gt.v, &gt.capacity);
    int n_bcf_samples = bcf_hdr_nsamples(dataset_header);
    if (!gt.v || nGT != 2*n_bcf_samples) return Status::Failure("genotyper::translate_genotypes bcf_get_genotypes");
    assert(record->p->n_sample == bcf_hdr_nsamples(dataset_header));

    S(depth.Load(dataset, dataset_header, record->p.get()));

    // for each shared sample, record the genotype call.
    for (const auto& ij : sample_mapping) {
        assert(2*ij.first < nGT);
        assert(ij.second < min_ref_depth.size());

        // TODO: are depth and allele_mapping checks inside-out????
        #define fill_allele(ofs)                                                  \
            if (gt[2*ij.first+ofs] != bcf_int32_vector_end &&                     \
                !bcf_gt_is_missing(gt[2*ij.first+(ofs)])) {                       \
                auto al = bcf_gt_allele(gt[2*ij.first+(ofs)]);                    \
                assert(al >= 0 && al < record->p->n_allele);                      \
                int rd = min_ref_depth[ij.second];                                \
                if (depth.get(ij.first, al) >= cfg.required_dp                    \
                    && (rd < 0 || rd >= cfg.required_dp)) {                       \
                    if (record->allele_mapping[al] >= 0) {                        \
                        genotypes[2*ij.second+(ofs)] =                            \
                            one_call(bcf_gt_unphased(record->allele_mapping[al]), \
                                     NoCallReason::N_A);                          \
                    } else {                                                      \
                        genotypes[2*ij.second+(ofs)].RNC =                        \
                            record->deletion_allele[al]                           \
                                ? NoCallReason::LostDeletion                      \
                                : NoCallReason::LostAllele;                       \
                    }                                                             \
                } else {                                                          \
                    genotypes[2*ij.second+(ofs)].RNC =                            \
                        NoCallReason::InsufficientDepth;                          \
                }                                                                 \
            }

        if (half_call) {
            genotypes[2*ij.second].RNC = NoCallReason::OverlappingVariants;
        } else {
            fill_allele(0)
        }
        fill_allele(1)
    }

    return Status::OK();
}

/// streamlined version of translate_genotypes for monoallelic sites
/// FIXME: not coded to deal with multi-sample gVCFs properly.
static Status translate_monoallelic(const genotyper_config& cfg, const unified_site& site,
                                    const string& dataset, const bcf_hdr_t* dataset_header,
                                    int bcf_nsamples, const map<int,int>& sample_mapping,
                                    const vector<shared_ptr<bcf1_t_plus>>& variant_records,
                                    AlleleDepthHelper& depth,
                                    vector<int>& min_ref_depth,
                                    vector<one_call>& genotypes) {
    assert(genotypes.size() == 2*min_ref_depth.size());
    assert(site.monoallelic);
    assert(site.alleles.size() == 2);
    Status s;
    bcf1_t_plus* record = nullptr;

    // scan the (potentially) multiple overlapping records to find the one with
    // the desired allele
    for (auto& a_record : variant_records) {
        assert(!a_record->is_ref);
        for (unsigned al = 1; al < a_record->allele_mapping.size(); al++) {
            if (a_record->allele_mapping[al] == 1) {
                if (record == nullptr) {
                    record = a_record.get();
                } else {
                    // uh, overlapping records with the same allele!?
                    for (int i = 0; i < bcf_nsamples; i++) {
                        genotypes[sample_mapping.at(i)*2].RNC =
                            genotypes[sample_mapping.at(i)*2+1].RNC =
                                NoCallReason::OverlappingVariants;
                    }
                    return Status::OK();
                }
            }
        }
    }

    if (record == nullptr) {
        // we have nothing to say here
        for (int i = 0; i < bcf_nsamples; i++) {
            genotypes[sample_mapping.at(i)*2].RNC =
                genotypes[sample_mapping.at(i)*2+1].RNC =
                    NoCallReason::MonoallelicSite;
        }
        return Status::OK();
    }

    // Now, translating genotypes from one variant BCF record.

    // get the genotype calls
    htsvecbox<int> gt;
    int nGT = bcf_get_genotypes(dataset_header, record->p.get(), &gt.v, &gt.capacity);
    int n_bcf_samples = bcf_hdr_nsamples(dataset_header);
    if (!gt.v || nGT != 2*n_bcf_samples) return Status::Failure("genotyper::translate_genotypes bcf_get_genotypes");
    assert(record->p->n_sample == bcf_hdr_nsamples(dataset_header));

    S(depth.Load(dataset, dataset_header, record->p.get()));

    // for each shared sample, record the genotype call.
    for (const auto& ij : sample_mapping) {
        assert(2*ij.first < nGT);
        assert(ij.second < min_ref_depth.size());

        #define fill_monoallelic(ofs)                                             \
            if (gt[2*ij.first+ofs] != bcf_int32_vector_end &&                     \
                !bcf_gt_is_missing(gt[2*ij.first+(ofs)])) {                       \
                auto al = bcf_gt_allele(gt[2*ij.first+(ofs)]);                    \
                assert(al >= 0 && al < record->p->n_allele);                      \
                if (record->allele_mapping[al] > 0) {                             \
                    if (depth.get(ij.first, al) >= cfg.required_dp) {             \
                        genotypes[2*ij.second+(ofs)] =                            \
                            one_call(bcf_gt_unphased(record->allele_mapping[al]), \
                                     NoCallReason::N_A);                          \
                    } else {                                                      \
                        genotypes[2*ij.second+(ofs)].RNC =                        \
                            NoCallReason::InsufficientDepth;                      \
                    }                                                             \
                } else {                                                          \
                    genotypes[2*ij.second+(ofs)].RNC =                            \
                        NoCallReason::MonoallelicSite;                            \
                }                                                                 \
            }

        fill_monoallelic(0)
        fill_monoallelic(1)
    }

    return Status::OK();
}

Status genotype_site(const genotyper_config& cfg, MetadataCache& cache, BCFData& data, const unified_site& site,
                     const std::string& sampleset, const vector<string>& samples,
                     const bcf_hdr_t* hdr, shared_ptr<bcf1_t>& ans,
                     bool residualsFlag, shared_ptr<string> &residual_rec,
                     atomic<bool>* ext_abort) {
    Status s;

    // Initialize a vector for the unified genotype calls for each sample,
    // starting with everything missing. We'll then loop through BCF records
    // overlapping this site and fill in the genotypes as we encounter them.
    vector<one_call> genotypes(2*samples.size());

    // Setup format field helpers
    vector<unique_ptr<FormatFieldHelper>> format_helpers;
    S(setup_format_helpers(format_helpers, cfg.liftover_fields, site, samples));

    shared_ptr<const set<string>> samples2, datasets;
    vector<unique_ptr<RangeBCFIterator>> iterators;
    S(data.sampleset_range(cache, sampleset, site.pos, nullptr,
                           samples2, datasets, iterators));
    assert(samples.size() == samples2->size());

    AlleleDepthHelper adh(cfg);
    vector<DatasetResidual> lost_calls_info;

    map<string,int> samples_index;
    for (int i = 0; i < samples.size(); i++) {
        assert(samples_index.find(samples[i]) == samples_index.end());
        samples_index[samples[i]] = i;
    }

    // for each pertinent dataset
    for (const auto& dataset : *datasets) {
        if (ext_abort && *ext_abort) {
            return Status::Aborted();
        }

        // load BCF records overlapping the site by "merging" the iterators
        shared_ptr<const bcf_hdr_t> dataset_header;
        vector<vector<shared_ptr<bcf1_t>>> recordss(iterators.size());
        vector<shared_ptr<bcf1_t>> records;

        for (const auto& iter : iterators) {
            string this_dataset;
            vector<shared_ptr<bcf1_t>> these_records;
            S(iter->next(this_dataset, dataset_header, these_records));
            if (dataset != this_dataset) {
                return Status::Failure("genotype_site: iterator returned unexpected dataset",
                                       this_dataset + " instead of " + dataset);
            }
            records.insert(records.end(), these_records.begin(), these_records.end());
        }

        assert(is_sorted(records.begin(), records.end(),
                         [] (shared_ptr<bcf1_t>& p1, shared_ptr<bcf1_t>& p2) {
                            return range(p1) < range(p2);
                         }));

        // index the samples shared between the sample set and the BCFs.
        // this could be cached on sampleset/dataset cross
        map<int,int> sample_mapping;
        int bcf_nsamples = bcf_hdr_nsamples(dataset_header.get());
        for (int i = 0; i < bcf_nsamples; i++) {
            string sample_i(bcf_hdr_int2id(dataset_header.get(), BCF_DT_SAMPLE, i));
            const auto p = samples_index.find(sample_i);
            if (p != samples_index.end()) {
                sample_mapping[i] = p->second;
            }
        }

        // pre-process the records
        vector<int> min_ref_depth(samples.size(), -1);
        vector<shared_ptr<bcf1_t_plus>> all_records, variant_records;
        NoCallReason rnc = NoCallReason::MissingData;
        S(prepare_dataset_records(cfg, site, dataset, dataset_header.get(), bcf_nsamples,
                                  sample_mapping, records, adh, rnc, min_ref_depth,
                                  all_records, variant_records));

        if (rnc != NoCallReason::N_A) {
            // no call for the samples in this dataset (several possible
            // reasons)
            for (const auto& p : sample_mapping) {
                genotypes[p.second*2].RNC =
                    genotypes[p.second*2+1].RNC = rnc;
            }
        } else if (!site.monoallelic) {
            // make genotype calls for the samples in this dataset
            S(translate_genotypes(cfg, site, dataset, dataset_header.get(), bcf_nsamples,
                                  sample_mapping, variant_records, adh, min_ref_depth,
                                  genotypes));
        } else {
            S(translate_monoallelic(cfg, site, dataset, dataset_header.get(), bcf_nsamples,
                                    sample_mapping, variant_records, adh, min_ref_depth,
                                    genotypes));
        }

        // Update FORMAT fields for this dataset.
        S(update_format_fields(dataset, dataset_header.get(), sample_mapping, site, format_helpers,
                               all_records, variant_records));
        // But if rnc = MissingData, PartialData, UnphasedVariants, or OverlappingVariants, then
        // we must censor the FORMAT fields as potentially unreliable/misleading.
        for (const auto& p : sample_mapping) {
            auto rnc1 = genotypes[p.second*2].RNC;
            auto rnc2 = genotypes[p.second*2+1].RNC;

            if (rnc1 == NoCallReason::MissingData || rnc1 == NoCallReason::PartialData) {
                assert(rnc1 == rnc2);
                for (const auto& fh : format_helpers) {
                    S(fh->censor(p.second, false));
                }
            } else if (site.monoallelic || rnc1 == NoCallReason::UnphasedVariants || rnc1 == NoCallReason::OverlappingVariants) {
                const bool half_call = site.monoallelic || ((genotypes[p.second*2].RNC == NoCallReason::N_A) != (genotypes[p.second*2+1].RNC == NoCallReason::N_A));

                for (const auto& fh : format_helpers) {
                    if (fh->field_info.name != "DP") {
                        S(fh->censor(p.second, half_call));
                    }
                }
            }
        }

        // Handle residuals
        if (residualsFlag) {
            // TODO: don't emit residuals for lost alleles which will be represented in
            // a separate monoallelic site
            const set<NoCallReason> non_residual_RNCs = { NoCallReason::N_A, NoCallReason::MissingData,
                                                          NoCallReason::PartialData, NoCallReason::InsufficientDepth,
                                                          NoCallReason::MonoallelicSite };

            bool any_lost_calls = false;
            for (int i = 0; i < bcf_nsamples; i++) {
                if (non_residual_RNCs.find(genotypes[sample_mapping.at(i)*2].RNC) == non_residual_RNCs.end() ||
                    non_residual_RNCs.find(genotypes[sample_mapping.at(i)*2 + 1].RNC) == non_residual_RNCs.end()) {
                    any_lost_calls = true;
                    break;
                }
            }

            if (any_lost_calls) {
                // missing call, keep it in memory
                DatasetResidual dsr;
                dsr.name = dataset;
                dsr.header = dataset_header;
                dsr.records = records;
                lost_calls_info.push_back(dsr);
            }
        }
    }

    // Clean up emission order of alleles
    for(size_t i=0; i < samples.size(); i++) {
        if (genotypes[2*i].allele != bcf_gt_missing && genotypes[2*i+1].allele == bcf_gt_missing ||
            genotypes[2*i].allele > genotypes[2*i+1].allele) {
            swap(genotypes[2*i], genotypes[2*i+1]);
        }
    }
    // Create the destination BCF record for this site.
    ans = shared_ptr<bcf1_t>(bcf_init(), &bcf_destroy);
    ans->rid = site.pos.rid;
    ans->pos = site.pos.beg;
    ans->rlen = site.pos.end - site.pos.beg;
    ans->qual = site.qual;

    // alleles
    vector<const char*> c_alleles;
    for (const auto& allele : site.alleles) {
        c_alleles.push_back(allele.c_str());
    }
    if (bcf_update_alleles(hdr, ans.get(), c_alleles.data(), c_alleles.size()) != 0) {
        return Status::Failure("bcf_update_alleles");
    }

    // GT
    vector<int32_t> gt;
    for (const auto& c : genotypes) {
        gt.push_back(c.allele);
    }
    assert(gt.size() == genotypes.size());
    if (bcf_update_genotypes(hdr, ans.get(), gt.data(), gt.size()) != 0) {
        return Status::Failure("bcf_update_genotypes");
    }

    // Lifted-over FORMAT fields (non-genotype based)
    for (auto& format_helper : format_helpers) {
        S(format_helper->update_record_format(hdr, ans.get()));
    }

    // RNC
    vector<const char*> rnc;
    for (const auto& c : genotypes) {
        char* v = (char*) "M";
        #define RNC_CASE(reason,code) case NoCallReason::reason: v = (char*) code ; break;
        switch (c.RNC) {
            RNC_CASE(N_A,".")
            RNC_CASE(PartialData,"P")
            RNC_CASE(InsufficientDepth,"D")
            RNC_CASE(LostDeletion,"-")
            RNC_CASE(LostAllele,"L")
            RNC_CASE(UnphasedVariants,"U")
            RNC_CASE(OverlappingVariants,"O")
            RNC_CASE(MonoallelicSite,"M")
            default:
                assert(c.RNC == NoCallReason::MissingData);
        }
        rnc.push_back(v);
    }
    assert (gt.size() == rnc.size());
    if (bcf_update_format_string(hdr, ans.get(), "RNC", rnc.data(), rnc.size()) != 0) {
        return Status::Failure("bcf_update_format_string RNC");
    }

    if (site.monoallelic && bcf_add_filter(hdr, ans.get(), bcf_hdr_id2int(hdr, BCF_DT_ID, "MONOALLELIC")) != 1) {
        return Status::Failure("bcf_add_filter MONOALLELIC");
    }

    if (residualsFlag &&
        !lost_calls_info.empty()) {
        // Write loss record to the residuals file, useful for offline debugging.
        residual_rec = make_shared<string>();
        S(residuals_gen_record(site, hdr, ans.get(), lost_calls_info,
                               cache, samples,
                               *residual_rec));
    }

    // Overwrite the output BCF record with a duplicate. Why? This forces htslib to
    // perform some internal serialization of the data (see the static bcf1_sync
    // function in vcf.c, which we can't call directly, but is called by bcf_dup).
    // htslib would otherwise do this serialization implicitly while writing the
    // record out to a file, but by doing it explicitly here, we get to do some of the
    // work in the current worker thread rather than the single thread responsible for
    // writing out the file.
    auto ans2 = shared_ptr<bcf1_t>(bcf_dup(ans.get()), &bcf_destroy);
    ans = move(ans2);

    return Status::OK();
}

}

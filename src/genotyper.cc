#include <assert.h>
#include <math.h>
#include <algorithm>
#include "genotyper.h"
#include "diploid.h"
#include "vcfutils.h"

using namespace std;

#include "genotyper_utils.h"

// Prevent dependency on unnecessarily new version of glibc/libm
// https://stackoverflow.com/a/5977518
__asm__(".symver logf,logf@GLIBC_2.2.5");

namespace GLnexus {

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
    assert(rng.rid == site.pos.rid);

    ans.p = record;

    ans.is_ref = is_gvcf_ref_record(record.get());

    auto nGT = bcf_get_genotypes(hdr, record.get(), &ans.gt.v, &ans.gt.capacity);
    if (record->n_sample == 1 && nGT == 1 && !ans.gt.empty()) {
        // special case for Strelka2 and other callers which emit some gVCF
        // records with GT=. or GT=0 or GT=1: rewrite these to look like ./.
        // and ./0 and ./1 as far as our genotyper is concerned.
        ans.gt.v = (int*) realloc(ans.gt.v, 2*sizeof(int));
        ans.gt.capacity = 2;
        swap(ans.gt[0], ans.gt[1]);
        ans.gt[0] = bcf_gt_missing;
        assert(bcf_gt_is_missing(ans.gt[0]));
        ans.was_haploid = true;
    } else if(nGT != 2*record->n_sample || !ans.gt.v) {
        return Status::Failure("genotyper::preprocess_record: unexpected result from bcf_get_genotypes");
    }

    ans.allele_mapping.assign(record->n_allele, -1);
    ans.allele_mapping[0] = 0;
    ans.deletion_allele.assign(record->n_allele, false);

    // map the bcf1_t alt alleles according to unification
    // checking for valid dna regex match
    string ref_al(record->d.allele[0]);
    for (int i = 1; i < record->n_allele; i++) {
        string al(record->d.allele[i]);
        if (is_dna(al)) {
            auto p = site.unification.find(allele(rng, al));
            if (p != site.unification.end()) {
                ans.allele_mapping[i] = p->second;
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

/// Revise genotypes in a variant record based on population allele frequencies.
/// With high-coverage sequencing in mind, where we expect the data should
/// usually overwhelm frequency-based priors, we use a light touch aiming for
/// two practical refinements of low-quality/depth genotype calls:
///   1) Low-quality calls of low-quality alleles (that didn't make it through
///      the unifier) can "round down" (ignored- in a principled way!)
///   2) Low-quality homozygous calls of rare alleles (e.g. GT=1/1 DP=2 AD=0,2)
///      shrink to the next most likely heterozygous genotype.
/// Mutates the vr.p pointer.
Status revise_genotypes(const genotyper_config& cfg, const unified_site& us,
                        const map<int, int>& sample_mapping,
                        const bcf_hdr_t* hdr, bcf1_t_plus& vr) {
    assert(!vr.is_ref);
    // Speed optimization: our prior on genotypes will be effectively flat
    // if there are no lost ALT alleles or homozygous-ALT genotypes called, so
    // exit early in that case. Calls of the <NON_REF> symbolic allele count as
    // lost.
    // Strictly speaking, the GQ of the existing called genotype might change
    // slightly if we did continue, but we judge this not to justify the
    // calculations/allocations we'd do to get there.
    bool needs_revision = false;
    for (const auto& sample : sample_mapping) {
        assert(sample.first < vr.p->n_sample);
        auto al1 = vr.gt.v[sample.first*2];
        if (!bcf_gt_is_missing(al1) && vr.allele_mapping.at(bcf_gt_allele(al1)) == -1) {
            needs_revision = true;
            break;
        }
        auto al2 = vr.gt.v[sample.first*2+1];
        if (!bcf_gt_is_missing(al2) && vr.allele_mapping.at(bcf_gt_allele(al2)) == -1) {
            needs_revision = true;
            break;
        }
        if (!bcf_gt_is_missing(al1) && !bcf_gt_is_missing(al2) &&
            bcf_gt_allele(al1) > 0 && bcf_gt_allele(al1) == bcf_gt_allele(al2)) {
            needs_revision = true;
            break;
        }
    }
    if (!needs_revision) {
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

    // construct "prior" over genotypes which penalizes lost ALT alleles and
    // homozygous-ALT genotypes (otherwise flat)
    const float lost_log_prior = log(std::max(us.lost_allele_frequency, cfg.min_assumed_allele_frequency));
    vector<double> gt_log_prior(nGT, 0.0);
    for (unsigned gt = 0; gt < gt_log_prior.size(); gt++) {
        auto als = diploid::gt_alleles(gt);
        if (vr.allele_mapping.at(als.first) == -1 || vr.allele_mapping.at(als.second) == -1) {
            gt_log_prior[gt] = lost_log_prior;
        } else if (als.first > 0 && als.first == als.second) {
            gt_log_prior[gt] = log(std::max(us.alleles[vr.allele_mapping[als.first]].frequency,
                                            cfg.min_assumed_allele_frequency));
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
            auto g_ll = sample_gll[g] + gt_log_prior[g];
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

    // Keep records not overlapping the site only if they actually have an ALT
    // allele that's in the unification (a non-aligned indel). Also, collect
    // ranges
    vector<range> record_rngs;
    vector<shared_ptr<bcf1_t>> relevant_records;
    for (const auto& record: records) {
        assert(bcf_nsamples == record->n_sample);
        range record_rng(record.get());
        bool keep = record_rng.overlaps(site.pos);
        for (int i = 1; !keep && i < record->n_allele; i++) {
            string al(record->d.allele[i]);
            if (is_dna(al)) {
                auto p = site.unification.find(allele(record_rng, al));
                if (p != site.unification.end()) {
                    keep = true;
                }
            }
        }
        if (keep) {
            record_rngs.push_back(record_rng);
            relevant_records.push_back(record);
        }
    }

    // check the records span the site, otherwise we need to produce
    // PartialData no-calls
    if (record_rngs.empty()) {
        return Status::OK(); // MissingData
    }
    if (!cfg.allow_partial_data && !site.pos.spanned_by(record_rngs)) {
        rnc = NoCallReason::PartialData;
        return Status::OK();
    }

    vector<shared_ptr<bcf1_t_plus>> ref_records;
    for (const auto& record : relevant_records) {
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

    // ex post facto check for reference confidence records whose GT is other
    // than 0/0 (probably ./.), which we'll translate to InputNonCalled
    // We exclude 'haploid' records from this treatment for now as observed
    // examples (e.g. in Strelka2 gVCFs) don't seem to require it, but this
    // may need to be configurable in the future.
    for (const auto& rp : all_records) {
        if (rp->is_ref && !rp->was_haploid) {
            for (unsigned i = 0; i < 2*rp->p->n_sample; i++) {
                if (bcf_gt_is_missing(rp->gt[i]) || bcf_gt_allele(rp->gt[i]) != 0) {
                    rnc = NoCallReason::InputNonCalled;
                    return Status::OK();
                }
            }
        }
    }

    // Success...
    rnc = NoCallReason::N_A;
    return Status::OK();
}

/// Based on the cluster of variant records and min_ref_depth produced by
/// prepare_dataset_records, fill genotypes for this dataset's samples with
/// appropriate calls (currently by translation of the input hard-calls).
/// Updates genotypes and may modify min_ref_depth.
///
/// variant_records_used will be filled to the variant records from which
/// the call(s) were actually made, if any.
///
/// FIXME: not coded to deal with multi-sample gVCFs properly.
static Status translate_genotypes(const genotyper_config& cfg, const unified_site& site,
                                  const string& dataset, const bcf_hdr_t* dataset_header,
                                  int bcf_nsamples, const map<int,int>& sample_mapping,
                                  const vector<shared_ptr<bcf1_t_plus>>& variant_records,
                                  AlleleDepthHelper& depth,
                                  vector<int>& min_ref_depth,
                                  vector<one_call>& genotypes,
                                  vector<shared_ptr<bcf1_t_plus>>& variant_records_used) {
    assert(genotypes.size() == 2*min_ref_depth.size());
    assert(!site.monoallelic);
    variant_records_used.clear();
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

    bcf1_t_plus *record = nullptr, *record2 = nullptr;

    int call_mode = -1;  // -1: no-call, record is null
                         //  0: call the 'left' allele of record
                         //  1: call the 'right' allele of record
                         //  2: full diploid call from record; record2 is null

    int call_mode2 = -1; // -1: do not use record2
                         //  0: call the 'left' allele of record2
                         //  1: call the 'right' allele of record2

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
        call_mode = 2;
        variant_records_used.push_back(records_non00[0]);
    } else {
        // complex situation: multiple non-0/0 records overlapping the unified site.
        //
        // If the records don't all share at least one position in common (i.e. their
        // ALT alleles aren't mutually exclusive on one chromosome), punt with
        // UnphasedVariants. We may improve upon this especially if there is phasing
        // info.
        //
        // If there are exactly two non-0/0 records each calling one recognized ALT
        // allele, combine these into one heterozygous call.
        //
        // Otherwise if all the non-0/0 records call one ALT allele, generate a half
        // call from the highest-quality record calling a recognized ALT allele.
        //
        // Otherwise, bail with OverlappingVariants.

        range intersection(records_non00[0]->p);
        for (auto& a_record : records_non00) {
            range record_rng(a_record->p);
            assert(record_rng.rid == intersection.rid);
            intersection.beg = max(record_rng.beg, intersection.beg);
            intersection.end = min(record_rng.end, intersection.end);
        }
        if (intersection.beg >= intersection.end) {
            for (int i = 0; i < bcf_nsamples; i++) {
                genotypes[sample_mapping.at(i)*2].RNC =
                    genotypes[sample_mapping.at(i)*2+1].RNC =
                        NoCallReason::UnphasedVariants;
            }
            return Status::OK();
        }

        vector<tuple<float,shared_ptr<bcf1_t_plus>,bool>> usable_half_calls;
        for (auto& a_record : records_non00) {
            range record_rng(a_record->p);
            assert(record_rng.rid == intersection.rid);
            intersection.beg = max(record_rng.beg, intersection.beg);
            intersection.end = min(record_rng.end, intersection.end);

            assert(a_record->gt.capacity >= 2);
            int al0 = bcf_gt_is_missing(a_record->gt[0]) ? -1 : bcf_gt_allele(a_record->gt[0]),
                al1 = bcf_gt_is_missing(a_record->gt[1]) ? -1 : bcf_gt_allele(a_record->gt[1]);
            if (al0 > 0 && al1 > 0) {
                // record calls multiple ALT alleles; we're going to bail with
                // OverlappingVariants
                usable_half_calls.clear();
                break;
            } else {
                bool usable0 = al0 > 0 && a_record->allele_mapping[al0] > 0;
                bool usable1 = al1 > 0 && a_record->allele_mapping[al1] > 0;
                if (usable0) {
                    assert(!usable1);
                    usable_half_calls.push_back(make_tuple(-1.0f*a_record->p->qual,a_record,false));
                }
                if (usable1) {
                    assert(!usable0);
                    usable_half_calls.push_back(make_tuple(-1.0f*a_record->p->qual,a_record,true));
                }
            }
        }

        if (usable_half_calls.empty()) {
            for (int i = 0; i < bcf_nsamples; i++) {
                genotypes[sample_mapping.at(i)*2].RNC =
                    genotypes[sample_mapping.at(i)*2+1].RNC =
                        NoCallReason::OverlappingVariants;
            }
            return Status::OK();
        }

        // sorts usable half-call records by decreasing qual
        sort(usable_half_calls.begin(), usable_half_calls.end());

        auto& r1 = get<1>(usable_half_calls[0]);
        record = r1.get();
        variant_records_used.push_back(r1);
        call_mode = get<2>(usable_half_calls[0]) ? 1 : 0;
        if (usable_half_calls.size() == 2) {
            auto& r2 = get<1>(usable_half_calls[1]);
            record2 = r2.get();
            variant_records_used.push_back(r2);
            call_mode2 = get<2>(usable_half_calls[1]) ? 1 : 0;
        }
    }

    // Now, translating genotypes from one variant BCF record.
    assert(record != nullptr);
    assert(call_mode >= 0);
    assert(call_mode2 == -1 || record2 != nullptr);
    assert(record->p->n_sample == bcf_hdr_nsamples(dataset_header));

    S(depth.Load(dataset, dataset_header, record->p.get()));

    // for each shared sample, record the genotype call.
    for (const auto& ij : sample_mapping) {
        assert(ij.first < record->p->n_sample);
        assert(ij.second < min_ref_depth.size());

        // TODO: are depth and allele_mapping checks inside-out????
        #define fill_allele(rec,depth,in_ofs,out_ofs)                             \
            assert(rec);                                                          \
            genotypes[2*ij.second+(out_ofs)].RNC = NoCallReason::InputNonCalled;  \
            if (rec->gt[2*ij.first+in_ofs] != bcf_int32_vector_end &&             \
                !bcf_gt_is_missing(rec->gt[2*ij.first+(in_ofs)])) {               \
                auto al = bcf_gt_allele(rec->gt[2*ij.first+(in_ofs)]);            \
                assert(al >= 0 && al < rec->p->n_allele);                         \
                int rd = min_ref_depth[ij.second];                                \
                if (depth.get(ij.first, al) >= cfg.required_dp                    \
                    && (rd < 0 || rd >= cfg.required_dp)) {                       \
                    if (rec->allele_mapping[al] >= 0) {                           \
                        genotypes[2*ij.second+(out_ofs)] =                        \
                            one_call(bcf_gt_unphased(rec->allele_mapping[al]),    \
                                     NoCallReason::N_A);                          \
                    } else {                                                      \
                        genotypes[2*ij.second+(out_ofs)].RNC =                    \
                            rec->deletion_allele[al]                              \
                                ? NoCallReason::LostDeletion                      \
                                : NoCallReason::LostAllele;                       \
                    }                                                             \
                } else {                                                          \
                    genotypes[2*ij.second+(out_ofs)].RNC =                        \
                        NoCallReason::InsufficientDepth;                          \
                }                                                                 \
            }


        switch (call_mode) {
            case 2:
                fill_allele(record,depth,0,0)
                fill_allele(record,depth,1,1)
                break;
            case 0:
            case 1:
                fill_allele(record,depth,call_mode,0)
                genotypes[2*ij.second].half_call = true;
                switch (call_mode2) {
                    case -1:
                        genotypes[2*ij.second+1].RNC = NoCallReason::OverlappingVariants;
                        break;
                    case 0:
                    case 1:
                        {
                            auto depth2 = NewAlleleDepthHelper(cfg);
                            S(depth2->Load(dataset, dataset_header, record2->p.get()));

                            fill_allele(record2,(*depth2),call_mode2,1);
                            assert(genotypes[2*ij.second+1].RNC != NoCallReason::MissingData);
                            genotypes[2*ij.second+1].half_call = true;
                        }
                        break;
                    default:
                        assert(false);
                }
                break;
            default:
                assert(false);
        }
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
                                    vector<one_call>& genotypes,
                                    vector<shared_ptr<bcf1_t_plus>>& variant_records_used) {
    assert(genotypes.size() == 2*min_ref_depth.size());
    assert(site.monoallelic);
    assert(site.alleles.size() == 2);
    variant_records_used.clear();
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
                    variant_records_used.push_back(a_record);
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

    assert(record->p->n_sample == bcf_hdr_nsamples(dataset_header));

    S(depth.Load(dataset, dataset_header, record->p.get()));

    // for each shared sample, record the genotype call.
    for (const auto& ij : sample_mapping) {
        assert(ij.first < record->p->n_sample);
        assert(ij.second < min_ref_depth.size());

        #define fill_monoallelic(ofs)                                             \
            genotypes[2*ij.second+(ofs)].RNC = NoCallReason::InputNonCalled;      \
            if (record->gt[2*ij.first+ofs] != bcf_int32_vector_end &&             \
                !bcf_gt_is_missing(record->gt[2*ij.first+(ofs)])) {               \
                auto al = bcf_gt_allele(record->gt[2*ij.first+(ofs)]);            \
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
    S(setup_format_helpers(format_helpers, cfg, site, samples));

    // query database for pertinent records across the samples -- the range
    // encompassing all the original alleles
    range query_range(site.pos);
    for (const auto& p : site.unification) {
        const range& pr = p.first.pos;
        assert(pr.rid == query_range.rid);
        query_range.beg = min(query_range.beg, pr.beg);
        query_range.end = max(query_range.end, pr.end);
    }
    shared_ptr<const set<string>> samples2, datasets;
    vector<unique_ptr<RangeBCFIterator>> iterators;
    S(data.sampleset_range(cache, sampleset, query_range, nullptr,
                           samples2, datasets, iterators));
    assert(samples.size() == samples2->size());

    auto adh = NewAlleleDepthHelper(cfg);
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
        if (sample_mapping.empty()) {
            continue;
        }

        // pre-process the records
        vector<int> min_ref_depth(samples.size(), -1);
        vector<shared_ptr<bcf1_t_plus>> all_records, variant_records, variant_records_used;
        NoCallReason rnc = NoCallReason::MissingData;
        S(prepare_dataset_records(cfg, site, dataset, dataset_header.get(), bcf_nsamples,
                                  sample_mapping, records, *adh, rnc, min_ref_depth,
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
                                  sample_mapping, variant_records, *adh, min_ref_depth,
                                  genotypes, variant_records_used));
        } else {
            S(translate_monoallelic(cfg, site, dataset, dataset_header.get(), bcf_nsamples,
                                    sample_mapping, variant_records, *adh, min_ref_depth,
                                    genotypes, variant_records_used));
        }

        // Update FORMAT fields for this dataset.
        if (!(cfg.squeeze && variant_records.empty() && !all_records.empty())) {
            S(update_format_fields(cfg, dataset, dataset_header.get(), sample_mapping, site,
                                format_helpers, all_records, variant_records_used));
            // But if rnc = MissingData, PartialData, UnphasedVariants, or OverlappingVariants, then
            // we must censor the FORMAT fields as potentially unreliable/misleading.
            for (const auto& p : sample_mapping) {
                auto rnc1 = genotypes[p.second*2].RNC;
                auto rnc2 = genotypes[p.second*2+1].RNC;
                bool half_call = site.monoallelic || genotypes[p.second*2].half_call || genotypes[p.second*2+1].half_call;

                if (rnc1 == NoCallReason::MissingData || rnc1 == NoCallReason::PartialData) {
                    assert(rnc1 == rnc2);
                    for (const auto& fh : format_helpers) {
                        S(fh->censor(p.second, false));
                    }
                } else if (rnc1 == NoCallReason::UnphasedVariants || rnc2 == NoCallReason::UnphasedVariants ||
                        rnc1 == NoCallReason::OverlappingVariants || rnc2 == NoCallReason::OverlappingVariants) {
                    for (const auto& fh : format_helpers) {
                        if (fh->field_info.name != "DP" && fh->field_info.name != "FT") { // whitelist
                            S(fh->censor(p.second, half_call));
                        }
                    }
                } else if (half_call) {
                    for (const auto& fh : format_helpers) {
                        if (fh->field_info.name != "DP" && fh->field_info.name != "GQ"
                            && fh->field_info.name != "FT") {
                            S(fh->censor(p.second, true));
                        }
                    }
                }
            }
        } else {
            // Short path if cfg.squeeze && variant_records.empty() && !all_records.empty():
            //   Update DP only and apply squeeze transform
            S(update_format_fields(cfg, dataset, dataset_header.get(), sample_mapping, site,
                                   format_helpers, all_records, variant_records_used, true));
            for (const auto& p : sample_mapping) {
                genotypes[p.second*2].RNC = NoCallReason::N_A;
                genotypes[p.second*2+1].RNC = NoCallReason::N_A;
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
        if ((genotypes[2*i].allele != bcf_gt_missing && genotypes[2*i+1].allele == bcf_gt_missing) ||
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
        c_alleles.push_back(allele.dna.c_str());
    }
    if (bcf_update_alleles(hdr, ans.get(), c_alleles.data(), c_alleles.size()) != 0) {
        return Status::Failure("bcf_update_alleles");
    }

    // populate ID column with a normalized representation of each ALT
    ostringstream anr;
    for (int i = 1; i < site.alleles.size(); i++) {
        const auto& norm = site.alleles[i].normalized;
        if (!site.pos.contains(norm.pos)) {
            return Status::Failure("logic error: unified allele normalized representation isn't contained within site", site.pos.str());
        }
        if (i > 1) {
            anr << ";";
        }
        anr << bcf_seqname(hdr, ans.get())
            << "_" << (norm.pos.beg+1)
            << "_" << site.alleles[0].dna.substr(norm.pos.beg - site.pos.beg, norm.pos.size())
            << "_" << norm.dna;
    }
    if (bcf_update_id(hdr, ans.get(), anr.str().c_str())) {
        return Status::Failure("bcf_update_id", anr.str());
    }

    // AF
    vector<float> af;
    bool output_af = true;
    for (int i = 1; i < site.alleles.size(); i++) {
        auto f = site.alleles[i].frequency;
        if (f == f) {
            af.push_back(f);
        } else {
            output_af = false;
            break;
        }
    }
    if (output_af && bcf_update_info_float(hdr, ans.get(), "AF", af.data(), af.size()) != 0) {
        return Status::Failure("bcf_update_info_int32 AQ");
    }

    // AQ
    vector<int32_t> aq;
    bool any_aq = false;
    for (int i = 1; i < site.alleles.size(); i++) {
        auto q = site.alleles[i].quality;
        aq.push_back(q);
        if (q) {
            any_aq = true;
        }
    }
    if (any_aq && bcf_update_info_int32(hdr, ans.get(), "AQ", aq.data(), aq.size()) != 0) {
        return Status::Failure("bcf_update_info_int32 AQ");
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
            RNC_CASE(MonoallelicSite,"1")
            RNC_CASE(InputNonCalled, "I")
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

    if (cfg.trim_uncalled_alleles) {
        if (bcf_trim_alleles(hdr, ans.get()) < 0) {
            return Status::Failure("bcf_trim_alleles");
        }
        if (ans->n_allele < 2) {
            ans.reset();
        }
    }

    // Overwrite the output BCF record with a duplicate. Why? This forces htslib to
    // perform some internal serialization of the data (see the static bcf1_sync
    // function in vcf.c, which we can't call directly, but is called by bcf_dup).
    // htslib would otherwise do this serialization implicitly while writing the
    // record out to a file, but by doing it explicitly here, we get to do some of the
    // work in the current worker thread rather than the single thread responsible for
    // writing out the file.
    if (ans) {
        auto ans2 = shared_ptr<bcf1_t>(bcf_dup(ans.get()), &bcf_destroy);
        ans = move(ans2);
    }

    return Status::OK();
}

}

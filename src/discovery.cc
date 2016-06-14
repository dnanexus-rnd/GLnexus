#include "discovery.h"
#include "diploid.h"

using namespace std;

namespace GLnexus {

// for each allele, find the max GQ of any hard-call of the allele across the given sample indices
Status alleles_maxGQ(const bcf_hdr_t* hdr, bcf1_t* record, const vector<unsigned> samples, vector<int>& ans) {
    auto dataset_nsamples = bcf_hdr_nsamples(hdr);
    ans.assign(record->n_allele,0);

    htsvecbox<int> gt;
    htsvecbox<int32_t> gq;
    int nGT = bcf_get_genotypes(hdr, record, &gt.v, &gt.capacity);
    if (!gt.v || nGT != 2*dataset_nsamples) return Status::Failure("discover_alleles bcf_get_genotypes");
    int nGQ = bcf_get_format_int32(hdr, record, "GQ", &gq.v, &gq.capacity);

    if (nGQ != dataset_nsamples) return Status::OK();

    for (unsigned sample : samples) {
        assert(sample < dataset_nsamples);
        #define update_maxGQ(ofs) \
            if (!bcf_gt_is_missing(gt[2*sample+(ofs)])) {  \
                int al = bcf_gt_allele(gt[2*sample+(ofs)]); \
                if (al < 0 || al >= record->n_allele) return Status::Failure("discover_alleles bcf_get_genotypes invalid allele"); \
                ans[al] = std::max(ans[al], gq[sample]); \
            }
        update_maxGQ(0)
        update_maxGQ(1)
    }

    return Status::OK();
}

Status discover_alleles_from_iterator(const set<string>& samples,
                                      const range& pos,
                                      RangeBCFIterator& iterator,
                                      discovered_alleles& final_dsals) {
    Status s;

    // get dataset BCF records
    string dataset;
    shared_ptr<const bcf_hdr_t> dataset_header;
    vector<shared_ptr<bcf1_t>> records;
    vector<vector<double>> copy_number;
    while ((s = iterator.next(dataset, dataset_header, records)).ok()) {
        discovered_alleles dsals;
        // determine which of the dataset's samples are in the desired sample set
        // TODO: FCMM-memoize this during the operation
        size_t dataset_nsamples = bcf_hdr_nsamples(dataset_header.get());
        vector<unsigned> dataset_relevant_samples;
        for (unsigned i = 0; i < dataset_nsamples; i++) {
            if (samples.find(string(bcf_hdr_int2id(dataset_header.get(), BCF_DT_SAMPLE, i))) != samples.end()) {
                dataset_relevant_samples.push_back(i);
            }
        }

        // for each BCF record
        for (const auto& record : records) {
            assert(!is_gvcf_ref_record(record.get()));
            range rng(record);
            assert(pos.overlaps(rng));
            if (!pos.contains(rng)) {
                // Skip records that dangle off the edge of the target range.
                // The problem is that such alleles could overlap other
                // alleles not overlapping the target range at all, and that
                // we therefore won't see in our query. So it wouldn't be safe
                // for us to make claims about the genotypes of the resulting
                // sites.
                continue;
            }

            // find the max GQ for each allele
            vector<int> maxGQ;
            S(alleles_maxGQ(dataset_header.get(), record.get(), dataset_relevant_samples, maxGQ));

            // calculate estimated allele copy numbers for this record
            // TODO: configurable bias00
            // TODO: ideally we'd compute them only for relevant samples
            S(diploid::estimate_allele_copy_number(dataset_header.get(), record.get(), 1.0, copy_number));
            #define round4(x) (float(round(x*10000.0)/10000.0))

            // FIXME -- minor potential bug -- double-counting copy number of
            // alleles that span multiple discovery ranges

            // create a discovered_alleles entry for each alt allele matching [ACGT]+
            // In particular this excludes gVCF <NON_REF> symbolic alleles
            bool any_alt = false;
            for (int i = 1; i < record->n_allele; i++) {
                double copy_number_i = 0.0;
                for (unsigned sample : dataset_relevant_samples) {
                    copy_number_i += copy_number[sample][i];
                }
                if (copy_number_i) {
                    string aldna(record->d.allele[i]);
                    transform(aldna.begin(), aldna.end(), aldna.begin(), ::toupper);
                    if (aldna.size() > 0 && regex_match(aldna, regex_dna)) {
                        discovered_allele_info ai = { false, round4(copy_number_i), maxGQ[i] };
                        dsals.insert(make_pair(allele(rng, aldna), ai));
                        any_alt = true;
                    }
                }
            }

            // create an entry for the ref allele, if we discovered at least one alt allele.
            string refdna(record->d.allele[0]);
            transform(refdna.begin(), refdna.end(), refdna.begin(), ::toupper);
            if (refdna.size() > 0 && regex_match(refdna, regex_dna)) {
                if (any_alt) {
                    double ref_copy_number = 0.0;
                    for (unsigned sample : dataset_relevant_samples) {
                        ref_copy_number += copy_number[sample][0];
                    }
                    discovered_allele_info ai = { true, round4(ref_copy_number), maxGQ[0] };
                    dsals.insert(make_pair(allele(rng, refdna), ai));
                }
            } else {
                ostringstream errmsg;
                errmsg << dataset << " " << refdna << "@" << rng.str();
                return Status::Invalid("invalid reference allele", errmsg.str());
            }
        }
        S(merge_discovered_alleles(dsals, final_dsals));
        records.clear();
    }

    if (s != StatusCode::NOT_FOUND) {
        return s;
    }
    return Status::OK();
}

Status discovered_alleles_refcheck(const discovered_alleles& als,
                                   const vector<pair<string,size_t>>& contigs) {
    set<range> ranges;
    multimap<range,pair<string,bool> > refcheck;
    for (const auto& it : als) {
        UNPAIR(it,allele,ai)
        refcheck.insert(make_pair(allele.pos, make_pair(allele.dna, ai.is_ref)));
        ranges.insert(allele.pos);
    }
    for (const auto& rng : ranges) {
        vector<string> refs;
        auto its = refcheck.equal_range(rng);
        for (auto it = its.first; it != its.second; it++) {
            UNPAIR(*it, refcheck_rng, refcheckp)
            UNPAIR(refcheckp, refcheck_dna, refcheck_is_ref)
            assert (refcheck_rng == rng);
            if (refcheck_is_ref) {
                refs.push_back(refcheck_dna);
            }
        }
        if (refs.size() > 1) {
            ostringstream errmsg;
            errmsg << rng.str(contigs);
            for (const auto& r : refs) {
                errmsg << ' ' << r;
            }
            return Status::Invalid("data sets contain inconsistent reference alleles", errmsg.str());
        } else if (refs.size() == 0) {
            return Status::Invalid("data sets contain no reference allele", rng.str(contigs));
        }
    }

    return Status::OK();
}

}

#include "discovery.h"
#include "diploid.h"

using namespace std;

namespace GLnexus {

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
        vector<top_AQ> topAQ;
        vector<zygosity_by_GQ> zGQ;
        for (const auto& record : records) {
            assert(!is_gvcf_ref_record(record.get()));
            range rng(record);
            assert(pos.overlaps(rng));
            bool filtered = (bcf_has_filter(dataset_header.get(), record.get(), ".") == 0);

            // find the max AQ for each allele
            S(diploid::bcf_alleles_topAQ(dataset_header.get(), record.get(), dataset_relevant_samples, topAQ));

            // find zygosity_by_GQ for each allele
            S(diploid::bcf_zygosity_by_GQ(dataset_header.get(), record.get(), dataset_relevant_samples, zGQ));

            // FIXME -- minor potential bug -- double-counting copy number of
            // alleles that span multiple discovery ranges

            // create a discovered_alleles entry for each alt allele matching [ACGT]+
            // In particular this excludes gVCF <NON_REF> symbolic alleles, and any
            // ALT alleles containing IUPAC degenerate letters. Also, we want to see
            // at least one copy of the allele called (no GQ threshold)
            bool any_alt = false;
            for (int i = 1; i < record->n_allele; i++) {
                string aldna(record->d.allele[i]);
                transform(aldna.begin(), aldna.end(), aldna.begin(), ::toupper);
                if (aldna.size() > 0 && is_dna(aldna) && zGQ[i].copy_number(0) > 0) {
                    discovered_allele_info ai;
                    ai.is_ref = false;
                    ai.all_filtered = filtered;
                    ai.topAQ = topAQ[i];
                    ai.zGQ = zGQ[i];
                    dsals.insert(make_pair(allele(rng, aldna), ai));
                    any_alt = true;
                }
            }

            // create an entry for the ref allele, if we discovered at least one alt allele.
            string refdna(record->d.allele[0]);
            transform(refdna.begin(), refdna.end(), refdna.begin(), ::toupper);
            if (refdna.size() > 0 && is_iupac_nucleotides(refdna)) {
                if (any_alt) {
                    discovered_allele_info ai;
                    ai.is_ref = true;
                    ai.all_filtered = filtered;
                    ai.topAQ = topAQ[0];
                    ai.zGQ = zGQ[0];
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

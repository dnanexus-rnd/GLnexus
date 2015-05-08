#include "service.h"
#include "data.h"
#include "alleles.h"
#include <algorithm>
#include <map>
#include <assert.h>

using namespace std;

namespace GLnexus{

Status Service::Start(Data *data, unique_ptr<Service>& svc) {
    svc.reset(new Service());
    svc->data_ = data;
    // TODO: contigs_ should be loaded from service config, and checked against Data?
    return data->contigs(svc->contigs_);
}

Status Service::discover_alleles(const string& sampleset, const range& pos, discovered_alleles& ans) {
    Status s;

    // Find the data sets containing the samples in the sample set.
    // TODO results could be cached
    shared_ptr<const set<string> > samples;
    set<string> datasets;
    s = data_->sampleset_samples(sampleset, samples);
    if (s.bad()) return s;
    for (const auto& it : *samples) {
        string dataset;
        s = data_->sample_dataset(it, dataset);
        if (s.bad()) return s;
        datasets.insert(dataset);
    }

    // extract alleles from each dataset
    ans.clear();
    for (const auto& dataset : datasets) {
        // get dataset BCF records
        shared_ptr<const bcf_hdr_t> dataset_header;
        vector<shared_ptr<bcf1_t>> records;
        s = data_->dataset_bcf(dataset, pos, dataset_header, records);
        if (s.bad()) return s;

        // for each BCF record
        discovered_alleles dsals;
        for (const auto& record : records) {
            range rng(*record);
            vector<float> obs_counts(record->n_allele, 0.0);

            // count hard-called allele observations
            // TODO: could use GLs for soft estimate
            // TODO: "max ref extension" distance for each allele
            int *gt = nullptr, gtsz = 0;
            int ngt = bcf_get_genotypes(dataset_header.get(), record.get(), &gt, &gtsz);
            for (int i = 0; i < ngt; i++) {
                if (gt[i] != bcf_int32_vector_end) {
                    int al_i = bcf_gt_allele(gt[i]);
                    if (al_i >= 0 && al_i < record->n_allele) {
                        obs_counts[al_i] += 1.0;
                    }
                }
            }
            if (gt) {
                free(gt);
            }

            // create a discovered_alleles entry
            for (int i = 0; i < record->n_allele; i++) {
                allele al(rng, record->d.allele[i]);
                bool is_ref = (i == 0);
                dsals.insert(make_pair(al, make_pair(is_ref, obs_counts[i])));
            }
        }

        // merge in this dataset's alleles
        s = merge_discovered_alleles(dsals, ans);
        if (s.bad()) return s;
    }

    // ex post facto check: exactly one reference allele at any given range
    // TODO: should be factored out because the distributed service will also need to check this
    set<range> ranges;
    multimap<range,pair<string,bool> > refcheck;
    for (const auto& it : ans) {
        UNPAIR(it,allele,ai)
        UNPAIR(ai,is_ref,obs_count)
        assert(obs_count == obs_count);
        refcheck.insert(make_pair(allele.pos, make_pair(allele.dna, is_ref)));
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
            errmsg << rng.str(contigs_);
            for (const auto& r : refs) {
                errmsg << ' ' << r;
            }
            return Status::Invalid("data sets contain inconsistent reference alleles", errmsg.str());
        } else if (refs.size() == 0) {
            return Status::Invalid("data sets contain no reference allele", rng.str(contigs_));
        }
    }

    return Status::OK();
}

Status Service::genotype_sites(const string& sampleset, const vector<unified_site>& sites, std::string& filename) {
    // for each site,
    // for each dataset,
    // load the dataset's overlapping records
    // load the dataset's GLs for pertinent samples
    // unify the alleles and call the genotypes in a new bcf1_t record.
    return Status::Failure();
}

}

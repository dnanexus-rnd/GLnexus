#include "service.h"
#include "data.h"
#include "alleles.h"
#include "genotyper.h"
#include <algorithm>
#include <map>
#include <assert.h>
#include <sstream>

using namespace std;

namespace GLnexus{

Status Service::Start(Data *data, unique_ptr<Service>& svc) {
    svc.reset(new Service());
    return DataCache::Start(data, svc->data_);
}

bool is_dna(const string& str) {
    return all_of(str.begin(), str.end(),
                  [](char ch) { return ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T'; });
}

Status Service::discover_alleles(const string& sampleset, const range& pos, discovered_alleles& ans) {
    // Find the data sets containing the samples in the sample set.
    shared_ptr<const set<string>> samples, datasets;
    Status s;
    S(data_->sampleset_datasets(sampleset, samples, datasets));

    // extract alleles from each dataset
    ans.clear();
    for (const auto& dataset : *datasets) {
        // get dataset BCF records
        shared_ptr<const bcf_hdr_t> dataset_header;
        vector<shared_ptr<bcf1_t>> records;
        S(data_->dataset_bcf(dataset, pos, dataset_header, records));

        // for each BCF record
        discovered_alleles dsals;
        for (const auto& record : records) {
            range rng(record);
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

            // create a discovered_alleles entry for each alt allele matching [ACGT]+
            // In particular this excludes gVCF <NON_REF> symbolic alleles
            bool any_alt = false;
            for (int i = 1; i < record->n_allele; i++) {
                string aldna(record->d.allele[i]);
                transform(aldna.begin(), aldna.end(), aldna.begin(), ::toupper);
                if (aldna.size() > 0 && is_dna(aldna)) {
                    discovered_allele_info ai = { false, obs_counts[i] };
                    dsals.insert(make_pair(allele(rng, aldna), ai));
                    any_alt = true;
                }
            }

            // create an entry for the ref allele, if we discovered at least one alt allele.
            string refdna(record->d.allele[0]);
            transform(refdna.begin(), refdna.end(), refdna.begin(), ::toupper);
            if (refdna.size() > 0 && is_dna(refdna)) {
                if (any_alt) {
                    discovered_allele_info ai = { true, obs_counts[0] };
                    dsals.insert(make_pair(allele(rng, refdna), ai));
                }
            } else {
                ostringstream errmsg;
                errmsg << dataset << " " << refdna << "@" << pos.str();
                return Status::Invalid("invalid reference allele", errmsg.str());
            }
        }

        // merge in this dataset's alleles
        S(merge_discovered_alleles(dsals, ans));
    }

    // ex post facto check: exactly one reference allele at any given range
    // TODO: should be factored out because the distributed service will also need to check this
    set<range> ranges;
    multimap<range,pair<string,bool> > refcheck;
    for (const auto& it : ans) {
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
            errmsg << rng.str(data_->contigs());
            for (const auto& r : refs) {
                errmsg << ' ' << r;
            }
            return Status::Invalid("data sets contain inconsistent reference alleles", errmsg.str());
        } else if (refs.size() == 0) {
            return Status::Invalid("data sets contain no reference allele", rng.str(data_->contigs()));
        }
    }

    return Status::OK();
}

Status Service::genotype_sites(const genotyper_config& cfg, const string& sampleset, const vector<unified_site>& sites, const string& filename) {
    Status s;
    shared_ptr<const set<string>> samples, datasets;
    S(data_->sampleset_datasets(sampleset, samples, datasets));

    // get a BCF header for this sample set
    // TODO: memoize
    shared_ptr<bcf_hdr_t> hdr(bcf_hdr_init("w"), &bcf_hdr_destroy);
    const char* hdrGT = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
    if (bcf_hdr_append(hdr.get(), hdrGT) != 0) {
        return Status::Failure("bcf_hdr_append", hdrGT);
    }
    for (const auto& ctg : data_->contigs()) {
        ostringstream stm;
        stm << "##contig=<ID=" << ctg.first << ",length=" << ctg.second << ">";
        if (bcf_hdr_append(hdr.get(), stm.str().c_str()) != 0) {
            return Status::Failure("bcf_hdr_append", stm.str());
        }
    }
    for (const auto& sample : *samples) {
        if (bcf_hdr_add_sample(hdr.get(), sample.c_str()) != 0) {
            return Status::Failure("bcf_hdr_add_sample", sample);
        }
    }
    if (bcf_hdr_sync(hdr.get()) != 0) {
        return Status::Failure("bcf_hdr_sync");
    }

    // open output BCF file
    unique_ptr<vcfFile, void(*)(vcfFile*)> outfile(bcf_open(filename.c_str(), "wb"), [](vcfFile* f) { bcf_close(f); });
    if (!outfile) {
        return Status::IOError("failed to open BCF file for writing", filename);
    }
    if (bcf_hdr_write(outfile.get(), hdr.get()) != 0) {
        return Status::IOError("bcf_hdr_write", filename);
    }

    // for each site
    for (const auto& site : sites) {
        // compute genotypes
        shared_ptr<bcf1_t> site_bcf;
        S(genotype_site(cfg, *data_, site, *samples, *datasets, hdr.get(), site_bcf));

        // write out a BCF record
        if (bcf_write(outfile.get(), hdr.get(), site_bcf.get()) != 0) {
            return Status::IOError("bcf_write", filename);
        }
    }

    if (bcf_close(outfile.release()) != 0) {
        return Status::IOError("bcf_close", filename);
    }

    return Status::OK();
}

}

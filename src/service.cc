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

Status Service::sampleset_datasets(const string& sampleset, shared_ptr<const set<string>>& ans) {
    // TODO cache this stuff
    shared_ptr<const set<string> > samples;
    auto datasets = make_shared<set<string>>();
    Status s;
    S(data_->sampleset_samples(sampleset, samples));
    for (const auto& it : *samples) {
        string dataset;
        S(data_->sample_dataset(it, dataset));
        datasets->insert(dataset);
    }
    ans = datasets;
    return Status::OK();
}

Status Service::discover_alleles(const string& sampleset, const range& pos, discovered_alleles& ans) {
    // Find the data sets containing the samples in the sample set.
    shared_ptr<const set<string>> datasets;
    Status s;
    S(sampleset_datasets(sampleset, datasets));

    // extract alleles from each dataset
    ans.clear();
    for (const auto& dataset : *datasets) {
        // read the header
        shared_ptr<const bcf_hdr_t> dataset_header;
        S(data_->dataset_bcf_header(dataset, dataset_header));

        // get dataset BCF records
        vector<shared_ptr<bcf1_t>> records;
        S(data_->dataset_bcf(dataset, pos, records));

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
        S(merge_discovered_alleles(dsals, ans));
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

Status Service::genotype_sites(const string& sampleset, const vector<unified_site>& sites, const string& filename) {
    Status s;
    shared_ptr<const set<string>> samples, datasets;
    S(data_->sampleset_samples(sampleset, samples));
    S(data_->sampleset_datasets(sampleset, datasets));

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
        S(genotype_site(*data_, site, *samples, *datasets, hdr.get(), site_bcf));

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

#include "service.h"
#include "data.h"
#include "unifier.h"
#include "genotyper.h"
#include <algorithm>
#include <iostream>
#include <map>
#include <assert.h>
#include <sstream>
#include "ctpl_stl.h"

using namespace std;

namespace GLnexus{

// pImpl idiom
struct Service::body {
    BCFData& data_;
    std::unique_ptr<MetadataCache> metadata_;
    ctpl::thread_pool threadpool_;

    body(BCFData& data) : data_(data) {}
};

Service::Service(BCFData& data) {
    body_ = make_unique<Service::body>(data);
    body_->threadpool_.resize(thread::hardware_concurrency());
}

Service::~Service() = default;

Status Service::Start(Metadata& metadata, BCFData& data, unique_ptr<Service>& svc) {
    svc.reset(new Service(data));
    return MetadataCache::Start(metadata, svc->body_->metadata_);
}

bool is_dna(const string& str) {
    return all_of(str.begin(), str.end(),
                  [](char ch) { return ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T'; });
}

// discover_alleles in one dataset -- used on the thread pool below
static Status discover_alleles_in_dataset(BCFData& data,
                                          const set<string>& samples,
                                          const string& dataset, const range& pos,
                                          discovered_alleles& dsals) {
    Status s;

    // get dataset BCF records
    shared_ptr<const bcf_hdr_t> dataset_header;
    vector<shared_ptr<bcf1_t>> records;
    S(data.dataset_range_and_header(dataset, pos, dataset_header, records));

    // determine which of the dataset's samples are in the desired sample set
    int dataset_nsamples = bcf_hdr_nsamples(dataset_header.get());
    vector<bool> dataset_sample_relevant;
    for (size_t i = 0; i < dataset_nsamples; i++) {
        string sample_i(bcf_hdr_int2id(dataset_header.get(), BCF_DT_SAMPLE, i));
        dataset_sample_relevant.push_back(samples.find(sample_i) != samples.end());
    }

    // for each BCF record
    dsals.clear();
    for (auto& record : records) {
        range rng(record);
        vector<float> obs_counts(record->n_allele, 0.0);

        // count hard-called alt allele observations for desired samples
        // TODO: could use GLs for soft estimate
        // TODO: "max ref extension" distance for each allele
        int *gt = nullptr, gtsz = 0;
        int ngt = bcf_get_genotypes(dataset_header.get(), record.get(), &gt, &gtsz);
        assert(ngt == 2*dataset_nsamples);
        for (int i = 0; i < ngt; i++) {
            if (gt[i] != bcf_int32_vector_end) {
                int al_i = bcf_gt_allele(gt[i]);
                if (al_i >= 0 && al_i < record->n_allele
                    && dataset_sample_relevant.at(i/2)) {
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
            if (obs_counts[i] > 0.0) { // TODO: threshold for soft estimates
                string aldna(record->d.allele[i]);
                transform(aldna.begin(), aldna.end(), aldna.begin(), ::toupper);
                if (aldna.size() > 0 && is_dna(aldna)) {
                    discovered_allele_info ai = { false, obs_counts[i] };
                    dsals.insert(make_pair(allele(rng, aldna), ai));
                    any_alt = true;
                }
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

    return Status::OK();
}

// verify that the discovered_alleles encodes the same reference allele at
// each distinct range
static Status discovered_alleles_refcheck(const discovered_alleles& als,
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

Status Service::discover_alleles(const string& sampleset, const range& pos, discovered_alleles& ans) {
    // Find the data sets containing the samples in the sample set.
    shared_ptr<const set<string>> samples, datasets;
    Status s;
    S(body_->metadata_->sampleset_datasets(sampleset, samples, datasets));

    // Enqueue processing of each dataset on the thread pool.
    // TODO: improve cache-friendliness for long ranges
    atomic<bool> abort(false);
    vector<future<Status>> statuses;
    vector<discovered_alleles> results(datasets->size());
    // ^^^ results to be filled by side-effect in the individual tasks below.
    // We assume that by virtue of preallocating, no mutex is necessary to
    // use it as follows because writes and reads of individual elements are
    // serialized by the futures.
    size_t i = 0;
    for (const auto& dataset : *datasets) {
        auto fut = body_->threadpool_.push([&, i, dataset](int tid){
            if (abort) {
                return Status::Invalid();
            }

            discovered_alleles dsals;
            Status ls = discover_alleles_in_dataset(body_->data_, *samples, dataset, pos, dsals);
            results[i] = move(dsals);
            return ls;
        });

        statuses.push_back(move(fut));
        i++;
    }
    assert(statuses.size() == datasets->size());

    // Retrieve the results and merge them into ans. Record the first error
    // that occurs, if any, but always wait for all tasks to finish.
    ans.clear();
    s = Status::OK();
    for (size_t i = 0; i < datasets->size(); i++) {
        // wait for task i to complete and find out its status
        Status s_i(statuses[i].get());
        discovered_alleles dsals = move(results[i]);

        if (s.ok() && s_i.ok()) {
            s = merge_discovered_alleles(dsals, ans);
        } else if (s.ok() && s_i.bad()) {
            // record the first error, and tell remaining tasks to abort
            s = move(s_i);
            abort = true;
        }
    }
    if (s.bad()) {
        return s;
    }

    return discovered_alleles_refcheck(ans, body_->metadata_->contigs());
}

Status Service::genotype_sites(const genotyper_config& cfg, const string& sampleset, const vector<unified_site>& sites, const string& filename, consolidated_loss& dlosses) {
    Status s;
    shared_ptr<const set<string>> samples, datasets;
    S(body_->metadata_->sampleset_datasets(sampleset, samples, datasets));
    vector<string> sample_names(samples->begin(), samples->end());

    // create a BCF header for this sample set
    // TODO: memoize
    vector<string> hdr_lines;
    hdr_lines.push_back("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    hdr_lines.push_back("##FORMAT=<ID=RNC,Number=G,Type=Character,Description=\"Reason for No Call in GT: . = N/A, M = Missing data, L = Lost/unrepresentable allele\">");
    for (const auto& ctg : body_->metadata_->contigs()) {
        ostringstream stm;
        stm << "##contig=<ID=" << ctg.first << ",length=" << ctg.second << ">";
        hdr_lines.push_back(stm.str());
    }

    shared_ptr<bcf_hdr_t> hdr(bcf_hdr_init("w"), &bcf_hdr_destroy);
    for (const auto& line : hdr_lines) {
        if (bcf_hdr_append(hdr.get(), line.c_str()) != 0) {
            return Status::Failure("bcf_hdr_append", line);
        }
    }
    for (const auto& sample : sample_names) {
        if (bcf_hdr_add_sample(hdr.get(), sample.c_str()) != 0) {
            return Status::Failure("bcf_hdr_add_sample", sample);
        }
    }
    if (bcf_hdr_sync(hdr.get()) != 0) {
        return Status::Failure("bcf_hdr_sync");
    }

    // safeguard against update_genotypes failure from improper header
    assert(bcf_hdr_nsamples(hdr) == samples->size());

    // open output BCF file
    unique_ptr<vcfFile, void(*)(vcfFile*)> outfile(bcf_open(filename.c_str(), "wb"), [](vcfFile* f) { bcf_close(f); });
    if (!outfile) {
        return Status::IOError("failed to open BCF file for writing", filename);
    }
    if (bcf_hdr_write(outfile.get(), hdr.get()) != 0) {
        return Status::IOError("bcf_hdr_write", filename);
    }

    // Enqueue processing of each site as a task on the thread pool.
    vector<future<Status>> statuses;
    vector<pair<shared_ptr<bcf1_t>,consolidated_loss>> results(sites.size());
    // ^^^ results to be filled by side-effect in the individual tasks below.
    // We assume that by virtue of preallocating, no mutex is necessary to
    // use it as follows because writes and reads of individual elements are
    // serialized by the futures.
    atomic<bool> abort(false);
    for (size_t i = 0; i < sites.size(); i++) {
        auto fut = body_->threadpool_.push([&, i](int tid){
            if (abort) {
                return Status::Invalid();
            }
            shared_ptr<bcf1_t> bcf;
            consolidated_loss losses_for_site;
            Status ls = genotype_site(cfg, body_->data_, sites[i], sample_names, *datasets, hdr.get(), bcf, losses_for_site);
            if (ls.bad()) {
                return ls;
            }
            results[i] = make_pair(move(bcf),losses_for_site);
            return ls;
        });
        statuses.push_back(move(fut));
    }
    assert(statuses.size() == sites.size());

    // Retrieve the resulting BCF records, and write them to the output file,
    // in the given order. Record the first error that occurs, if any, but
    // always wait for all tasks to finish.
    s = Status::OK();
    for (size_t i = 0; i < sites.size(); i++) {
        // wait for task i to complete and find out its status
        Status s_i(statuses[i].get());
        // always retrieve the result BCF record, if any, to ensure we'll free
        // the memory it takes ASAP
        shared_ptr<bcf1_t> bcf_i = move(results[i].first);
        assert(!results[i].first);
        consolidated_loss losses_for_site = results[i].second;

        if (s.ok() && s_i.ok()) {
            // if everything's OK, proceed to write the record
            assert(bcf_i);
            if (bcf_write(outfile.get(), hdr.get(), bcf_i.get()) != 0) {
                s = Status::IOError("bcf_write", filename);
            }
            merge_loss_stats(losses_for_site, dlosses);
        } else if (s.ok() && s_i.bad()) {
            // record the first error, and tell remaining tasks to abort
            s = move(s_i);
            abort = true;
        }
    }
    if (s.bad()) {
        return s;
    }
    // TODO: for very large sample sets, bucket cache-friendliness might be
    // improved by genotyping in grid squares of N>1 sites and M>1 samples

    // close the output file
    return bcf_close(outfile.release()) == 0
                ? Status::OK() : Status::IOError("bcf_close", filename);
}

}

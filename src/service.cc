#include "service.h"
#include "data.h"
#include "unifier.h"
#include "genotyper.h"
#include "residuals.h"
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iostream>
#include <map>
#include <assert.h>
#include "ctpl_stl.h"

using namespace std;

namespace GLnexus{

// pImpl idiom
struct Service::body {
    BCFData& data_;
    std::unique_ptr<MetadataCache> metadata_;

    // thread pool for executing discover_alleles and genotype_sites operations
    ctpl::thread_pool threadpool_;

    // "meta" thread pool for driving multiple concurrent discover_alleles and
    // "genotype_sites operations
    ctpl::thread_pool metapool_;

    body(BCFData& data) : data_(data) {}
};

Service::Service(BCFData& data) {
    body_ = make_unique<Service::body>(data);
    body_->threadpool_.resize(thread::hardware_concurrency());
    body_->metapool_.resize(thread::hardware_concurrency());
}

Service::~Service() = default;

Status Service::Start(Metadata& metadata, BCFData& data, unique_ptr<Service>& svc) {
    svc.reset(new Service(data));
    return MetadataCache::Start(metadata, svc->body_->metadata_);
}

// used on the thread pool below
static Status discover_alleles_thread(const set<string>& samples,
                                      const range& pos,
                                      RangeBCFIterator& iterator,
                                      discovered_alleles& final_dsals) {
    Status s;

    // get dataset BCF records
    string dataset;
    shared_ptr<const bcf_hdr_t> dataset_header;
    vector<shared_ptr<bcf1_t>> records;
    while ((s = iterator.next(dataset, dataset_header, records)).ok()) {
        discovered_alleles dsals;
        // determine which of the dataset's samples are in the desired sample set
        // TODO: FCMM-memoize this during the operation
        int dataset_nsamples = bcf_hdr_nsamples(dataset_header.get());
        vector<bool> dataset_sample_relevant;
        for (size_t i = 0; i < dataset_nsamples; i++) {
            string sample_i(bcf_hdr_int2id(dataset_header.get(), BCF_DT_SAMPLE, i));
            dataset_sample_relevant.push_back(samples.find(sample_i) != samples.end());
        }

        // for each BCF record
        for (const auto& record : records) {
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

            // FIXME -- minor potential bug -- double-counting observations of
            // alleles that span multiple discovery ranges

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
                errmsg << dataset << " " << refdna << "@" << rng.str();
                return Status::Invalid("invalid reference allele", errmsg.str());
            }
        }
        S(merge_discovered_alleles(dsals, final_dsals));
    }

    if (s != StatusCode::NOT_FOUND) {
        return s;
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
    vector<unique_ptr<RangeBCFIterator>> iterators;
    Status s;
    S(body_->data_.sampleset_range(*(body_->metadata_), sampleset, pos,
                                   samples, datasets, iterators));

    // Enqueue processing of each dataset on the thread pool.
    // TODO: improve cache-friendliness for long ranges
    atomic<bool> abort(false);
    vector<future<Status>> statuses;
    vector<discovered_alleles> results(iterators.size());
    // ^^^ results to be filled by side-effect in the individual tasks below.
    // We assume that by virtue of preallocating, no mutex is necessary to
    // use it as follows because writes and reads of individual elements are
    // serialized by the futures.
    size_t i = 0;
    for (const auto& iterator : iterators) {
        RangeBCFIterator* raw_iter = iterator.get();
        auto fut = body_->threadpool_.push([&, i, raw_iter](int tid){
            if (abort) {
                return Status::Invalid();
            }

            discovered_alleles dsals;
            Status ls = discover_alleles_thread(*samples, pos, *raw_iter, dsals);
            results[i] = move(dsals);
            return ls;
        });

        statuses.push_back(move(fut));
        i++;
    }
    assert(statuses.size() == iterators.size());

    // Retrieve the results and merge them into ans. Record the first error
    // that occurs, if any, but always wait for all tasks to finish.
    ans.clear();
    s = Status::OK();
    for (size_t i = 0; i < iterators.size(); i++) {
        // wait for task i to complete and find out its status
        Status s_i(statuses[i].get());
        discovered_alleles dsals = move(results[i]);

        if (s.ok() && s_i.ok()) {
            s = merge_discovered_alleles(dsals, ans);
            if (s.bad()) {
                abort = true;
            }
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

Status Service::discover_alleles(const string& sampleset, const vector<range>& ranges,
                                 vector<discovered_alleles>& ans) {
    atomic<bool> abort(false);
    vector<future<Status>> statuses;
    vector<discovered_alleles> results(ranges.size());

    size_t i = 0;
    for (const auto& range : ranges) {
        auto fut = body_->metapool_.push([&, i, range](int tid){
            if (abort) {
                return Status::Invalid();
            }

            discovered_alleles dsals;
            Status ls = discover_alleles(sampleset, range, dsals);
            results[i] = move(dsals);
            return ls;
        });

        statuses.push_back(move(fut));
        i++;
    }

    ans.clear();
    Status s = Status::OK();
    for (i = 0; i < ranges.size(); i++) {
        // wait for task i to complete and find out its status
        Status s_i(statuses[i].get());
        discovered_alleles dsals = move(results[i]);

        if (s.ok() && s_i.ok()) {
            ans.push_back(move(dsals));
        } else if (s.ok() && s_i.bad()) {
            // record the first error, and tell remaining tasks to abort
            s = move(s_i);
            abort = true;
        }
    }
    assert(s.bad() || ans.size() == ranges.size());

    return s;
}

static Status prepare_bcf_header(const vector<pair<string,size_t> >& contigs,
                                 const vector<string>& samples,
                                 shared_ptr<bcf_hdr_t>& ans) {
    vector<string> hdr_lines;
    hdr_lines.push_back("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    hdr_lines.push_back("##FORMAT=<ID=RNC,Number=G,Type=Character,Description=\"Reason for No Call in GT: . = N/A, M = Missing data, L = Lost/unrepresentable allele\">");
    for (const auto& ctg : contigs) {
        ostringstream stm;
        stm << "##contig=<ID=" << ctg.first << ",length=" << ctg.second << ">";
        hdr_lines.push_back(stm.str());
    }

    shared_ptr<bcf_hdr_t> hdr(bcf_hdr_init("w1"), &bcf_hdr_destroy);
    for (const auto& line : hdr_lines) {
        if (bcf_hdr_append(hdr.get(), line.c_str()) != 0) {
            return Status::Failure("bcf_hdr_append", line);
        }
    }
    for (const auto& sample : samples) {
        if (bcf_hdr_add_sample(hdr.get(), sample.c_str()) != 0) {
            return Status::Failure("bcf_hdr_add_sample", sample);
        }
    }
    if (bcf_hdr_sync(hdr.get()) != 0) {
        return Status::Failure("bcf_hdr_sync");
    }

     // safeguard against update_genotypes failure from improper header
    assert(bcf_hdr_nsamples(hdr) == samples.size());
    ans = hdr;
    return Status::OK();
}

class BCFFileSink {
    bool open_ = true;
    const string& filename_;
    bcf_hdr_t* header_;
    vcfFile *outfile_;;

    BCFFileSink(const std::string& filename, bcf_hdr_t* hdr, vcfFile* outfile)
        : filename_(filename), header_(hdr), outfile_(outfile)
        {}

public:
    static Status Open(const genotyper_config& cfg,
                       const string& filename,
                       bcf_hdr_t* hdr,
                       unique_ptr<BCFFileSink>& ans) {

        vcfFile* outfile;
        if (cfg.output_format == GLnexusOutputFormat::VCF) {
            // open as (uncompressed) vcf
            outfile = vcf_open(filename.c_str(), "w");
        } else if (cfg.output_format == GLnexusOutputFormat::BCF) {
            // open as bcf
            outfile = bcf_open(filename.c_str(), "wb");
        } else {
            return Status::Invalid("BCFFileSink::Open: Invalid output format");
        }
        if (!outfile) {
            return Status::IOError("failed to open BCF file for writing", filename);
        }
        if (bcf_hdr_write(outfile, hdr) != 0) {
            bcf_close(outfile);
            return Status::IOError("bcf_hdr_write", filename);
        }

        ans.reset(new BCFFileSink(filename, hdr, outfile));
        return Status::OK();
    }

    virtual ~BCFFileSink() {
        if (open_) {
            bcf_close(outfile_);
        }
    }

    virtual Status write(bcf1_t* record) {
        if (!open_) return Status::Invalid("BCFFilkSink::write() called on closed writer");
        return bcf_write(outfile_, header_, record) == 0
                ? Status::OK() : Status::IOError("bcf_write", filename_);

    }

    virtual Status close() {
        if (!open_) return Status::Invalid("BCFFileSink::close() called on closed writer");
        open_ = false;
        return bcf_close(outfile_) == 0
                ? Status::OK() : Status::IOError("bcf_close", filename_);
    }
};

// Figure out if there were any losses reported for this site
static bool any_losses(consolidated_loss &losses_for_site) {
    for (auto ls_pair : losses_for_site) {
        loss_stats &ls = ls_pair.second;
        if (ls.n_calls_lost > 0) return true;
    }
    return false;
}

Status Service::genotype_sites(const genotyper_config& cfg, const string& sampleset, const vector<unified_site>& sites, const string& filename, consolidated_loss& dlosses) {
    Status s;
    shared_ptr<const set<string>> samples;
    S(body_->metadata_->sampleset_samples(sampleset, samples));
    vector<string> sample_names(samples->begin(), samples->end());

    // create a BCF header for this sample set
    // TODO: make optional
    shared_ptr<bcf_hdr_t> hdr;
    S(prepare_bcf_header(body_->metadata_->contigs(), sample_names, hdr));

    // open output BCF file
    unique_ptr<BCFFileSink> bcf_out;
    S(BCFFileSink::Open(cfg, filename, hdr.get(), bcf_out));

    // set up the residuals file if requested
    unique_ptr<Residuals> residuals = nullptr;
    if (cfg.output_residuals &&
        !cfg.residuals_file.empty()) {
        s = Residuals::Open(cfg.residuals_file, *(body_->metadata_), body_->data_,
                            sampleset, sample_names, residuals);
        if (s.bad())
            return Status::Invalid("Problem opening the residuals file ", cfg.residuals_file);
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
            Status ls = genotype_site(cfg, *(body_->metadata_), body_->data_, sites[i],
                                      sampleset, sample_names, hdr.get(), bcf, losses_for_site);
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
        consolidated_loss losses_for_site = move(results[i].second);

        if (s.ok() && s_i.ok()) {
            // if everything's OK, proceed to write the record
            assert(bcf_i);
            merge_loss_stats(losses_for_site, dlosses);
            s = bcf_out->write(bcf_i.get());
            if (s.bad()) {
                abort = true;
            }

            if (residuals != nullptr &&
                any_losses(losses_for_site)) {
                // write out a residuals loss record
                s = residuals->write_record(sites[i], hdr.get(), bcf_i.get());
                if (s.bad()) {
                    abort = true;
                }
            }
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
    return bcf_out->close();
}

}

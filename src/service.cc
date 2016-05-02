#include "service.h"
#include "data.h"
#include "unifier.h"
#include "genotyper.h"
#include "residuals.h"
#include "diploid.h"
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iostream>
#include <map>
#include <assert.h>
#include <tuple>
#include "ctpl_stl.h"

using namespace std;

namespace GLnexus{

// pImpl idiom
struct Service::body {
    service_config cfg_;
    BCFData& data_;
    std::unique_ptr<MetadataCache> metadata_;

    // thread pool for executing discover_alleles and genotype_sites operations
    ctpl::thread_pool threadpool_;

    // "meta" thread pool for driving multiple concurrent discover_alleles and
    // "genotype_sites operations
    ctpl::thread_pool metapool_;

    body(BCFData& data) : data_(data) {}
};

Service::Service(const service_config& cfg, BCFData& data) {
    body_ = make_unique<Service::body>(data);
    body_->cfg_ = cfg;
    if (body_->cfg_.threads == 0) {
        body_->cfg_.threads = thread::hardware_concurrency();
    }
    body_->threadpool_.resize(body_->cfg_.threads);
    body_->metapool_.resize(body_->cfg_.threads);
}

Service::~Service() = default;

Status Service::Start(const service_config& cfg, Metadata& metadata, BCFData& data,
                      unique_ptr<Service>& svc) {
    svc.reset(new Service(cfg, data));
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
    vector<vector<float>> copy_number;
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
            assert(record->n_allele >= 3);
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

            // calculate estimated allele copy numbers for this record
            // TODO: ideally we'd compute them only for relevant samples
            S(diploid::estimate_allele_copy_number(dataset_header.get(), record.get(), copy_number));
            #define round4(x) (roundf(x*10000.0f)/10000.0f)

            // FIXME -- minor potential bug -- double-counting copy number of
            // alleles that span multiple discovery ranges

            // create a discovered_alleles entry for each alt allele matching [ACGT]+
            // In particular this excludes gVCF <NON_REF> symbolic alleles
            bool any_alt = false;
            for (int i = 1; i < record->n_allele; i++) {
                float copy_number_i = 0.0;
                for (unsigned sample : dataset_relevant_samples) {
                    copy_number_i += copy_number[sample][i];
                }
                if (copy_number_i >= 0.5) { // TODO: configurable threshold
                    string aldna(record->d.allele[i]);
                    transform(aldna.begin(), aldna.end(), aldna.begin(), ::toupper);
                    if (aldna.size() > 0 && regex_match(aldna, regex_dna)) {
                        discovered_allele_info ai = { false, round4(copy_number_i) };
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
                    float ref_copy_number = 0.0;
                    for (unsigned sample : dataset_relevant_samples) {
                        ref_copy_number += copy_number[sample][0];
                    }
                    discovered_allele_info ai = { true, round4(ref_copy_number) };
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

Status Service::discover_alleles(const string& sampleset, const range& pos,
                                 discovered_alleles& ans, atomic<bool>* ext_abort) {
    // Find the data sets containing the samples in the sample set.
    shared_ptr<const set<string>> samples, datasets;
    vector<unique_ptr<RangeBCFIterator>> iterators;
    Status s;

    // Query for (iterators to) records overlapping pos in all the data sets.
    // We query with min_alleles=3 to get variant records only (excluding
    // reference confidence records which have 2 alleles)
    S(body_->data_.sampleset_range(*(body_->metadata_), sampleset, pos, 3,
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
            if (abort || (ext_abort && *ext_abort)) {
                abort = true;
                return Status::Aborted();
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
                                 vector<discovered_alleles>& ans, atomic<bool>* ext_abort) {
    atomic<bool> abort(false);
    vector<future<Status>> statuses;
    vector<discovered_alleles> results(ranges.size());

    size_t i = 0;
    for (const auto& range : ranges) {
        auto fut = body_->metapool_.push([&, i, range](int tid){
            if (abort || (ext_abort && *ext_abort)) {
                abort = true;
                return Status::Aborted();
            }

            discovered_alleles dsals;
            Status ls = discover_alleles(sampleset, range, dsals, &abort);
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
                                 const vector<string>& samples, const vector<retained_format_field> format_fields,
                                 shared_ptr<bcf_hdr_t>& ans) {
    vector<string> hdr_lines;
    hdr_lines.push_back("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    hdr_lines.push_back("##FORMAT=<ID=RNC,Number=G,Type=Character,Description=\"Reason for No Call in GT: . = n/a, M = Missing data, P = Partial data, D = insufficient Depth of coverage, - = unrepresentable overlapping deletion, L = Lost/unrepresentable allele (other than deletion), U = multiple Unphased variants present, O = multiple Overlapping variants present\">");
    for (auto& format_field : format_fields) {
        hdr_lines.push_back(format_field.description);
    }
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
            outfile = bcf_open(filename.c_str(), "wb1");
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

Status Service::genotype_sites(const genotyper_config& cfg, const string& sampleset,
                               const vector<unified_site>& sites,
                               const string& filename, consolidated_loss& dlosses,
                               atomic<bool>* ext_abort) {
    Status s;
    shared_ptr<const set<string>> samples;
    S(body_->metadata_->sampleset_samples(sampleset, samples));
    vector<string> sample_names(samples->begin(), samples->end());

    // create a BCF header for this sample set
    // TODO: make optional
    shared_ptr<bcf_hdr_t> hdr;
    S(prepare_bcf_header(body_->metadata_->contigs(), sample_names, cfg.liftover_fields, hdr));

    // open output BCF file
    unique_ptr<BCFFileSink> bcf_out;
    S(BCFFileSink::Open(cfg, filename, hdr.get(), bcf_out));

    // set up the residuals file
    unique_ptr<ResidualsFile> residualsFile = nullptr;
    string res_filename;
    if (filename != "-" && filename.find(".") > 0) {
        int lastindex = filename.find_last_of(".");
        string rawname = filename.substr(0, lastindex);
        res_filename = rawname + ".residuals.yml";
    } else {
        res_filename = "/tmp/residuals.yml";
    }
    if (cfg.output_residuals) {
        S(ResidualsFile::Open(res_filename, residualsFile));
    }

    // Enqueue processing of each site as a task on the thread pool.
    vector<future<Status>> statuses;
    vector<tuple<shared_ptr<bcf1_t>,consolidated_loss,shared_ptr<string>>> results(sites.size());
    // ^^^ results to be filled by side-effect in the individual tasks below.
    // We assume that by virtue of preallocating, no mutex is necessary to
    // use it as follows because writes and reads of individual elements are
    // serialized by the futures.
    atomic<size_t> results_retrieved(0);
    atomic<bool> abort(false);
    for (size_t i = 0; i < sites.size(); i++) {
        auto fut = body_->threadpool_.push([&, i](int tid){
            while (i > results_retrieved+4*body_->cfg_.threads) {
                // throttle worker thread if the results retrieval, below, is falling
                // too far behind. Otherwise memory usage would be unbounded because
                // the results have to be retrieved and written out before they can
                // be deallocated.
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }
            if (abort || (ext_abort && *ext_abort)) {
                abort = true;
                return Status::Aborted();
            }
            shared_ptr<string> residual_rec = nullptr;
            shared_ptr<bcf1_t> bcf;
            consolidated_loss losses_for_site;
            Status ls = genotype_site(cfg, *(body_->metadata_), body_->data_, sites[i],
                                      sampleset, sample_names, hdr.get(), bcf, losses_for_site,
                                      residualsFile != nullptr, residual_rec,
                                      &abort);
            if (ls.bad()) {
                return ls;
            }

            results[i] = make_tuple(move(bcf), losses_for_site, residual_rec);
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
        shared_ptr<bcf1_t> bcf_i = move(std::get<0>(results[i]));
        assert(std::get<0>(results[i]) == nullptr);
        consolidated_loss losses_for_site = move(std::get<1>(results[i]));
        shared_ptr<string> residual_rec =  move(std::get<2>(results[i]));
        assert(std::get<2>(results[i]) == nullptr);

        if (s.ok() && s_i.ok()) {
            // if everything's OK, proceed to write the record
            assert(bcf_i);
            merge_loss_stats(losses_for_site, dlosses);
            s = bcf_out->write(bcf_i.get());
            if (s.bad()) {
                abort = true;
            }
            else if (residual_rec != nullptr) {
                // We have a residuals record, write it to disk.
                s = residualsFile->write_record(*residual_rec);
                if (s.bad()) {
                    abort = true;
                }
            }
        } else if (s.ok() && s_i.bad()) {
            // record the first error, and tell remaining tasks to abort
            s = move(s_i);
            abort = true;
        }
        results_retrieved++;
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

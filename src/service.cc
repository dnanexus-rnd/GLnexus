#include "service.h"
#include "data.h"
#include "discovery.h"
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

    atomic<uint64_t> threads_stalled_ms_;

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
    body_->threads_stalled_ms_ = 0;
}

Service::~Service() = default;

Status Service::Start(const service_config& cfg, Metadata& metadata, BCFData& data,
                      unique_ptr<Service>& svc) {
    svc.reset(new Service(cfg, data));
    return MetadataCache::Start(metadata, svc->body_->metadata_);
}

Status Service::discover_alleles(const string& sampleset, const range& pos,
                                 unsigned& N, discovered_alleles& ans,
                                 atomic<bool>* ext_abort) {
    // Find the data sets containing the samples in the sample set.
    shared_ptr<const set<string>> samples, datasets;
    vector<unique_ptr<RangeBCFIterator>> iterators;
    Status s;
    N = 0;

    // Query for (iterators to) records overlapping pos in all the data sets.
    // We query for variant records only (excluding reference confidence records
    // which have only a symbolic ALT allele)
    bcf_predicate predicate = [](const bcf_hdr_t* hdr, bcf1_t* bcf, bool &retval) {
        if (bcf_unpack(bcf, BCF_UN_STR)) {
            return Status::IOError("bcf_unpack");
        }
        retval = !is_gvcf_ref_record(bcf);
        return Status::OK();
    };
    S(body_->data_.sampleset_range(*(body_->metadata_), sampleset, pos, predicate,
                                   samples, datasets, iterators));
    N = samples->size();

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
            Status ls = discover_alleles_from_iterator(*samples, pos, *raw_iter, dsals);
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
                                 unsigned& N, vector<discovered_alleles>& ans, atomic<bool>* ext_abort) {
    atomic<bool> abort(false);
    vector<future<Status>> statuses;
    vector<discovered_alleles> results(ranges.size());
    N = 0;

    size_t i = 0;
    for (const auto& range : ranges) {
        auto fut = body_->metapool_.push([&, i, range](int tid){
            if (abort || (ext_abort && *ext_abort)) {
                abort = true;
                return Status::Aborted();
            }

            discovered_alleles dsals;
            unsigned tmpN;
            Status ls = discover_alleles(sampleset, range, tmpN, dsals, &abort);
            if (ls.ok()) {
                if (i == 0) {
                    // tmpN should be the same across all ranges
                    N = tmpN;
                }
                results[i] = move(dsals);
            }
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
                                 const vector<retained_format_field> format_fields,
                                 const vector<string>& extra_header_lines,
                                 shared_ptr<bcf_hdr_t>& ans) {
    vector<string> hdr_lines;
    hdr_lines.push_back("##GLnexusVersion=" + string(GIT_REVISION));
    for (const auto& line : extra_header_lines) {
        hdr_lines.push_back(line);
    }
    hdr_lines.push_back("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency estimate for each alternate allele\">");
    hdr_lines.push_back("##INFO=<ID=AQ,Number=A,Type=Integer,Description=\"Allele Quality score reflecting evidence for each alternate allele (Phred scale)\">");
    hdr_lines.push_back("##FILTER=<ID=MONOALLELIC,Description=\"Site represents one ALT allele in a region with multiple variants that could not be unified into non-overlapping multi-allelic sites\">");
    hdr_lines.push_back("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    hdr_lines.push_back("##FORMAT=<ID=RNC,Number=2,Type=Character,Description=\"Reason for No Call in GT: . = n/a, M = Missing data, P = Partial data, I = gVCF input site is non-called, D = insufficient Depth of coverage, - = unrepresentable overlapping deletion, L = Lost/unrepresentable allele (other than deletion), U = multiple Unphased variants present, O = multiple Overlapping variants present, 1 = site is Monoallelic, no assertion about presence of REF or ALT allele\">");
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
                               const string& filename,
                               atomic<bool>* ext_abort) {
    Status s;
    shared_ptr<const set<string>> samples;
    S(body_->metadata_->sampleset_samples(sampleset, samples));
    vector<string> sample_names(samples->begin(), samples->end());

    // create a BCF header for this sample set
    // TODO: make optional
    shared_ptr<bcf_hdr_t> hdr;
    S(prepare_bcf_header(body_->metadata_->contigs(), sample_names, cfg.liftover_fields,
                         body_->cfg_.extra_header_lines, hdr));

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
    vector<tuple<shared_ptr<bcf1_t>,shared_ptr<string>>> results(sites.size());
    // ^^^ results to be filled by side-effect in the individual tasks below.
    // We assume that by virtue of preallocating, no mutex is necessary to
    // use it as follows because writes and reads of individual elements are
    // serialized by the futures.
    atomic<size_t> results_retrieved(0);
    atomic<bool> abort(false);
    for (size_t i = 0; i < sites.size(); i++) {
        auto fut = body_->threadpool_.push([&, i](int tid){
            if (abort || (ext_abort && *ext_abort)) {
                abort = true;
                return Status::Aborted();
            }

            uint64_t stalled_ms = 0;
            while (i > results_retrieved+4*body_->cfg_.threads) {
                // throttle worker thread if the results retrieval, below, is falling
                // too far behind. Otherwise memory usage would be unbounded because
                // the results have to be retrieved and written out before they can
                // be deallocated.
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
                stalled_ms += 10;
            }
            if (i < body_->cfg_.threads && results_retrieved == 0) {
                // throttle startup so that database cache can burn in
                std::this_thread::sleep_for(std::chrono::milliseconds(i*10));
                stalled_ms += i*10;
            }
            if (stalled_ms) body_->threads_stalled_ms_ += stalled_ms;

            shared_ptr<string> residual_rec = nullptr;
            shared_ptr<bcf1_t> bcf;
            Status ls = genotype_site(cfg, *(body_->metadata_), body_->data_, sites[i],
                                      sampleset, sample_names, hdr.get(), bcf,
                                      residualsFile != nullptr, residual_rec,
                                      &abort);
            if (ls.bad()) {
                return ls;
            }

            results[i] = make_tuple(move(bcf), residual_rec);
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
        shared_ptr<string> residual_rec =  move(std::get<1>(results[i]));
        assert(std::get<1>(results[i]) == nullptr);

        if (s.ok() && s_i.ok()) {
            // if everything's OK, proceed to write the record
            if (bcf_i) {
                s = bcf_out->write(bcf_i.get());
            }
            if (s.bad()) {
                abort = true;
            } else if (residual_rec != nullptr) {
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

uint64_t Service::threads_stalled_ms() const { return body_->threads_stalled_ms_; }

}

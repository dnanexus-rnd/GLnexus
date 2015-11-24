#include <assert.h>
#include <algorithm>
#include "residuals.h"
#include <kstring.h>
#include <vcf.h>

using namespace std;

namespace GLnexus {

// convert a BCF record into a string
static shared_ptr<string> bcf1_to_string(const bcf_hdr_t *hdr, const bcf1_t *bcf) {
    kstring_t kstr;
    //memset((void)&kstr, 0, sizeof(kstring_t));
    kstr.l = 0; // I -think- this is a safe way to initialize a kstring_t
    kstr.m = 0;
    kstr.s = NULL;

    vcf_format(hdr, bcf, &kstr);
    auto retval = make_shared<string>(kstr.s, kstr.l);

    // cleanup
    if (kstr.s != NULL)
        free(kstr.s);

    return retval;
}

// destructor
Residuals::~Residuals() {
    // close the residuals file
    ofs_.close();
}


Status Residuals::Open(std::string filename,
                       const MetadataCache& cache, BCFData& data,
                       const std::string& sampleset, const std::vector<std::string>& samples,
                       unique_ptr<Residuals> &ans) {
    ans = make_unique<Residuals>(filename, cache, data, sampleset, samples);
    ans->ofs_.open(filename, std::ofstream::out | std::ofstream::app);
    return Status::OK();
}


Status Residuals::write_record_(const unified_site& site,
                                const bcf_hdr_t *gl_hdr,
                                bcf1_t *gl_call) {
    Status s;

    // write the original bcf records
    ofs_ << "- original BCF records" << endl;

    shared_ptr<const set<string>> samples2, datasets;
    vector<unique_ptr<RangeBCFIterator>> iterators;
    S(data_.sampleset_range(cache_, sampleset_, site.pos,
                           samples2, datasets, iterators));
    assert(samples_.size() == samples2->size());

    // for each pertinent dataset
    for (const auto& dataset : *datasets) {
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

        ofs_ << "   - dataset: " << dataset << endl;
        for (auto rec : records) {
            ofs_ << "     " << bcf1_to_string(dataset_header.get(), rec.get()) << endl;
        }
    }

    // write the unified site
    ofs_ << "- unified site" << endl;
    const auto& contigs = cache_.contigs();
    YAML::Emitter yaml_out;
    S(site.yaml(contigs, yaml_out));
    ofs_ << yaml_out.c_str();

    // write the output bcf, the calls that GLnexus made.
    ofs_ << "- GLnexus call(s)" << endl;
    ofs_ <<  bcf1_to_string(gl_hdr, gl_call) << endl;

    // separator between calls
    ofs_ << endl;
    return Status::OK();
}

Status Residuals::write_record(const unified_site& site,
                               const bcf_hdr_t *gl_hdr,
                               bcf1_t *gl_call) {
    std::lock_guard<std::mutex> lock(mutex_);
    return write_record_(site, gl_hdr, gl_call);
}

}

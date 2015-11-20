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
    kstr.l = 0; // I -think- this is a safe way to initialize a kstring_t

    vcf_format(hdr, bcf, &kstr);
    auto retval = make_shared<string>(kstr.s, kstr.l);

    // cleanup
    if (kstr.s != NULL)
        free(kstr.s);

    return retval;
}

Status write_residuals_record(ofstream &ofs,
                               const MetadataCache& cache, BCFData& data,
                               const unified_site& site,
                               const std::string& sampleset, const vector<string>& samples,
                               const bcf_hdr_t *gl_hdr, shared_ptr<bcf1_t>& gl_call) {
    Status s;

    // write the original bcf records
    ofs << "- original BCF records" << endl;

    shared_ptr<const set<string>> samples2, datasets;
    vector<unique_ptr<RangeBCFIterator>> iterators;
    S(data.sampleset_range(cache, sampleset, site.pos,
                           samples2, datasets, iterators));
    assert(samples.size() == samples2->size());

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

        ofs << "   - dataset: " << dataset << endl;
        for (auto rec : records) {
            ofs << "     " << bcf1_to_string(dataset_header.get(), rec.get()) << endl;
        }
    }

    // write the unified site
    ofs << "- unified site" << endl;
    const auto& contigs = cache.contigs();
    YAML::Emitter yaml_out;
    S(site.yaml(contigs, yaml_out));
    ofs << yaml_out.c_str();

    // write the output bcf, the calls that GLnexus made.
    ofs << "- GLnexus call(s)" << endl;
    ofs <<  bcf1_to_string(gl_hdr, gl_call.get()) << endl;

    // separator between calls
    ofs << endl;
    return Status::OK();
}

}

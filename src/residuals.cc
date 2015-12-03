#include <assert.h>
#include <algorithm>
#include <kstring.h>
#include <vcf.h>
#include "BCFSerialize.h"
#include "residuals.h"

using namespace std;

namespace GLnexus {

static string concat(const std::vector<std::string>& samples) {
    return std::accumulate(samples.begin(), samples.end(), std::string(),
                    [](const std::string& a, const std::string& b) -> std::string {
                        return a + (a.length() > 0 ? "\t" : "") + b;
                    } );
}


// return the list of samples in this dataset, in the format
// of one string, with a tab as a delimiter. This matches BCF headers.
static string samples_from_dataset(const bcf_hdr_t *hdr) {
    vector<string> samples;

    int nsamples = bcf_hdr_nsamples(hdr);
    for (int i=0; i < nsamples; i++) {
        auto sample_name = string(bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i));
        samples.push_back(sample_name);
    }

    return concat(samples);
}

// destructor
Residuals::~Residuals() {}

Status Residuals::Open(const MetadataCache& cache, BCFData& data,
                       const std::string& sampleset, const std::vector<std::string>& samples,
                       unique_ptr<Residuals> &ans) {
    ans = make_unique<Residuals>(cache, data, sampleset, samples);
    return Status::OK();
}


Status Residuals::gen_record(const unified_site& site,
                             const bcf_hdr_t *gl_hdr,
                             const bcf1_t *gl_call,
                             std::string &ynode) {
    Status s;
    YAML::Emitter out;

    out << YAML::BeginMap;

    // write the original records
    out << YAML::Key << "input";
    out << YAML::BeginMap;
    out << YAML::Key << "body";
    out << YAML::BeginSeq;
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
                return Status::Failure("residuals: iterator returned unexpected dataset",
                                       this_dataset + " instead of " + dataset);
            }
            records.insert(records.end(), these_records.begin(), these_records.end());
        }

        for (auto rec : records) {
            out << YAML::BeginMap;
            out << YAML::Key << dataset;
            out << YAML::Literal << samples_from_dataset(dataset_header.get()) + "\n" + *(bcf1_to_string(dataset_header.get(), rec.get()));
            out << YAML::EndMap;
        }
    }
    out << YAML::EndSeq;
    out << YAML::EndMap;

    // write the unified site
    out << YAML::Key << "unified_site";
    const auto& contigs = cache_.contigs();
    S(site.yaml(contigs, out));


    // write the output bcf, the calls that GLnexus made. Also, add the samples matching the calls.
    out << YAML::Key << "output_vcf";
    out << YAML::Literal << concat(samples_) + "\n" + *(bcf1_to_string(gl_hdr, gl_call));

    out << YAML::EndMap;

    // Emit YAML format
    ynode = string(out.c_str());
    return Status::OK();
}

// destructor
ResidualsFile::~ResidualsFile() {
    // write EOF marker
    ofs_ << "..." << endl;

    ofs_.close();
}


static bool is_file_exist(const string &fileName) {
    std::ifstream infile(fileName);
    return infile.good();
}

Status ResidualsFile::Open(std::string filename,
                                  std::unique_ptr<ResidualsFile> &ans) {
    ans = make_unique<ResidualsFile>(filename);

    if (is_file_exist(filename)) {
        // Erase the existing file, it is a result of a previous run.
        // This avoids confusing results of old runs, with the current run, in case of failure.
        int rc = remove(filename.c_str());
        if (rc != 0)
            return Status::Invalid("Residuals file cannot be removed ", filename);
    }
    ans->ofs_.open(filename, std::ofstream::out | std::ofstream::app);
    if (ans->ofs_.fail())
        return Status::Invalid("Error opening file for append", filename);
    return Status::OK();
}

Status ResidualsFile::write_record(std::string &rec) {
    if (!ofs_.good())
        return Status::Invalid("File cannot be written to ", filename_);

    ofs_ << "---" << endl;
    ofs_ << rec << endl;
    ofs_ << endl;
    return Status::OK();
}

}

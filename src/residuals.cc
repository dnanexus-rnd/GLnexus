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


Status residuals_gen_record(const unified_site& site,
                            const bcf_hdr_t *gl_hdr,
                            const bcf1_t *gl_call,
                            const std::vector<DatasetSiteInfo> &sites,
                            const MetadataCache& cache,
                            const std::vector<std::string>& samples,
                            std::string &ynode) {
    Status s;
    YAML::Emitter out;
    out << YAML::BeginMap;
    out << YAML::Key << "input";
    out << YAML::BeginMap;
    out << YAML::Key << "body";
    out << YAML::BeginSeq;

    // for each dataset
    for (const auto& ds_info : sites) {
        ostringstream records_text;
        // line with sample name(s)
        for (int i = 0; i < bcf_hdr_nsamples(ds_info.header); i++) {
            if (i > 0) {
                records_text << '\t';
            }
            records_text << bcf_hdr_int2id(ds_info.header, BCF_DT_SAMPLE, i);
        }
        // gVCF records, one per line
        for (const auto& rec : ds_info.records) {
            records_text << "\n" << *(bcf1_to_string(ds_info.header.get(), rec.get()));
        }
        out << YAML::BeginMap;
        out << YAML::Key << ds_info.name;
        out << YAML::Value << YAML::Literal << records_text.str();
        out << YAML::EndMap;
    }
    out << YAML::EndSeq;
    out << YAML::EndMap;

    // write the unified site
    out << YAML::Key << "unified_site";
    const auto& contigs = cache.contigs();
    S(site.yaml(contigs, out));


    // write the output bcf, the calls that GLnexus made. Also, add the samples matching the calls.
    out << YAML::Key << "output_vcf";
    out << YAML::Literal << concat(samples) + "\n" + *(bcf1_to_string(gl_hdr, gl_call));

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


Status ResidualsFile::Open(std::string filename,
                           std::unique_ptr<ResidualsFile> &ans) {
    ans = make_unique<ResidualsFile>(filename);

    // Note: we open in truncate mode, to erase the existing file, if any.
    // This avoids confusing results of old runs, with the current run, in case of failure.
    ans->ofs_.open(filename, std::ofstream::out | std::ofstream::trunc);
    if (ans->ofs_.fail())
        return Status::Invalid("Error opening file for truncate ", filename);
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

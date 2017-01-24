#pragma once

#include <string>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <assert.h>
#include <vcf.h>
#include <mutex>
#include <regex>
#include <atomic>
#include <algorithm>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
// suppressed warnings due to use of deprecated auto_ptr in yaml-cpp
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

#define UNPAIR(p,nm1,nm2) auto nm1 = (p).first; auto nm2 = (p).second;
template<typename T> inline void ignore_retval(T) {}

namespace GLnexus {

enum class StatusCode {
    OK,
    FAILURE,           // unspecified failure
    INVALID,           // invalid input/data
    NOT_FOUND,
    EXISTS,            // conflict with something that already exists
    IO_ERROR,
    NOT_IMPLEMENTED,
    ABORTED            // aborted per external request to do so
};

/// Function status (return) codes.

/// GLnexus functions generally return Status to indicate success or failure.
/// They throw exceptions only in truly "exceptional" circumstances (out of
/// memory or assertion violations)
class Status {
    StatusCode code_;
    const char *msg_;
    std::unique_ptr<std::string> noun_;

public:

    /// Default Status constructor.

    /// If provided, msg SHOULD be a string literal (not subject to moving or deallocation)
    Status(StatusCode code = StatusCode::OK, const char* msg = nullptr) noexcept
        : code_(code), msg_(msg) {}


    /// Extended Status constructor.

    /// If provided, msg SHOULD BE a string literal (not subject to moving or deallocation).
    /// The 'noun' will be copied and thus needs not be a literal.
    Status(StatusCode code, const char* msg, const std::string& noun)
        : code_(code), msg_(msg) {
        noun_ = std::make_unique<std::string>(noun);
    }

    // copy constructor
    Status(const Status &s) noexcept
        : code_(s.code_), msg_(s.msg_){
        if (s.noun_ != nullptr)
            noun_ = std::make_unique<std::string>(*s.noun_);
    }

    // assignment constructor
    Status& operator=(const Status &s) noexcept {
        // check for self-assignment
        if (&s == this)
            return *this;
        code_ = s.code_;
        msg_ = s.msg_;
        if (s.noun_ != nullptr)
            noun_ = std::make_unique<std::string>(*s.noun_);
        return *this;
    }

    bool ok() const noexcept { return code_ == StatusCode::OK; }
    bool bad() const noexcept { return !ok(); }
    operator StatusCode() const noexcept { return code_; }
    operator int() const noexcept { return static_cast<int>(code_); }

    static Status OK() noexcept { return Status(StatusCode::OK); }

    #define STATUS_SUGAR(nm,code) \
        static Status nm(const char* msg = nullptr) noexcept { return Status(code, msg); } \
        static Status nm(const char *msg, const std::string& noun) { return Status(code, msg, noun); }

    STATUS_SUGAR(Failure,StatusCode::FAILURE)
    STATUS_SUGAR(Invalid,StatusCode::INVALID)
    STATUS_SUGAR(NotFound,StatusCode::NOT_FOUND)
    STATUS_SUGAR(Exists,StatusCode::EXISTS)
    STATUS_SUGAR(IOError,StatusCode::IO_ERROR)
    STATUS_SUGAR(NotImplemented,StatusCode::NOT_IMPLEMENTED)
    STATUS_SUGAR(Aborted,StatusCode::ABORTED)

    std::string str() const {
        std::ostringstream ans;
        switch (code_) {
            case StatusCode::OK: ans << "OK"; break;
            case StatusCode::INVALID: ans << "Invalid"; break;
            case StatusCode::NOT_FOUND: ans << "NotFound"; break;
            case StatusCode::EXISTS: ans << "Exists"; break;
            case StatusCode::IO_ERROR: ans << "IOError"; break;
            case StatusCode::NOT_IMPLEMENTED: ans << "NotImplemented"; break;
            case StatusCode::ABORTED: ans << "Aborted"; break;
            default: ans << "Failure";
        }
        if (msg_) {
            ans << ": " << msg_;
            if (noun_) {
                ans << " (" << *noun_ << ")";
            }
        }
        return ans.str();
    }
};

// Convenience macro for re-raising a bad status when no recovery/cleanup is needed
#define S(st) s = st; if (s.bad()) return s;

/// common regular expressions
extern std::regex regex_dna, regex_id;

/// Genomic range (chromosome id, begin coordinate, end coordinate)
struct range {
    int rid=-1, beg=-1, end=-1;

    range(int rid_, int beg_, int end_) noexcept : rid(rid_), beg(beg_), end(end_) {
        if (beg_ > end_) {
            throw std::invalid_argument("invalid range (beginning > end)");
        }
    }

    /// Get the genomic range covered by a bcf1_t record. The end position is
    /// determined by the END INFO field if present (for structural variants and
    /// gVCF reference coverage records), or from the length of the reference
    /// allele otherwise.
    range(const bcf1_t* bcf) noexcept {
        assert(bcf != nullptr);
        rid = bcf->rid;
        beg = bcf->pos;
        end = bcf->pos+bcf->rlen;
    }

    range(const std::shared_ptr<const bcf1_t>& bcf) noexcept : range(bcf.get()) {}

    size_t size() const noexcept { return end-beg; }

    bool operator==(const range& r) const noexcept { return rid == r.rid && beg == r.beg && end == r.end; }
    bool operator!=(const range& r) const noexcept { return !(*this == r); }
    bool operator<(const range& r) const noexcept {
        return rid < r.rid || (rid == r.rid && (beg < r.beg || (beg == r.beg && end < r.end)));
    }
    bool operator<=(const range& r) const noexcept { return *this < r || *this == r; }

    bool overlaps(const range& r) const noexcept { return rid == r.rid && end > r.beg && beg < r.end; }
    bool within(const range& r) const noexcept { return rid == r.rid && beg >= r.beg && end <= r.end; }
    bool contains(const range& r) const noexcept { return r.within(*this); }

    bool spanned_by(const std::vector<range>& record_rngs) const noexcept {
        // Return false trivially when empty range given
        if (record_rngs.empty()) return false;

        std::vector<range> ranges(record_rngs);
        std::vector<range> merged_ranges;

        std::sort(ranges.begin(), ranges.end());
        range curr = ranges[0];
        for (auto& rng : ranges) {
            // Discontinuous region, start as a new range
            if(!curr.merge_contiguous(rng)) {
                merged_ranges.push_back(curr);
                curr = rng;
            }
        }
        merged_ranges.push_back(curr);

        for (auto& rng : merged_ranges) {
            if (rng.contains(*this)) return true;
        }

        return false;
    }

    bool contiguous(const range& r) const noexcept { return rid == r.rid && (r.beg == end || r.end == beg); }
    bool contigous_or_overlap(const range& r) const noexcept { return r.overlaps(*this) || r.contiguous(*this); }

    std::unique_ptr<range> intersect(const range& r) const {
        if (!overlaps(r)) return nullptr;
        return std::make_unique<range>(rid, std::max(beg,r.beg), std::min(end,r.end));
    }

    bool merge_contiguous(const range& r) {
        // Return false if the range to merge is not contiguous
        if (!contigous_or_overlap(r)) return false;

        // Return true with edited range
        beg = std::min(beg, r.beg);
        end = std::max(end, r.end);
        return true;
    }

    std::string str(const std::vector<std::pair<std::string,size_t> >& contigs) const {
        std::ostringstream os;
        if (rid >= 0 && rid < (int)contigs.size()) {
            os << std::get<0>(contigs[rid]);
        } else {
            os << '<' << rid << '>';
        }
        os << ':' << (beg+1) << '-' << end;
        return os.str();
    }
    std::string str() const {
        return str({});
    }
};

Status range_yaml(const std::vector<std::pair<std::string,size_t> >& contigs,
                  const range& r,
                  YAML::Emitter& yaml,
                  bool omit_ref = false);

Status range_of_yaml(const YAML::Node& yaml,
                     const std::vector<std::pair<std::string,size_t> >& contigs,
                     range& ans,
                     int default_rid = -1);

struct allele {
    range pos;
    std::string dna;

    allele(const range& pos_, const std::string& dna_) : pos(pos_), dna(dna_) {
        // Note; dna.size() may not equal pos.size(), for indel alleles
        if (!std::regex_match(dna, regex_dna)) throw std::invalid_argument("allele(): invalid DNA " + dna);
    }

    /// Equality is based on identity of position and allele
    bool operator==(const allele& rhs) const noexcept { return pos == rhs.pos && dna == rhs.dna; }
    bool operator!=(const allele& rhs) const noexcept { return !(*this == rhs); }

    /// Order is by pos and then allele.
    bool operator<(const allele& rhs) const noexcept { return pos < rhs.pos || (pos == rhs.pos && dna < rhs.dna); }
    bool operator<=(const allele& rhs) const noexcept { return *this < rhs || *this == rhs; }

    std::string str() const {
        std::ostringstream os;
        os << pos.str() << "(" << dna << ")";
        return os.str();
    }
};

// *Allele Quality (AQ)* of an allele in a VCF genotype call is defined in terms of
// the genotype likelihoods as follows: (the maximum likelihood of any genotype
// containing an allele / the maximum likelihood of any genotype not containing that
// allele), expressed on phred scale, truncated below zero.
//
// top_AQ is used to store the highest COUNT observations (descending order) of
// AQ for an allele across all genotype calls in the cohort.
const int MAX_AQ = 9999;
struct top_AQ {
    static const unsigned COUNT = 10;
    int V[COUNT] __attribute__ ((aligned));

    // This is a temporary buffer, used when adding observations. It
    // does not need to be serialized.
    std::vector<int> addbuf;

    top_AQ() {
        clear();
    }

    top_AQ(int AQ1) {
        clear();
        V[0] = AQ1;
    }

    void clear() {
        memset(&V, 0, sizeof(int)*COUNT);
        addbuf.clear();
    }

    void add(const int* rhs, const size_t rhs_count) {
        addbuf.resize(COUNT+rhs_count);
        memcpy(addbuf.data(), &V, COUNT*sizeof(int));
        memcpy(addbuf.data()+COUNT, rhs, rhs_count*sizeof(int));
        std::partial_sort(addbuf.begin(), addbuf.begin()+COUNT, addbuf.end(), std::greater<int>());
        memcpy(&V, addbuf.data(), COUNT*sizeof(int));
    }

    void operator+=(const top_AQ& rhs) {
        add(rhs.V, COUNT);
    }

    void operator+=(const std::vector<int>& rhs) {
        if (rhs.size()) {
            add(rhs.data(), rhs.size());
        }
    }

    bool operator==(const top_AQ& rhs) const {
        return memcmp(V, rhs.V, sizeof(int)*COUNT) == 0;
    }
};

// zygosity_by_GQ: holds information about how many times an allele is observed
// during the discovery process, stratified by genotype quality.
// A 10x2 matrix (M), the GQ_BANDS correspond to ten bands of phred-scaled GQ:
// 0 <= GQ < 10, 10 <= GQ < 20, ..., 90 <= GQ. The two columns correspond to
// allele zygosity (heterozygotes and homozygotes, or alelle copy number 1 & 2).
// The entries are how many genotype calls with the corresponding zygosity of
// the allele were observed in the cohort with GQ in the corresponding band.
//
// For example, z.M[5][1] is the number of homozygous calls with 50 <= GQ < 60.
struct zygosity_by_GQ {
    static const unsigned GQ_BANDS = 10;
    static const unsigned PLOIDY = 2;

    unsigned M[GQ_BANDS][PLOIDY] __attribute__ ((aligned));

    zygosity_by_GQ() {
        clear();
    }

    zygosity_by_GQ(unsigned zygosity, int GQ, unsigned count=1) {
        clear();
        add(zygosity, GQ, count);
    }

    void clear() {
        memset(&M, 0, sizeof(unsigned)*GQ_BANDS*PLOIDY);
    }

    void add(unsigned zygosity, int GQ, unsigned count=1) {
        assert(zygosity >= 1 && zygosity <= PLOIDY);
        unsigned i = std::min(unsigned(std::max(GQ, 0))/10U,GQ_BANDS-1U);
        M[i][zygosity-1] += count;
    }

    bool operator==(const zygosity_by_GQ& rhs) const {
        return memcmp(M, rhs.M, sizeof(int)*GQ_BANDS*PLOIDY) == 0;
    }

    void operator+=(const zygosity_by_GQ& rhs) {
        for (unsigned i = 0; i < GQ_BANDS; i++) {
            for (unsigned j = 0; j < PLOIDY; j++) {
                M[i][j] += rhs.M[i][j];
            }
        }
    }

    // estimate allele copy number in called genotypes with GQ >= minGQ
    unsigned copy_number(int minGQ = 0) const {
        unsigned ans = 0;
        unsigned i_lo = std::min(unsigned(std::max(minGQ, 0))/10U,GQ_BANDS-1U);

        for (unsigned i = i_lo; i < GQ_BANDS; i++) {
             for (unsigned j = 0; j < PLOIDY; j++) {
                ans += M[i][j]*(j+1);
            }
        }

        return ans;
    }
};

struct discovered_allele_info {
    bool is_ref = false;

    // top_AQ statistics are used to adjudicate allele existence
    top_AQ topAQ;

    // zygosity_by_GQ statsitics are used to estimate allele copy number
    zygosity_by_GQ zGQ;

    bool operator==(const discovered_allele_info& rhs) const noexcept {
        return is_ref == rhs.is_ref && topAQ == rhs.topAQ && zGQ == rhs.zGQ;
    }
    bool operator!=(const discovered_allele_info& rhs) const noexcept { return !(*this == rhs); }

    std::string str() const {
        std::ostringstream os;
        os << "[ is_ref: " << std::boolalpha << is_ref << " maxAQ: " << topAQ.V[0] << " copy number: " << zGQ.copy_number() << "]";
        return os.str();
    }
};
using discovered_alleles = std::map<allele,discovered_allele_info>;
Status merge_discovered_alleles(const discovered_alleles& src, discovered_alleles& dest);

Status yaml_of_one_discovered_allele(const allele& allele,
                                     const discovered_allele_info& ainfo,
                                     const std::vector<std::pair<std::string,size_t> >& contigs,
                                     YAML::Emitter& out);
Status one_discovered_allele_of_yaml(const YAML::Node&,
                                     const std::vector<std::pair<std::string,size_t> >& contigs,
                                     allele& allele,
                                     discovered_allele_info& ainfo);

Status yaml_of_discovered_alleles(const discovered_alleles&,
                                  const std::vector<std::pair<std::string,size_t> >& contigs,
                                  YAML::Emitter&);
Status discovered_alleles_of_yaml(const YAML::Node&,
                                  const std::vector<std::pair<std::string,size_t> >& contigs,
                                  discovered_alleles&);

 // Serialization of data structures with cap'n proto (https://capnproto.org/index.html)
//
// The issue this module tries to solve is that YAML serialization is slow for
// big data structures, and discovered-alleles tends to be large. Capnp serialization
// is very fast.
//
// write discovered_alleles structure to a file, with cap'n proto serialization
Status capnp_of_discovered_alleles(unsigned int sample_count,
                                   const std::vector<std::pair<std::string,size_t> >& contigs,
                                   const discovered_alleles &dsals,
                                   const std::string &filename);

// write discovered_alleles structure to a file descriptor.
//
// Note: this will not work with a C++ stream, only a low level file descriptor.
Status capnp_of_discovered_alleles_fd(unsigned int sample_count,
                                      const std::vector<std::pair<std::string,size_t> >& contigs,
                                      const discovered_alleles &dsals,
                                      int fd);

// read discovered_alleles structure from a file, as serialized by cap'n proto
Status discovered_alleles_of_capnp(const std::string &filename,
                                   unsigned int &sample_count,
                                   std::vector<std::pair<std::string,size_t> >& contigs,
                                   discovered_alleles &dsals);

// Verify that we can serialize and deserialize a discovered-alleles
// structure. Temporary results are written to [filename].
//
// Note: this is a debugging function
Status capnp_discover_alleles_verify(unsigned int sample_count,
                                     const std::vector<std::pair<std::string,size_t> >& contigs,
                                     const discovered_alleles &dsals,
                                     const std::string &filename);


struct unified_site {
    range pos;

    /// Optional: the sequencing target range (e.g. exon) containing this site
    range containing_target;


    /// Alleles at the position.

    /// Each allele is a string over [ACTGN]+. The first allele is the
    /// reference. The order of the remaining alleles is arbitrary.
    std::vector<std::string> alleles;


    /// Mapping of overlapping alleles (reference begin, reference end, allele
    /// DNA) onto the unified alleles (by index).
    std::map<allele,int> unification;

    std::vector<float> allele_frequencies;
    float lost_allele_frequency = 0.0f;

    // variant QUAL score (as in VCF)
    int qual = 0;

    bool operator==(const unified_site& rhs) const noexcept {
        if (!(pos == rhs.pos && alleles == rhs.alleles && unification == rhs.unification
              && allele_frequencies.size() == rhs.allele_frequencies.size()
              && lost_allele_frequency == rhs.lost_allele_frequency
              && qual == rhs.qual)) {
            return false;
        }
        // nan-tolerant comparison of allele_frequencies
        for (unsigned i = 0; i < allele_frequencies.size(); i++) {
            if (allele_frequencies[i] == allele_frequencies[i]) {
                if (allele_frequencies[i] != rhs.allele_frequencies[i]) {
                    return false;
                }
            } else {
                if (rhs.allele_frequencies[i] == rhs.allele_frequencies[i]) {
                    return false;
                }
            }
        }
        return true;
    }
    bool operator<(const unified_site& rhs) const noexcept{
        if (pos != rhs.pos) return pos < rhs.pos;
        if (alleles != rhs.alleles) return alleles < rhs.alleles;
        if (unification != rhs.unification) return unification < rhs.unification;
        return allele_frequencies < rhs.allele_frequencies;
    }

    unified_site(const range& pos_) noexcept : pos(pos_), containing_target(-1,-1,-1) {}

    Status yaml(const std::vector<std::pair<std::string,size_t> >& contigs,
                YAML::Emitter& out) const;
    static Status of_yaml(const YAML::Node& yaml,
                          const std::vector<std::pair<std::string,size_t> >& contigs,
                          unified_site& ans);
};

// Statistics collected during range queries
struct StatsRangeQuery {
    int64_t nBCFRecordsRead;    // how many BCF records were read from the DB
    int64_t nBCFRecordsInRange; // how many were in the requested range

    // constructor
    StatsRangeQuery() {
        nBCFRecordsRead = 0;
        nBCFRecordsInRange = 0;
    }

    // copy constructor
    StatsRangeQuery(const StatsRangeQuery &srq) {
        nBCFRecordsRead = srq.nBCFRecordsRead;
        nBCFRecordsInRange = srq.nBCFRecordsInRange;
    }

    // Addition
    StatsRangeQuery& operator+=(const StatsRangeQuery& srq) {
        nBCFRecordsRead += srq.nBCFRecordsRead;
        nBCFRecordsInRange += srq.nBCFRecordsInRange;
        return *this;
    }

    // return a human readable string
    std::string str() {
        std::ostringstream os;
        os << "Num BCF records read " << std::to_string(nBCFRecordsRead)
           << "  query hits " << std::to_string(nBCFRecordsInRange);
        return os.str();
    }
};

enum class UnifierPreference { Common, Small };

struct unifier_config {
    // AQ phred score thresholds: the unifier will include alleles having any
    // observation with AQ > min_AQ1, or having multiple observations with
    // AQ > min_AQ2 (min_AQ1 >= min_AQ2).
    //
    // All else equal, increasing min_AQ will increase specificity and reduce
    // sensitivity, and also speed up the genotyper (as fewer weak sites will
    // be considered)
    int min_AQ1 = 0, min_AQ2 = 0;

    // GQ phred score threshold for an input genotype call to "count" towards
    // copy number estimates for the constituent alleles.
    // Suggested value: = min_AQ2
    int min_GQ = 0;

    // Keep only alleles with at least this estimated copy number discovered
    // in the cohort.
    int min_allele_copy_number = 1;

    /// Maximum number of alleles per unified site; excess alleles will be
    /// pruned. If zero, then no specific limit is enforced.
    size_t max_alleles_per_site = 0;

    /// The unifier may need to prune alleles if they're too numerous and/or
    /// overlap in conflicting ways. The preference controls which alleles the
    /// unifier will try hardest to keep: common alleles (default), or alleles
    /// editing the smallest portion of the reference (least likely to
    /// conflict with other alleles).
    UnifierPreference preference = UnifierPreference::Common;

    bool operator==(const unifier_config& rhs) const noexcept {
        return min_allele_copy_number == rhs.min_allele_copy_number &&
            max_alleles_per_site == rhs.max_alleles_per_site &&
            preference == rhs.preference;
    }

    Status yaml(YAML::Emitter& out) const;
    static Status of_yaml(const YAML::Node& yaml, unifier_config& ans);
};

enum class GLnexusOutputFormat {
    /// Compressed bcf (default option)
    BCF,

    /// Uncompressed vcf (for ease of comparison in small cases)
    VCF,
};

enum class RetainedFieldType {
    INT,
    FLOAT,
};

enum class FieldCombinationMethod {
    MIN,
    MAX,
};

enum class RetainedFieldNumber {
    // In vcf specification, corresponds to Number=<numeric>
    BASIC,

    // In vcf specification, corresponds to Number=A
    ALT,

    // In vcf specification, corresponds to Number=G
    GENOTYPE,

    // Not a standard vcf specification, corresponds to 1 value per allele (REF + ALT)
    // Example, GATK's AD field (which is represented in vcf spec as Number=.)
    ALLELES,

};

enum class DefaultValueFiller {
    // Fill with "missing" values as default (applicable to most)
    MISSING = 0,

    // Fill with 0 values as default (applicable to AD)
    ZERO,

    // Do not fill with default values (FORMAT field will be dropped
    // when >=1 sample(s) have no values according to VCF spec)
    NONE
};

struct retained_format_field {
    // Original names for this lifted-over field, found in input gvcf
    // Acepts a vector as confidence records and vcf records may use
    // different field names for similar 'semantic' information
    std::vector<std::string> orig_names;

    // name of the field to be used in output
    std::string name;

    // description of format field to be inserted into header
    std::string description;

    // Value type of retained field (int, float)
    RetainedFieldType type;

    // NUMBER of INFO/FORMAT field,
    // translated from VCF spec
    RetainedFieldNumber number;

    // Applicable when number == BASIC
    int count;

    // Handling of "missing" values by using default
    // values, or leaving as empty, as instructed
    DefaultValueFiller default_type;

    // Handling of "combining" the same field
    // from multiple records
    FieldCombinationMethod combi_method;

    // Ignore non-variant records (as produced by genotyper::find_variant_records)
    bool ignore_non_variants;

    // Constructor
    retained_format_field(std::vector<std::string> orig_names_, std::string name_, RetainedFieldType type_,
        FieldCombinationMethod combi_method_, RetainedFieldNumber number_, int count_=0,
        DefaultValueFiller default_type_=DefaultValueFiller::MISSING, bool ignore_non_variants_=false)
        : orig_names(orig_names_), name(name_), type(type_), number(number_), count(count_), default_type(default_type_), combi_method(combi_method_), ignore_non_variants(ignore_non_variants_) {
        // Keep the names in sorted order, so that the comparison operator
        // will compare orig_name vectors element-wise.
        std::sort(orig_names.begin(), orig_names.end());
    };

    Status yaml(YAML::Emitter &out) const;
    static Status of_yaml(const YAML::Node& yaml, std::unique_ptr<retained_format_field>& ans);
};

struct genotyper_config {
    /// Use genotype likelihoods and unified allele frequencies to revise
    /// genotype calls
    bool revise_genotypes = false;

    /// Minimum assumed allele frequency to use in genotype revision; increasing
    /// this typically increases sensitivity to borderline ALT allele calls.
    /// Suggested value: 1/(2N) but not less than 0.0001
    float min_assumed_allele_frequency = 0.0001;

    /// Require any allele call to be supported by at least this depth
    size_t required_dp = 0;

    /// FORMAT field to consult for per-allele depth in VCF records
    std::string allele_dp_format = "AD";

    /// FORMAT field to consult for reference depth in gVCF reference records
    std::string ref_dp_format = "MIN_DP";

    // Should the genotyper write a record describing each call loss?
    // If true, the output is recorded in YAML format in the
    // a file named [BCF/VCF output file].residuals.yml
    bool output_residuals = false;

    /// Output format (default = bcf), choices = "BCF", "VCF"
    GLnexusOutputFormat output_format = GLnexusOutputFormat::BCF;

    // FORMAT fields from the original gvcfs to be lifted over to the output
    std::vector<retained_format_field> liftover_fields;

    genotyper_config() = default;

    genotyper_config(GLnexusOutputFormat _output_format) : output_format(_output_format) {}

    Status yaml(YAML::Emitter& out) const;
    static Status of_yaml(const YAML::Node& yaml, genotyper_config& ans);
};

// convenience wrapper for a self-freeing vector with an exposed 'capacity' --
// used with htslib functions that reuse/realloc the buffer
template<class T> struct htsvecbox {
    T *v = nullptr;
    int capacity = 0; // in number of elements, NOT bytes
    bool empty() const { return v == nullptr; }
    T& operator[](unsigned i) { return v[i]; }
    void clear() {
        free(v);
        v = nullptr;
        capacity = 0;
    }
    ~htsvecbox() {
        free(v);
    }
};

// test string against <.*>
bool is_symbolic_allele(const char*);

/// Determine whether the given record is a gVCF reference confidence record
/// (or else a "normal" record with at least one specific ALT allele)
bool is_gvcf_ref_record(const bcf1_t* record);

// Predicate function used for filtering BCF records, as they are read from the database.
// [retval] is set to true, for any record that passes the test.
//
// Note: the BCF record may be provided in packed form.  The function
// can unpack it, and return bad status in case of error (e.g., data
// corruption).
typedef Status (*bcf_predicate)(const bcf_hdr_t*, bcf1_t*, bool &retval);

} //namespace GLnexus

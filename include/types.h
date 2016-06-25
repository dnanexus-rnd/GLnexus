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

// zygosity_by_GQ: holds information about how many times an allele is observed
// during the discovery process, stratified by genotype quality.
// A 10x2 matrix (M), the rows correspond to ten bands of phred-scaled GQ:
// 0 <= GQ < 10, 10 <= GQ < 20, ..., 90 <= GQ. The two columns correspond to
// allele zygosity (heterozygotes and homozygotes, or alelle copy number 1 & 2).
// The entries are how many genotype calls with the corresponding zygosity of
// the allele were observed in the cohort with GQ in the corresponding band.
//
// For example, z.M[5][1] is the number of homozygous calls with 50 <= GQ < 60.
struct zygosity_by_GQ {
    static const unsigned ROWS = 10;
    static const unsigned COLS = 2;

    unsigned M[ROWS][COLS] __attribute__ ((aligned));

    zygosity_by_GQ() {
        memset(&M, 0, sizeof(int)*ROWS*COLS);
    }

    void add(unsigned zygosity, int GQ) {
        assert(zygosity >= 1 && zygosity <= COLS);
        unsigned i = std::min(unsigned(std::max(GQ, 0))/10U,ROWS-1U);
        M[i][zygosity-1]++;
    }

    bool operator==(const zygosity_by_GQ& rhs) const {
        return memcmp(M, rhs.M, sizeof(int)*ROWS*COLS) == 0;
    }

    void operator+=(const zygosity_by_GQ& rhs) {
        for (unsigned i = 0; i < ROWS; i++) {
            for (unsigned j = 0; j < COLS; j++) {
                M[i][j] += rhs.M[i][j];
            }
        }
    }

    // estimate allele copy number in called genotypes with GQ >= minGQ 
    unsigned copy_number(int minGQ = 0) const {
        unsigned ans = 0;
        unsigned i_lo = std::min(unsigned(std::max(minGQ, 0))/10U,ROWS-1U);

        for (unsigned i = i_lo; i < ROWS; i++) {
             for (unsigned j = 0; j < COLS; j++) {
                ans += M[i][j]*(j+1);
            }
        }

        return ans;
    }
};

struct discovered_allele_info {
    bool is_ref = false;

    // *Allele Quality (AQ)* of an allele in a VCF genotype call is defined in terms of
    // the genotype likelihoods as follows: (the maximum likelihood of any genotype
    // containing an allele / the maximum likelihood of any genotype not containing that
    // allele), expressed on phred scale, truncated below zero.
    //
    // maxAQ is the maximum AQ observed for this allele across all genotype calls in the
    // cohort.
    int maxAQ = 0;

    zygosity_by_GQ zGQ;

    bool operator==(const discovered_allele_info& rhs) const noexcept {
        return is_ref == rhs.is_ref && maxAQ == rhs.maxAQ && zGQ == rhs.zGQ;
    }
    bool operator!=(const discovered_allele_info& rhs) const noexcept { return !(*this == rhs); }
};
using discovered_alleles = std::map<allele,discovered_allele_info>;
Status merge_discovered_alleles(const discovered_alleles& src, discovered_alleles& dest);

Status yaml_of_discovered_alleles(const discovered_alleles&,
                                  const std::vector<std::pair<std::string,size_t> >&,
                                  YAML::Emitter&);
Status discovered_alleles_of_yaml(const YAML::Node&,
                                  const std::vector<std::pair<std::string,size_t> >&,
                                  discovered_alleles&);

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

    std::vector<float> copy_number;
    //std::vector<float> genotype_prior;

    bool operator==(const unified_site& rhs) const noexcept {
        return pos == rhs.pos && alleles == rhs.alleles && unification == rhs.unification
               && copy_number == rhs.copy_number;
    }
    bool operator<(const unified_site& rhs) const noexcept{
        if (pos != rhs.pos) return pos < rhs.pos;
        if (alleles != rhs.alleles) return alleles < rhs.alleles;
        if (unification != rhs.unification) return unification < rhs.unification;
        return copy_number < rhs.copy_number;
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
    // Allele Quality threshold for consideration of an allele (phred-scaled)
    int minAQ = 0;

    // Keep only alleles with at least this estimated copy number discovered
    // in the cohort. The estimated copy number is a soft estimate based on
    // the genotype likelihoods, so setting this somewhere between 0 and 1 can
    // filter out weak singleton observations.
    float min_allele_copy_number = 0.5;

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

struct retained_format_field {
    // Original names for this lifted-over field, found in input gvcf
    // Acepts a vector as confidence records and vcf records may use
    // different field names for similar 'semantic' information
    std::vector<std::string> orig_names;

    // name of the field to be used in output
    std::string name;

    // description of format field to be inserted into header
    std::string description;

    RetainedFieldType type;

    RetainedFieldNumber number;

    // We use a bool here to sidestep the trouble of setting
    // the appropriate type for float, int default values
    bool default_to_zero;

    // Applicable when number ==BASIC
    int count;

    FieldCombinationMethod combi_method;

    retained_format_field(std::vector<std::string> orig_names_, std::string name_, RetainedFieldType type_,
        FieldCombinationMethod combi_method_, RetainedFieldNumber number_, int count_=0, bool default_to_zero_=false)
        : orig_names(orig_names_), name(name_), type(type_), number(number_), default_to_zero(default_to_zero_), count(count_), combi_method(combi_method_) {
        // Keep the names in sorted order, so that the comparison operator
        // will compare orig_name vectors element-wise.
        std::sort(orig_names.begin(), orig_names.end());
    }

    Status yaml(YAML::Emitter &out) const;
    static Status of_yaml(const YAML::Node& yaml, std::unique_ptr<retained_format_field>& ans);
};

struct genotyper_config {
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

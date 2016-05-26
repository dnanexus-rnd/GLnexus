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

    /// Order is by pos and then allele.
    bool operator<(const allele& rhs) const noexcept { return pos < rhs.pos || (pos == rhs.pos && dna < rhs.dna); }
    bool operator<=(const allele& rhs) const noexcept { return *this < rhs || *this == rhs; }

    std::string str() const {
        std::ostringstream os;
        os << pos.str() << "(" << dna << ")";
        return os.str();
    }
};

struct discovered_allele_info {
    bool is_ref;
    float copy_number;

    bool operator==(const discovered_allele_info& rhs) const noexcept { return is_ref == rhs.is_ref && copy_number == rhs.copy_number; }

    std::string str() const {
        std::ostringstream os;
        os << "[ is_ref: " << std::boolalpha << is_ref << " copy number: " << copy_number << "]";
        return os.str();
    }
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


struct loss_stats {
    int n_calls_total=0, n_bp_total=0;
    int n_gvcf_calls_total=0, n_gvcf_bp_total=0;
    int n_no_calls_total = 0;

    // captures loss of both vcf entries and gvcf entries
    int n_calls_lost=0, n_bp_lost=0;

    // loss specific to gvcf entries
    int n_gvcf_calls_lost=0, n_gvcf_bp_lost=0;

    // Merges another loss_stats and increment the count
    // variables accordingly
    void operator+=(const loss_stats& loss) {
        n_calls_total += loss.n_calls_total;
        n_bp_total += loss.n_bp_total;

        n_gvcf_calls_total += loss.n_gvcf_calls_total;
        n_gvcf_bp_total += loss.n_gvcf_bp_total;

        n_calls_lost += loss.n_calls_lost;
        n_bp_lost += loss.n_bp_lost;

        n_gvcf_calls_lost += loss.n_gvcf_calls_lost;
        n_gvcf_bp_lost += loss.n_gvcf_bp_lost;

        n_no_calls_total += loss.n_no_calls_total;
    }

    std::string str() const noexcept {
        std::ostringstream ans;

        ans << n_no_calls_total << " no call(s).\n";

        // stop here if no no calls
        if (!n_no_calls_total)
            return ans.str();

        ans << "This is made up of a loss of " << n_calls_lost << " original call(s) which cover " << n_bp_lost << " bp.\n";
        ans << "The loss is " <<  std::setprecision(3) << prop_calls_lost() << "% of " << n_calls_total << " calls; or " << prop_bp_lost() << "% of " << n_bp_total << " bp processed from the original dataset(s).\n";

        ans << "Looking at gvcf entries, there is a loss of " << n_gvcf_calls_lost << " call(s) which cover " << n_gvcf_bp_lost << " bp.\n";
        ans << "The loss is " << prop_gvcf_calls_lost() << "% of " << n_gvcf_calls_total << " calls; or " << prop_gvcf_bp_lost() << "% of " << n_gvcf_bp_total << " bp of gvcf entries processed.\n";


        return ans.str();
    }

    // Returns proportion of calls lost as a percentage
    float prop_calls_lost() const noexcept {
        if (!n_calls_total) return 0;
        return n_calls_lost / (float) n_calls_total * 100;
    }

    // Returns proportion of bp coverage lost as a percentage
    float prop_bp_lost() const noexcept {
        if (!n_bp_total) return 0;
        return n_bp_lost / (float) n_bp_total * 100;
    }

    // Returns proportion of gvcf bp coverage lost as a percentage
    float prop_gvcf_bp_lost() const noexcept {
        if (!n_gvcf_bp_total) return 0;
        return n_gvcf_bp_lost / (float) n_gvcf_bp_total * 100;
    }

    // Returns proportion of gvcf calls lost as a percentage
    float prop_gvcf_calls_lost() const noexcept {
        if (!n_gvcf_calls_total) return 0;
        return n_gvcf_calls_lost / (float) n_gvcf_calls_total * 100;
    }
};

// per-sample loss_stats
using consolidated_loss = std::map<std::string, loss_stats>;
Status merge_loss_stats(const consolidated_loss& src, consolidated_loss& dest);

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
    // Keep only alleles with at least this estimated copy number discovered
    // in the cohort. The estimated copy number is a soft estimate based on
    // the genotype likelihoods, so setting this somewhere between 0 and 1 can
    // filter out weak singleton observations.
    float min_allele_copy_number = 0.0;

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
        : orig_names(orig_names_), name(name_), type(type_), number(number_), default_to_zero(default_to_zero_), count(count_), combi_method(combi_method_) {}

    bool operator==(const retained_format_field& rhs) const noexcept {
        // compare the liftover_fields vector
        std::vector<std::string>::const_iterator it, itrhs;
        for (it = orig_names.begin(), itrhs = rhs.orig_names.begin();
             it != orig_names.end() && itrhs != rhs.orig_names.end();
             ++it, ++itrhs) {
            if ((*it) != (*itrhs)) return false;
        }

        return name == rhs.name &&
            description == rhs.description &&
            type == rhs.type &&
            number == rhs.number &&
            default_to_zero == rhs.default_to_zero &&
            count == rhs.count &&
            combi_method == rhs.combi_method;
    }

    Status yaml(YAML::Emitter &out) const;
    static Status of_yaml(const YAML::Node& yaml, std::unique_ptr<retained_format_field>& ans);
};

struct genotyper_config {
    /// Require any allele call to be supported by at least this depth
    size_t required_dp = 0;

    /// FORMAT field to consult for per-allele depth in VCF records
    std::string allele_dp_format = "AD";

    /// The symbolic allele used in gVCF reference confidence models
    std::string ref_symbolic_allele = "<NON_REF>";

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

    bool operator==(const genotyper_config& rhs) const noexcept {
        // compare the liftover_fields vector
        std::vector<retained_format_field>::const_iterator it, itrhs;
        for (it = liftover_fields.begin(), itrhs = rhs.liftover_fields.begin();
             it != liftover_fields.end() && itrhs != rhs.liftover_fields.end();
             ++it, ++itrhs) {
            if (!(*it == *itrhs)) return false;
        }

        return (required_dp == rhs.required_dp &&
                allele_dp_format == rhs.allele_dp_format &&
                ref_symbolic_allele == rhs.ref_symbolic_allele &&
                ref_dp_format == rhs.ref_dp_format &&
                output_residuals == rhs.output_residuals &&
                output_format == rhs.output_format);
    }

    Status yaml(YAML::Emitter& out) const;
    static Status of_yaml(const YAML::Node& yaml, genotyper_config& ans);
};

// convenience wrapper for a self-freeing vector with an exposed 'capacity' --
// used with htslib functions that reuse/realloc the buffer
template<class T> struct htsvecbox {
    T *v = nullptr;
    int capacity = 0;
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

// Predicate function used for filtering BCF records, as they are read from the database.
// [retval] is set to true, for any record that passes the test.
//
// Note: the BCF record may be provided in packed form.  The function
// can unpack it, and return bad status in case of error (e.g., data
// corruption).
typedef Status (*bcf_predicate)(const bcf_hdr_t*, bcf1_t*, bool &retval);

} //namespace GLnexus

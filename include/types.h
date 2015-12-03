#ifndef GLNEXUS_TYPES_H
#define GLNEXUS_TYPES_H

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
#include "yaml-cpp/yaml.h"

#define UNPAIR(p,nm1,nm2) auto nm1 = (p).first; auto nm2 = (p).second;
template<typename T> inline void ignore_retval(T) {}

namespace GLnexus {

enum class StatusCode { OK, FAILURE, INVALID, NOT_FOUND, EXISTS, IO_ERROR, NOT_IMPLEMENTED };

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

    std::string str() const {
        std::ostringstream ans;
        switch (code_) {
            case StatusCode::OK: ans << "OK"; break;
            case StatusCode::INVALID: ans << "Invalid"; break;
            case StatusCode::NOT_FOUND: ans << "NotFound"; break;
            case StatusCode::EXISTS: ans << "Exists"; break;
            case StatusCode::IO_ERROR: ans << "IOError"; break;
            case StatusCode::NOT_IMPLEMENTED: ans << "NotImplemented"; break;
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

/// [AGCTN]+
bool is_dna(const std::string&);

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

    std::unique_ptr<range> intersect(const range& r) const {
        if (!overlaps(r)) return nullptr;
        return std::make_unique<range>(rid, std::max(beg,r.beg), std::min(end,r.end));
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
        if (!is_dna(dna)) throw std::invalid_argument("allele(): invalid DNA " + dna);
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
    float observation_count;

    bool operator==(const discovered_allele_info& rhs) const noexcept { return is_ref == rhs.is_ref && observation_count == rhs.observation_count; }

    std::string str() const {
        std::ostringstream os;
        os << "[ is_ref: " << std::boolalpha << is_ref << " observation count: " << std::setprecision(1) << observation_count << "]";
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


    /// Alleles at the position.

    /// Each allele is a string over [ACTGN]+. The first allele is the
    /// reference. The order of the remaining alleles is arbitrary.
    std::vector<std::string> alleles;


    /// Mapping of overlapping alleles (reference begin, reference end, allele
    /// DNA) onto the unified alleles (by index).
    std::map<allele,int> unification;

    std::vector<float> observation_count;
    //std::vector<float> genotype_prior;

    bool operator==(const unified_site& rhs) const noexcept {
        return pos == rhs.pos && alleles == rhs.alleles && unification == rhs.unification
               && observation_count == rhs.observation_count;
    }
    bool operator<(const unified_site& rhs) const noexcept{
        if (pos != rhs.pos) return pos < rhs.pos;
        if (alleles != rhs.alleles) return alleles < rhs.alleles;
        if (unification != rhs.unification) return unification < rhs.unification;
        return observation_count < rhs.observation_count;
    }

    unified_site(const range& pos_) noexcept : pos(pos_) {}

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

} //namespace GLnexus

#endif

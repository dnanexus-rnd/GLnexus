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
        // TODO: validate allele matches [ACTGN]+
    }

    /// Equality is based on identity of position and allele
    bool operator==(const allele& rhs) const noexcept { return pos == rhs.pos && dna == rhs.dna; }

    /// Order is by pos and then allele.
    bool operator<(const allele& rhs) const noexcept { return pos < rhs.pos || (pos == rhs.pos && dna < rhs.dna); }
    bool operator<=(const allele& rhs) const noexcept { return *this < rhs || *this == rhs; }
};

struct discovered_allele_info {
    bool is_ref;
    float observation_count;
};
using discovered_alleles = std::map<allele,discovered_allele_info>;
Status merge_discovered_alleles(const discovered_alleles& src, discovered_alleles& dest);

struct unified_site {
    range pos;


    /// Alleles at the position.

    /// Each allele is a string over [ACTGN]+. The first allele is the
    /// reference. The order of the remaining alleles is arbitrary.
    std::vector<std::string> alleles;


    /// Mapping of overlapping alleles (begin coordinate, DNA) onto the unified alleles (by index).
    std::map<std::pair<int,std::string>,int> unification;

    std::vector<float> observation_count;
    //std::vector<float> genotype_prior;

    bool operator==(const unified_site& rhs) const noexcept{ return pos == rhs.pos && alleles == alleles && observation_count == observation_count; }
    bool operator<(const unified_site& rhs) const noexcept{ return pos < rhs.pos; }
    bool operator<=(const unified_site& rhs) const noexcept{ return pos <= rhs.pos; }

    unified_site(const range& pos_) noexcept : pos(pos_) {}
};

// Simple wrapper struct to store information about an original call
// for a unified site, used by loss_stats
struct orig_call {
    range pos;
    bool is_gvcf;

    bool operator==(const orig_call& rhs) const noexcept { return pos == rhs.pos && is_gvcf == rhs.is_gvcf; }

    bool operator<(const orig_call& rhs) const noexcept { return pos < rhs.pos; }
    bool operator<=(const orig_call& rhs) const noexcept { return pos <= rhs.pos; }
};

class loss_stats {

public:

    // Constructor takes a sample name; loss_stats should be
    // uniquely associated with each sample
    loss_stats(const unified_site site_) noexcept : site(site_) {orig_calls_for_site = std::map<range, int>(); }


    // Getter operations
    int get_n_calls_total() const noexcept { return n_calls_total; }
    int get_n_bp_total() const noexcept { return n_bp_total; }
    int get_n_calls_lost() const noexcept { return n_calls_lost; }
    int get_n_bp_lost() const noexcept {return n_bp_lost; }
    int get_n_no_calls_total() const noexcept {return n_no_calls_total; }

    // Called for each bcf record that is associated with the unified_site
    // examined to update "denominator" of total calls/bp covered in orig
    // dataset
    Status add_call_for_site(const range call, int n_calls) noexcept {
        if (is_finalized)
            return Status::Invalid("calling add_call_for_site for a finalized loss_stats");

        auto call_within_site_p = call.intersect(site.pos);

        if (call_within_site_p) {
            range call_within_site = *call_within_site_p;
            orig_calls_for_site[call_within_site] += n_calls;
            n_calls_total += n_calls;
            n_bp_total += n_calls * (call_within_site.size());
        }
        return Status::OK();
    }

    // Called after joint genotyping for a unified_site to update the loss
    // associated with that site.
    Status finalize_loss_for_site(int n_no_calls) noexcept {
        if (is_finalized)
            return Status::Invalid("calling finalize_loss_for_site when loss_stats is already finalized.");

        if (n_no_calls) {
            n_no_calls_total += n_no_calls;

            for (auto& kv : orig_calls_for_site) {
                // calls_lost will be 0 if n_orig_calls and n_no_calls are both 1 (which is interpreted as no loss of calls)
                int n_calls_lost_for_site = (kv.second * n_no_calls) / 2;
                n_calls_lost += n_calls_lost_for_site;

                // assumes that range in calls_for_site has already been
                // restricted to intersection with site
                n_bp_lost += kv.first.size() * n_calls_lost_for_site;
            }
        }

        // Clear map
        orig_calls_for_site.clear();
        is_finalized = true;

        return Status::OK();
    }

    // Merges another loss_stats and increment the count
    // variables accordingly
    Status merge_loss(const loss_stats& loss) {
        if (!is_finalized || !loss.is_finalized) {
            return Status::Failure("loss_summary: trying to add a loss stats which has not been finalized.");
        }

        // Update summary loss stats
        n_calls_total += loss.get_n_calls_total();
        n_bp_total += loss.get_n_bp_total();

        n_calls_lost += loss.get_n_calls_lost();
        n_bp_lost += loss.get_n_bp_lost();

        n_no_calls_total += loss.get_n_no_calls_total();

        return Status::OK();
    }

    std::string str() const noexcept {
        std::ostringstream ans;

        ans << n_no_calls_total << " no call(s).\n";

        // stop here if no no calls
        if (!n_no_calls_total)
            return ans.str();

        ans << "This is made up of a loss of " << n_calls_lost << " original call(s) which cover " << n_bp_lost << " bp.\n";
        ans << "The loss is " <<  std::setprecision(3) << prop_calls_lost() << "% of " << n_calls_total << " calls; or " << prop_bp_lost() << "% of " << n_bp_total << " bp processed from the original dataset(s).\n";

        return ans.str();
    }

    bool is_finalized = false;

private:

    unified_site site;

    // Original calls (identified by effective range within site) and
    // count of calls. Calls which are different in the bcf record but
    // share the same effective range within the site will be collapsed
    // into the same key.
    std::map<range, int> orig_calls_for_site;

    int n_calls_total=0, n_bp_total=0;
    int n_calls_lost=0, n_bp_lost=0;
    int n_no_calls_total = 0;

    // Returns proportion of calls lost as a percentage
    float prop_calls_lost() const noexcept {
        return n_calls_lost / (float) n_calls_total * 100;
    }

    // Returns proportion of bp coverage lost as a percentage
    float prop_bp_lost() const noexcept {
        return n_bp_lost / (float) n_bp_total * 100;
    }
};

using consolidated_loss = std::map<std::string, loss_stats>;

Status merge_loss_stats(const consolidated_loss& src, consolidated_loss& dest);

// class loss_summary {

// public:
//     loss_summary(const std::string sample_name_) noexcept : sample_name(sample_name_) {}

//     std::string sample_name;

//     int get_n_calls_total() const noexcept { return n_calls_total; }
//     int get_n_bp_total() const noexcept { return n_bp_total; }
//     int get_n_calls_lost() const noexcept { return n_calls_lost; }
//     int get_n_bp_lost() const noexcept {return n_bp_lost; }
//     int get_n_no_calls_total() const noexcept {return n_no_calls_total; }

//     // Add a loss stats object to the total loss summary statistics
//     Status add_loss_stats(const loss_stats& loss) noexcept {
//         // Ensure that sample_name match
//         if (loss.sample_name != sample_name) {
//             return Status::Failure("loss_summary: trying to add a loss stats with sample_name differing from loss_summary.");
//         }
//         if (!loss.is_finalized) {
//             return Status::Failure("loss_summary: trying to add a loss stats which has not been finalized.");
//         }

//         // Update summary loss stats
//         n_calls_total += loss.n_calls_total;
//         n_bp_total += loss.n_bp_total;

//         n_calls_lost += loss.n_calls_lost;
//         n_bp_lost += loss.n_bp_lost;

//         n_no_calls_total += loss.n_no_calls_total;


//         return Status::OK();
//     }

//     // String representation of the loss observed for the sample
//     std::string str() const noexcept {
//         std::ostringstream ans;

//         ans << sample_name << ": ";
//         ans << n_no_calls_total << " no call(s).\n";

//         // stop here if no no calls
//         if (!n_no_calls_total)
//             return ans.str();

//         ans << "This is made up of a loss of " << n_calls_lost << " original call(s) which cover " << n_bp_lost << " bp.\n";
//         ans << "The loss is " <<  std::setprecision(3) << prop_calls_lost() << "% of " << n_calls_total << " calls; or " << prop_bp_lost() << "% of " << n_bp_total << " bp processed from the original dataset(s).\n";

//         return ans.str();
//     }


// private:
//     int n_calls_total=0, n_bp_total=0;
//     int n_calls_lost=0, n_bp_lost=0;
//     int n_no_calls_total = 0;

//     // Returns proportion of calls lost as a percentage
//     float prop_calls_lost() const noexcept {
//         return n_calls_lost / (float) n_calls_total * 100;
//     }

//     // Returns proportion of bp coverage lost as a percentage
//     float prop_bp_lost() const noexcept {
//         return n_bp_lost / (float) n_bp_total * 100;
//     }
// };

} //namespace GLnexus

#endif

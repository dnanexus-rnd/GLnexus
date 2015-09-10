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
#include <assert.h>
#include <vcf.h>

#define UNPAIR(p,nm1,nm2) auto nm1 = (p).first; auto nm2 = (p).second;

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

class loss_stats {

public:

    // Constructor takes a sample name; loss_stats should be
    // uniquely associated with each sample
    loss_stats(const std::string sample_name_) noexcept : sample_name(sample_name_) {}

    // Prepares the map with the given unified_site, requires that this method
    // is not called again when an entry for the site exists
    Status initializeSite(const unified_site site) noexcept {
        if (orig_calls_for_sites.find(site) != orig_calls_for_sites.end())
            return Status::Failure("loss_stats: initializeSite called for existing site");

        orig_calls_for_sites.insert(std::make_pair(site, std::map<range, int>() ));
        return Status::OK();
    }

    // Called for each bcf record that is associated with the unified_site
    // examined to update "denominator" of total calls/bp covered in orig
    // dataset
    Status addCallForSite(const unified_site site, const range call, int n_calls) noexcept {
        if (orig_calls_for_sites.find(site) == orig_calls_for_sites.end())
            return Status::Failure("loss_stats: addCallForSite called before initializeSite or after finalizeLossForSite.");

        auto& orig_calls_for_site = orig_calls_for_sites[site];

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
    Status finalizeLossForSite(const unified_site site, int n_no_calls) noexcept {
        if (orig_calls_for_sites.find(site) == orig_calls_for_sites.end())
            return Status::Failure("loss_stats: finalizeLossForSite called before initializeSite or after finalizeLossForSite.");

        if (n_no_calls) {
            n_no_calls_total += n_no_calls;

            auto& calls_for_site = orig_calls_for_sites[site];
            for (auto& kv : calls_for_site) {
                // calls_lost will be 0 if n_orig_calls and n_no_calls are both 1 (which is interpreted as no loss of calls)
                int n_calls_lost_for_site = (kv.second * n_no_calls) / 2;
                n_calls_lost += n_calls_lost_for_site;

                // assumes that range in calls_for_site has already been
                // restricted to intersection with site
                n_bp_lost += kv.first.size();
            }
        }
        orig_calls_for_sites.erase(site);

        return Status::OK();
    }

    // String representation of the loss observed for the sample
    std::string str() const noexcept {
        std::ostringstream ans;

        ans << sample_name << ": ";
        ans << n_no_calls_total << " no call(s).\n";

        // stop here if no no calls
        if (!n_no_calls_total)
            return ans.str();

        ans << "This is made up of a loss of " << n_calls_lost << " original call(s) which cover " << n_bp_lost << " bp.\n";
        ans << "The loss is " <<  std::setprecision(3) << propCallsLost() << "% of " << n_calls_total << " calls; or " << propBpLost() << "% of " << n_bp_total << " bp processed from the original dataset(s).\n";

        return ans.str();
    }

private:
    std::string sample_name;

    // Map of unified sites to a vector of original calls that is associated
    // with that unified site. The calls associated with a site will exist
    // for the duration between initializeSite and and finalizeLossForSite
    std::map<unified_site, std::map<range, int>> orig_calls_for_sites;

    int n_calls_total=0, n_bp_total=0;
    int n_calls_lost=0, n_bp_lost=0;
    int n_no_calls_total = 0;

    // Returns proportion of calls lost as a percentage
    float propCallsLost() const noexcept {
        return n_calls_lost / (float) n_calls_total * 100;
    }

    // Returns proportion of bp coverage lost as a percentage
    float propBpLost() const noexcept {
        return n_bp_lost / (float) n_bp_total * 100;
    }

};

} //namespace GLnexus

#endif

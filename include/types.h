#ifndef GLNEXUS_TYPES_H
#define GLNEXUS_TYPES_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <stdexcept>
#include <sstream>
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

    // TODO: if the bcf1_t INFO has an END field, that should be used instead of the reference
    // length to determine the end coordinate.
    range(const bcf1_t* bcf) : range(bcf->rid, bcf->pos, bcf->pos+bcf->rlen) {}
    range(const bcf1_t& bcf) : range(&bcf) {}
    range(const std::shared_ptr<const bcf1_t>& bcf) : range(bcf.get()) {}

    size_t size() const noexcept { return end-beg; }

    bool operator==(const range& r) const noexcept { return rid == r.rid && beg == r.beg && end == r.end; }
    bool operator<(const range& r) const noexcept { 
        return rid < r.rid || (rid == r.rid && (beg < r.beg || (beg == r.beg && end < r.end)));
    }
    bool operator<=(const range& r) const noexcept { return *this < r || *this == r; }

    bool overlaps(const range& r) const noexcept { return rid == r.rid && end > r.beg && beg < r.end; }
    bool within(const range& r) const noexcept { return rid == r.rid && beg >= r.beg && end <= r.end; }

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
        os << ':' << (beg+1) << '-' << (end+1);
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

/// table of allele and (is_ref,observation_count)
using discovered_alleles = std::map<allele,std::pair<bool,float> >;
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

    unified_site(const range& pos_) noexcept : pos(pos_) {}
};

} //namespace GLnexus

#endif

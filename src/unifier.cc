#include <algorithm>
#include <assert.h>
#include <math.h>
#include "unifier.h"
#include <iostream>

using namespace std;

namespace GLnexus {

using discovered_allele = pair<allele,discovered_allele_info>;

// ALT alleles discovered across numerous samples may represent the same edit
// to the reference genome in different ways, specifically if they have
// different amounts of reference padding on either end. (Individual gVCF ALT
// alleles may be reference-padded to unify their representation with other
// ALT alleles observed in the same sample.) So we introduce a notion of
// "minimized" ALT allele, which strips any reference padding so that
// equivalent ALT alleles can be combined, while still remembering the
// original representations.
struct minimized_allele_info {
    set<allele> originals;
    bool all_filtered = false;
    top_AQ topAQ;
    unsigned copy_number = 0;
    range in_target = range(-1,-1,-1);

    string str() const {
        ostringstream os;
        os << "Minimized from originals: " << endl;
        for (auto& al : originals) {
            os << "  " << al.str() << endl;
        }
        os << "Max AQ: " << topAQ.V[0] << endl;
        os << "Copy number: " << copy_number << endl;
        return os.str();
    }

    void operator+=(const minimized_allele_info& rhs) {
        originals.insert(rhs.originals.begin(), rhs.originals.end());
        all_filtered = all_filtered && rhs.all_filtered;
        topAQ += rhs.topAQ;
        copy_number += rhs.copy_number;
        if (!in_target.overlaps(rhs.in_target) ||
            in_target.size() < rhs.in_target.size()) {
            in_target = rhs.in_target;
        }
    }
};
using minimized_alleles = map<allele,minimized_allele_info>;
using minimized_allele = pair<allele,minimized_allele_info>;

// orders on minimized alleles corresponding to the UnifierPreferences

bool minimized_allele_delta_lt(const minimized_allele& p1, const minimized_allele& p2) {
    auto d1 = std::max(p1.first.dna.size(),p1.first.pos.size()) - std::min(p1.first.dna.size(),p1.first.pos.size());
    auto d2 = std::max(p2.first.dna.size(),p2.first.pos.size()) - std::min(p2.first.dna.size(),p2.first.pos.size());
    return d1 < d2;
}

// UnifierPreference::Small
bool minimized_allele_small_lt(const minimized_allele& p1, const minimized_allele& p2) {
    // smallest stretch of the reference
    if (p1.first.pos.size() != p2.first.pos.size()) {
        return p1.first.pos.size() < p2.first.pos.size();
    }
    if (p2.second.copy_number != p1.second.copy_number) {
        return p2.second.copy_number < p1.second.copy_number;
    }
    return p2.second.topAQ.V[0] < p1.second.topAQ.V[0];
}

// UnifierPreference::Common
bool minimized_allele_common_lt(const minimized_allele& p1, const minimized_allele& p2) {
    if (p2.second.copy_number != p1.second.copy_number) {
        return p2.second.copy_number < p1.second.copy_number;
    }
    if (p2.second.topAQ.V[0] != p1.second.topAQ.V[0]) {
        return p2.second.topAQ.V[0] < p1.second.topAQ.V[0];
    }
    return minimized_allele_delta_lt(p1, p2);
}

std::function<bool(const minimized_allele&, const minimized_allele&)> minimized_allele_lt(UnifierPreference pref) {
    if (pref == UnifierPreference::Small) {
        return minimized_allele_small_lt;
    }
    return minimized_allele_common_lt;
}

/// Minimize an ALT allele by trimming bases equal to the reference from
/// either end. This assumes indel alleles are already left-aligned.
Status minimize_allele(const allele& ref, allele& alt) {
    if (ref.pos != alt.pos) return Status::Invalid("minimize_allele: wrong ref allele provided");

    size_t l,r;

    // trim from the right first, to maintain indel left-alignment
    for (r=0; r<alt.dna.size()-1 && r<ref.dna.size()-1
              && alt.dna[alt.dna.size()-r-1] == ref.dna[ref.dna.size()-r-1]; r++);

    // then trim from the left if possible
    for (l=0; l<alt.dna.size()-r-1 && l<ref.dna.size()-r-1 && alt.dna[l] == ref.dna[l]; l++);
    assert(l+r < alt.dna.size());
    assert(l+r < ref.dna.size());

    alt.pos.beg += l;
    alt.pos.end -= r;
    assert(alt.pos.beg <= alt.pos.end);
    alt.dna = alt.dna.substr(l, alt.dna.size()-l-r);
    assert(alt.dna.size());

    return Status::OK();
}

// Separate discovered alleles into the REF alleles and minimized ALT alleles
Status minimize_alleles(const unifier_config& cfg, const discovered_alleles& src,
                        map<range,discovered_allele>& refs, minimized_alleles& alts) {
    Status s;

    // separate the ref and alt alleles
    discovered_alleles dalts;
    refs.clear();
    for (const auto& dal : src) {
        if (dal.second.is_ref) {
            const auto p = refs.find(dal.first.pos);
            if (p != refs.end() && dal.first.dna != p->second.first.dna) {
                return Status::Invalid("detected inconsistent REF alleles", dal.first.pos.str());
            }
            refs.insert(make_pair(dal.first.pos,dal));
        } else {
            dalts.insert(dal);
        }
    }

    // go through each alt allele
    alts.clear();
    for (const auto& dal : dalts) {
        const allele& alt = dal.first;

        // find corresponding reference allele
        const auto rp = refs.find(alt.pos);
        if (rp == refs.end()) return Status::Invalid("minimize_alleles: missing REF allele for ", alt.pos.str());

        // minimize the alt allele
        allele min_alt = alt;
        S(minimize_allele(rp->second.first, min_alt));

        unsigned copy_number = dal.second.zGQ.copy_number(cfg.min_GQ);

        // add it to alts, combining originals, copy_number, and topAQ with any
        // previously observed occurrences of the same minimized alt allele.
        minimized_allele_info info;
        info.originals.insert(dal.first);
        info.all_filtered = dal.second.all_filtered;
        info.topAQ = dal.second.topAQ;
        info.copy_number = copy_number;
        info.in_target = dal.second.in_target;
        auto ap = alts.find(min_alt);
        if (ap == alts.end()) {
            alts[min_alt] = info;
        } else {
            ap->second += info;
        }
    }

    return Status::OK();
}

// find the range in ranges overlapping pos, if any. (ranges assumed to be non-overlapping)
Status find_target_range(const std::set<range>& ranges, const range& pos, range& ans) {
    if (ranges.size() == 0) {
        return Status::NotFound();
    }

    // The returned value here is the first element that is
    // greater or equal to [pos].
    auto it = ranges.lower_bound(pos);
    if (it == ranges.end() || !it->overlaps(pos)) {
        // we landed one range after the one we need
        if (it != ranges.begin()) {
            it = std::prev(it);
        }
    }

    if (it->overlaps(pos)) {
        // we got the right range
        ans = *it;
        return Status::OK();
    }
    return Status::NotFound();
}


// Partition alleles into non-overlapping "active regions" where all the
// alleles in an active region are transitively connected through overlap
// or adjacency. (However, individual pairs of alleles within an active
// region might be separated.)
// The input is cleared by side-effect to save memory.
template<class discovered_or_minimized_alleles>
auto partition(discovered_or_minimized_alleles& alleles) {
    map<range,discovered_or_minimized_alleles> ans;

    range rng(-1,-1,-1);
    discovered_or_minimized_alleles als;

    for (auto pit = alleles.begin(); pit != alleles.end(); alleles.erase(pit++)) {
        const auto& it = *pit;
        assert(rng <= it.first.pos);
        if (rng.rid != it.first.pos.rid || rng.end < it.first.pos.beg) {
            if (rng.rid != -1) {
                assert(!als.empty());
                ans[rng] = als;
            }
            rng = it.first.pos;
            als.clear();
        }

        assert(rng.end >= it.first.pos.beg);
        assert(rng <= it.first.pos);
        rng.end = max(rng.end, it.first.pos.end);
        assert(it.first.pos.within(rng));
        als.insert(it);
    }

    if (rng.rid != -1) {
        assert(!als.empty());
        ans[rng] = als;
    }

    return ans;
}

bool check_filtered(const unifier_config& cfg, const minimized_allele& al) {
    return !(cfg.drop_filtered && al.second.all_filtered);
}

bool check_AQ(const unifier_config& cfg, const minimized_allele& al) {
    return (al.second.topAQ.V[0] >= cfg.min_AQ1 || al.second.topAQ.V[1] >= cfg.min_AQ2);
}

bool check_copy_number(const unifier_config& cfg, const minimized_allele& al) {
    return al.second.copy_number >= cfg.min_allele_copy_number;
}

// Given an active region, decompose it into "sites" by heuristically pruning
// rare or lengthy alleles to avoid excessive collapsing
auto prune_alleles(const unifier_config& cfg, const minimized_alleles& alleles, minimized_alleles& pruned) {
    vector<minimized_allele> valleles;
    valleles.reserve(alleles.size());
    pruned.clear();

    // filter alleles with insufficient copy number or AQ
    for (auto al : alleles) {
        // fix-up pass: ensure copy_number is >=1 for any allele with sufficient AQ.
        // this might not be the case up until this point, in the rare case when all
        // individuals carrying an allele have weak, homozygous-alt genotype calls
        if (cfg.min_AQ1 && al.second.topAQ.V[0] >= cfg.min_AQ1) {
            al.second.copy_number = std::max(al.second.copy_number, 1U);
        }
        if (cfg.min_AQ2 && al.second.topAQ.V[1] >= cfg.min_AQ2) {
            al.second.copy_number = std::max(al.second.copy_number, 2U);
        }

        if (check_filtered(cfg, al) && check_AQ(cfg, al) && check_copy_number(cfg, al)) {
            valleles.push_back(al);
        } else {
            pruned.insert(al);
        }
    }

    // sort the alt alleles by decreasing copy number (+ some tiebreakers)
    // or possibly by increasing range size, according to configuration
    sort(valleles.begin(), valleles.end(), minimized_allele_lt(cfg.preference));

    // Build the sites by considering each alt allele in that order, and
    // accepting or rejecting it as follows. If the allele overlaps exactly
    // one existing site, then accept it. Also accept the allele if it doesn't
    // overlap any already- accepted allele, creating a new site. Reject
    // alleles that overlap multiple existing sites.
    unsigned kept_allele_count = 0;
    map<range,minimized_alleles> sites;
    for (const minimized_allele& mal : valleles) {
        // find existing site(s) overlapping this allele
        auto related_site = sites.end();
        bool multiple_sites = false;
        for (auto site = sites.begin(); site != sites.end(); site++) {
            if (mal.first.pos.overlaps(site->first)) {
                if (related_site == sites.end()) {
                    related_site = site;
                } else {
                    multiple_sites = true;
                    break;
                }
            }
        }
        if (multiple_sites) {
            // reject since this allele overlaps multiple existing sites
            pruned.insert(mal);
            continue;
        }
        if (related_site != sites.end()) {
            // enforce max_alleles_per_site; NB, we count the ref allele toward
            // this limit, while related_site.size() is the # of ALT alleles.
            if (cfg.max_alleles_per_site > 1 && related_site->second.size()+1 >= cfg.max_alleles_per_site) {
                pruned.insert(mal);
                continue;
            }
            // merge this allele into the related site
            range urng(mal.first.pos.rid,
                       min(mal.first.pos.beg, related_site->first.beg),
                       max(mal.first.pos.end, related_site->first.end));
            minimized_alleles mals(move(related_site->second));
            mals.insert(mal);
            sites.erase(related_site);
            sites[urng] = move(mals);
        } else {
            // no related site: insert a new site
            sites[mal.first.pos] = minimized_alleles({mal});
        }
        kept_allele_count++;
    }
    assert(kept_allele_count+pruned.size() == alleles.size());

    return sites;
}

// Given reference alleles covering pos, reconstruct the reference allele for
// pos exactly.
Status unify_ref(const range& pos, const discovered_alleles& refs, allele& ref) {
    ref.pos = pos;
    ref.dna.assign(pos.size(), char(0));

    for (const auto& ref1 : refs) {
        UNPAIR(ref1, al, info);
        assert(pos.overlaps(al.pos));
        assert(info.is_ref);
        assert(al.pos.size() == al.dna.size());
        for (int i = 0, j = (al.pos.beg-pos.beg); i < al.dna.size(); i++, j++) {
            if (j >= 0 && j < ref.dna.size()) {
                if (ref.dna[j] && ref.dna[j] != al.dna[i]) {
                    return Status::Invalid("detected inconsistent REF alleles", al.pos.str());
                }
                ref.dna[j] = al.dna[i];
            }
        }
    }
    if (any_of(ref.dna.begin(), ref.dna.end(), [] (const char c) { return c == 0; })) {
        return Status::Invalid("Incomplete REF allele coverage in unification (this should not happen!)", pos.str());
    }

    return Status::OK();
}

// Given an ALT allele and a unified REF allele containing it, pad the ALT
// allele as needed with reference bases so that it covers exactly the
// same range as the ref allele.
Status pad_alt_allele(const allele& ref, allele& alt) {
    if (!ref.pos.contains(alt.pos)) return Status::Invalid("BUG: unified REF allele doesn't contain all ALT alleles");
    size_t l = alt.pos.beg - ref.pos.beg;
    size_t r = ref.pos.end - alt.pos.end;

    alt.pos = ref.pos;
    alt.dna = ref.dna.substr(0, l) + alt.dna + ref.dna.substr(ref.dna.size()-r, r);
    return Status::OK();
}

// Given minimized alleles in an active region, detect equivalent ALT alleles
// with different positions (i.e. indels with some non-left-aligned
// representation) and collapse them together. This effectively realigns the
// indels, but only in the limited case where they fall in the same active
// region as so far defined. A non-left-aligned indel could be distant from
// its left-aligned position e.g. in lengthy tandem repeats, and this will
// not deal with those.
Status unify_nonaligned_alts(const allele& ref, const minimized_alleles& alts,
                             minimized_alleles& aligned_alts) {
    Status s;
    map<allele,minimized_allele> unified_alts;
    for (const auto& p : alts) {
        UNPAIR(p, alt, alt_info);
        allele unified_alt(alt);
        S(pad_alt_allele(ref, unified_alt));
        assert(unified_alt.pos == ref.pos);

        auto q = unified_alts.find(unified_alt);
        if (q == unified_alts.end()) {
            unified_alts.insert(make_pair(unified_alt,p));
        } else {
            q->second.second += alt_info;
        }
    }
    aligned_alts.clear();
    for (const auto& p : unified_alts) {
        aligned_alts.insert(move(p.second));
    }
    return Status::OK();
}

// Delineate sites given all discovered alleles, potentially pruning some as
// described above. The result maps the range of each site to a tuple of the
// discovered REF alleles, the minimized ALT alleles, and any pruned alleles
// overlapping the site. It is possible (even likely) that one pruned allele
// overlaps multiple sites, and thus appears in more than one of those entries.
// all_pruned_alleles is a unique list of the pruned alleles, each with the
// corresponding reference allele.
// The input alleles is cleared by side-effect to save memory.
Status delineate_sites(const unifier_config& cfg, discovered_alleles& alleles,
                       map<range,tuple<discovered_alleles,minimized_alleles,minimized_alleles>>& ans,
                       vector<pair<minimized_allele,discovered_allele>>& all_pruned_alleles) {
    Status s;

    // Start by coarsely partitioning "active regions" (allele clusters)
    auto active_regions = partition<discovered_alleles>(alleles);
    // alleles has been cleared by side effect

    // Decompose each active region into sites
    ans.clear();
    all_pruned_alleles.clear();
    for (auto par = active_regions.begin(); par != active_regions.end(); active_regions.erase(par++)) {
        const auto& active_region = *par;

        // minimize the alt alleles
        map<range,discovered_allele> refs_by_range;
        minimized_alleles alts, pruned;
        S(minimize_alleles(cfg, active_region.second, refs_by_range, alts));

        // exclude alleles not overlapping the discovery target range, if any,
        // after minimization.
        for (auto it = alts.begin(); it != alts.end(); ) {
            auto trg = it->second.in_target;
            if (trg.rid < 0 || trg.overlaps(it->first.pos)) {
                ++it;
            } else {
                alts.erase(it++);
            }
        }

        // reconstruct the active region's reference allele
        discovered_alleles refs;
        for (const auto& p : refs_by_range) {
            refs.insert(p.second);
        }
        allele active_region_ref(active_region.first,"A");
        S(unify_ref(active_region.first, refs, active_region_ref));

        // detect equivalent alt alleles at different positions, and collapse them
        minimized_alleles aligned_alts;
        S(unify_nonaligned_alts(active_region_ref, alts, aligned_alts));
        alts = move(aligned_alts);

        // prune alt alleles as necessary to yield sites
        const auto sites = prune_alleles(cfg, alts, pruned);

        for (const auto& site : sites) {
            // find the ref alleles overlapping this site
            discovered_alleles site_refs;
            for (const auto& ref : refs) {
                if (ref.first.pos.overlaps(site.first)) {
                    site_refs.insert(ref);
                }
            }

            // and the pruned alleles
            minimized_alleles site_pruned;
            for (const auto& p : pruned) {
                if (p.first.pos.overlaps(site.first)) {
                    site_pruned.insert(p);
                }
            }

            assert(ans.find(site.first) == ans.end());
            ans[site.first] = make_tuple(site_refs,site.second,site_pruned);
        }

        for (const auto& pa : pruned) {
            // we need to supply a discovered reference allele to go along with
            // the pruned alt allele; find the shortest one which contains the
            // alt. note, we realigned the alt so this might not cover the
            // original!
            const discovered_allele *shortest_containing_ref = nullptr;
            for (const auto& ref : refs_by_range) {
                if (ref.first.contains(pa.first.pos) &&
                    (!shortest_containing_ref || ref.first.size() < shortest_containing_ref->first.pos.size())) {
                        shortest_containing_ref = &ref.second;
                }
            }
            if (!shortest_containing_ref) {
                return Status::Invalid("delineate_sites: missing REF allele for ", pa.first.str());
            }
            all_pruned_alleles.push_back(make_pair(pa, *shortest_containing_ref));
        }
    }
    return Status::OK();
}

/// Unify the alleles at one site.
Status unify_alleles(const unifier_config& cfg, unsigned N, const range& pos,
                     const discovered_alleles& refs, const minimized_alleles& alts,
                     const minimized_alleles& pruned, unified_site& ans) {
    Status s;

    // collapse the refs to get the reference allele for this site
    allele ref(pos, "A");
    S(unify_ref(pos, refs, ref));

    // For presentation in the unified site, sort the alt alleles by
    // decreasing copy number count (+ some tiebreakers). Note, this may be a
    // different sort order than used in prune_allele earlier.
    assert(alts.size() > 0);
    assert(cfg.max_alleles_per_site < 2 || alts.size() < cfg.max_alleles_per_site);
    vector<minimized_allele> valts(alts.begin(), alts.end());
    sort(valts.begin(), valts.end(), minimized_allele_common_lt);

    // fill out the unification
    unified_site us(pos);
    // We don't attempt to specify the ref allele frequency for now; we don't quite
    // have enough information because the copy numbers in discovered_alleles omit
    // nearly all homozygous ref genotypes.
    us.alleles.push_back(unified_allele(ref.pos, ref.dna));
    for (int i = 1; i <= valts.size(); i++) {
        const minimized_allele& p = valts[i-1];
        UNPAIR(p, alt, alt_info);

        allele unified_alt(alt);
        S(pad_alt_allele(ref, unified_alt));
        assert(unified_alt.pos == ref.pos);
        unified_allele ua(unified_alt.pos, unified_alt.dna);
        ua.normalized = alt;
        ua.quality = alt_info.topAQ.V[0];
        float freq = float(alt_info.copy_number)/(N*zygosity_by_GQ::PLOIDY);
        ua.frequency = roundf(std::max(1.0f,freq*1e6f))/1e6f;
        us.alleles.push_back(ua);

        for (const auto& original : alt_info.originals) {
            us.unification[original] = i;
        }
    }
    us.fill_implicit_unification();

    // Sum frequency of pruned alleles. In general this may be an overestimate when
    // some of those alleles co-occur in cis; we accept this possible error for now.
    // TODO: consider collecting some LD information during allele discovery.
    // Furthermore, if we have any lost alleles, we floor at 1/(2N) regardless of
    // min_GQ.
    float lost_allele_frequency = 0.0;
    for (const auto& p : pruned) {
        lost_allele_frequency += p.second.copy_number;
    }
    lost_allele_frequency /= N*zygosity_by_GQ::PLOIDY;
    lost_allele_frequency = ceilf(lost_allele_frequency*1e6f)/1e6f;
    lost_allele_frequency = std::min(lost_allele_frequency, 1.0f);
    us.lost_allele_frequency = lost_allele_frequency;

    // set variant QUAL score to maximum AQ across all ALT alleles
    // this is a sort of conservative bound on "prob(no variant)"
    // e.g. evidence from multiple samples doesn't stack
    us.qual = 0;
    for (const auto& p : alts) {
        us.qual = std::max(us.qual, p.second.topAQ.V[0]);
    }

    // We expect the in_target of allthe site's alleles to be the same, but anyway
    // look for the largest.
    for (const auto& p : alts) {
        if (p.second.in_target.size() > us.in_target.size()) {
            us.in_target = p.second.in_target;
        }
    }

    ans = std::move(us);
    return Status::OK();
}

Status unified_sites(const unifier_config& cfg,
                     unsigned N, discovered_alleles& alleles,
                     vector<unified_site>& ans,
                     unifier_stats& stats_out) {
    Status s;
    unifier_stats stats;

    /* desperate-straits debugging:
    for (const auto& allele : alleles) {
        if (allele.first.pos.overlaps(range(7, 16188960, 16188970))) {
            cerr << allele.first.str();
            if (allele.second.is_ref) {
                cerr << " *";
            }
            cerr << endl;
        }
    }
    */

    map<range,tuple<discovered_alleles,minimized_alleles,minimized_alleles>> sites;
    vector<pair<minimized_allele,discovered_allele>> all_pruned_alleles;
    S(delineate_sites(cfg, alleles, sites, all_pruned_alleles));
    // at this point, alleles has been cleared to save memory usage

    for (auto psite = sites.begin(); psite != sites.end(); sites.erase(psite++)) {
        UNPAIR(*psite, pos, site_alleles);
        const auto& ref_alleles = get<0>(site_alleles);
        const auto& alt_alleles = get<1>(site_alleles);
        auto pruned_alleles = get<2>(site_alleles);
        unified_site us(pos);
        S(unify_alleles(cfg, N, pos, ref_alleles, alt_alleles, pruned_alleles, us));
        ans.push_back(us);
        stats.unified_alleles += alt_alleles.size();
    }

    auto k = ans.size();
    for (const auto& pa : all_pruned_alleles) {
        if (check_filtered(cfg, pa.first) && check_AQ(cfg, pa.first) && check_copy_number(cfg, pa.first)) {
            if (cfg.monoallelic_sites_for_lost_alleles) {
                unified_site ms(pa.first.first.pos);
                S(unify_alleles(cfg, N, pa.first.first.pos, discovered_alleles{pa.second},
                                minimized_alleles{pa.first}, minimized_alleles(), ms));
                ms.monoallelic = true;
                ms.in_target = pa.first.second.in_target;
                ans.push_back(ms);
            }
            stats.lost_alleles++;
        } else {
            stats.filtered_alleles++;
        }
    }
    // merge any the newly added monoallelic sites in position order with the others
    // (in linear time)
    std::inplace_merge(ans.begin(), ans.begin()+k, ans.end());
    assert(std::is_sorted(ans.begin(), ans.end()));

    stats_out = stats;
    return Status::OK();
}

}

#include <algorithm>
#include <assert.h>
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
    top_AQ topAQ;
    unsigned copy_number = 0;

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
                        discovered_alleles& refs, minimized_alleles& alts) {
    Status s;

    // separate the ref and alt alleles
    map<range,discovered_allele> refs_by_range;
    discovered_alleles dalts;
    refs.clear();
    for (const auto& dal : src) {
        if (dal.second.is_ref) {
            refs.insert(dal);
            refs_by_range.insert(make_pair(dal.first.pos, dal));
        } else {
            dalts.insert(dal);
        }
    }

    // go through each alt allele
    alts.clear();
    for (const auto& dal : dalts) {
        const allele& alt = dal.first;

        // find corresponding reference allele
        const auto rp = refs_by_range.find(alt.pos);
        if (rp == refs_by_range.end()) return Status::Invalid("minimize_alleles: missing REF allele for ", alt.pos.str());

        // minimize the alt allele
        allele min_alt = alt;
        S(minimize_allele(rp->second.first, min_alt));

        unsigned copy_number = dal.second.zGQ.copy_number(cfg.min_GQ);

        // add it to alts, combining originals, copy_number, and topAQ with any
        // previously observed occurrences of the same minimized alt allele.
        auto ap = alts.find(min_alt);
        if (ap == alts.end()) {
            minimized_allele_info info;
            info.originals.insert(dal.first);
            info.topAQ = dal.second.topAQ;
            info.copy_number = copy_number;
            alts[min_alt] = move(info);
        } else {
            ap->second.originals.insert(dal.first);
            ap->second.topAQ += dal.second.topAQ;
            ap->second.copy_number += copy_number;
        }
    }

    return Status::OK();
}


// Partition alleles into non-overlapping "active regions" where all the
// alleles in an active region are transitively connected through overlap.
// (However, individual pairs of alleles within an active region might be non-
// overlapping.)
// The input is cleared by side-effect to save memory.
template<class discovered_or_minimized_alleles>
auto partition(discovered_or_minimized_alleles& alleles) {
    map<range,discovered_or_minimized_alleles> ans;

    range rng(-1,-1,-1);
    discovered_or_minimized_alleles als;

    for (auto pit = alleles.begin(); pit != alleles.end(); alleles.erase(pit++)) {
        const auto& it = *pit;
        if (!rng.overlaps(it.first.pos)) {
            if (rng.rid != -1) {
                assert(!als.empty());
                ans[rng] = als;
            }
            rng = it.first.pos;
            als.clear();
        }

        assert(rng.overlaps(it.first.pos));
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

// Given a cluster of related/overlapping ALT alleles, decompose it into
// "sites" by heuristically pruning rare or lengthy alleles to avoid excessive
// collapsing
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

        if ((al.second.topAQ.V[0] >= cfg.min_AQ1 || al.second.topAQ.V[1] >= cfg.min_AQ2) &&
            al.second.copy_number >= cfg.min_allele_copy_number) {
            valleles.push_back(al);
        } else {
            pruned.insert(al);
        }
    }

    // sort the alt alleles by decreasing copy number (+ some tiebreakers)
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

// Delineate sites given all discovered alleles, potentially pruning some as
// described above. The result maps the range of each site to a pair of the
// discovered REF alleles and the minimized ALT alleles.
// The input alleles is cleared by side-effect to save memory.
Status delineate_sites(const unifier_config& cfg, discovered_alleles& alleles,
                       map<range,tuple<discovered_alleles,minimized_alleles,minimized_alleles>>& ans) {
    Status s;

    // Start by coarsely partitioning "active regions" (allele clusters)
    auto active_regions = partition<discovered_alleles>(alleles);
    // alleles has been cleared by side effect

    // Decompose each active region into sites
    ans.clear();
    for (auto par = active_regions.begin(); par != active_regions.end(); active_regions.erase(par++)) {
        const auto& active_region = *par;

        // minimize the alt alleles
        discovered_alleles refs;
        minimized_alleles alts, pruned;
        S(minimize_alleles(cfg, active_region.second, refs, alts));

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
    }
    return Status::OK();
}

// Given the alleles at a site, collapse the discovered REF alleles to get one
// reference allele covering the site's range.
Status unify_ref(const range& pos, const discovered_alleles& refs, allele& ref) {
    ref.pos = pos;
    ref.dna.assign(pos.size(), char(0));

    for (const auto& ref1 : refs) {
        UNPAIR(ref1, al, info);
        assert(pos.overlaps(al.pos));
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

// Given the unified REF allele and an ALT allele in the same site, pad the
// ALT allele with reference bases so that it covers exactly the site's range.
Status pad_alt_allele(const allele& ref, allele& alt) {
    if (!ref.pos.contains(alt.pos)) return Status::Invalid("BUG: unified REF allele doesn't contain all ALT alleles");
    size_t l = alt.pos.beg - ref.pos.beg;
    size_t r = ref.pos.end - alt.pos.end;

    alt.pos = ref.pos;
    alt.dna = ref.dna.substr(0, l) + alt.dna + ref.dna.substr(ref.dna.size()-r, r);
    return Status::OK();
}

/// Unify the alleles at one site.
/// pruned is in+out; we need to know which alleles were pruned in site delineation, and,
/// we might prune more alleles in the process.
Status unify_alleles(const unifier_config& cfg, unsigned N, const range& pos,
                     const discovered_alleles& refs, const minimized_alleles& alts,
                     minimized_alleles& pruned, unified_site& ans) {
    Status s;

    // collapse the refs to get the reference allele for this site
    allele ref(pos, "A");
    S(unify_ref(pos, refs, ref));

    // For presentation in the unified site, sort the alt alleles by
    // decreasing copy number count (+ some tiebreakers). Note, this may be a
    // different sort order than used in prune_allele earlier.
    assert(alts.size() > 0);
    vector<minimized_allele> valts(alts.begin(), alts.end());
    sort(valts.begin(), valts.end(), minimized_allele_common_lt);

    // enforce max_alleles_per_site. NB, valts does not include the ref
    // allele, so we're truncating valts to cfg.max_alleles_per_site-1
    if (cfg.max_alleles_per_site > 1 && valts.size() >= cfg.max_alleles_per_site) {
        // Out of fairness, we also prune alt alleles with the same
        // copy number as the first truncated one, except we always need
        // to keep at least one alt allele of course.
        unsigned int trunc0 = 1, trunc1 = 1;
        for (; trunc1 < cfg.max_alleles_per_site; trunc1++) {
            assert(valts[trunc0].second.copy_number >= valts[trunc1].second.copy_number);
            while(valts[trunc0].second.copy_number > valts[trunc1].second.copy_number)
                trunc0++;
        }
        pruned.insert(valts.begin()+trunc0, valts.end());
        valts.assign(valts.begin(), valts.begin()+trunc0);
        assert(valts.size() < cfg.max_alleles_per_site);
    }

    // fill out the unification
    unified_site us(pos);
    us.alleles.push_back(ref.dna);
    // We don't attempt to specify the ref allele frequency for now; we don't quite
    // have enough information because the copy numbers in discovered_alleles omit
    // nearly all homozygous ref genotypes.
    us.allele_frequencies.push_back(NAN);
    for (int i = 1; i <= valts.size(); i++) {
        const minimized_allele& p = valts[i-1];
        UNPAIR(p, alt, alt_info);

        allele unified_alt(alt);
        S(pad_alt_allele(ref, unified_alt));
        assert(unified_alt.pos == ref.pos);

        us.alleles.push_back(unified_alt.dna);
        // TODO: configurable ploidy setting
        float freq = float(alt_info.copy_number)/(N*zygosity_by_GQ::PLOIDY);
        us.allele_frequencies.push_back(ceilf(freq*1e6f)/1e6f);
        for (const auto& original : alt_info.originals) {
            us.unification[original] = i;
        }
    }

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

    ans = std::move(us);
    return Status::OK();
}

Status unified_sites(const unifier_config& cfg,
                     unsigned N, discovered_alleles& alleles,
                     vector<unified_site>& ans) {
    Status s;

    map<range,tuple<discovered_alleles,minimized_alleles,minimized_alleles>> sites;
    S(delineate_sites(cfg, alleles, sites));
    // at this point, alleles has been cleared to save memory usage

    for (auto psite = sites.begin(); psite != sites.end(); sites.erase(psite++)) {
        UNPAIR(*psite, pos, site_alleles);
        const auto& ref_alleles = get<0>(site_alleles);
        const auto& alt_alleles = get<1>(site_alleles);
        auto pruned_alleles = get<2>(site_alleles);
        unified_site us(pos);
        S(unify_alleles(cfg, N, pos, ref_alleles, alt_alleles, pruned_alleles, us));
        ans.push_back(us);
    }
    return Status::OK();
}

}

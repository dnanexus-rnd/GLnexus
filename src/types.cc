#include "types.h"
#include <algorithm>
#include <regex>

// cap'n proto headers
#include <capnp/message.h>
#include <capnp/serialize-packed.h>
#include <defs.capnp.h>

// For file descriptors
//#include <type_traits>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

using namespace std;

namespace GLnexus {

// We assume the chromosome ploidy is two in our serialization format.
static_assert(zygosity_by_GQ::PLOIDY == 2, "PLOIDY needs to be two");

// This limit should be much larger than the message size. It is a security,
// preventing the demarshaller from using too much memory, and getting into
// infinite loops.
const uint64_t CAPNP_TRAVERSAL_LIMIT = 20 * 1024L * 1024L * 1024L;

regex regex_id ("[-_a-zA-Z0-9\\.]{1,100}");

// We'd prefer is_dna and is_iupac_nucleotides be done with regex as originally.
// However we ran into some appalling libstdc++ regex bugs with long strings
//   https://gcc.gnu.org/bugzilla/show_bug.cgi?id=61582
// So we replaced with these hard-coded functions, which themselves could be
// improved to use lookup tables.
bool is_dna(const std::string& s) {
    size_t sz = s.size();
    if (sz == 0) {
        return false;
    }
    for (size_t i = 0; i < sz; i++) {
        switch (s[i]) {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
                continue;
            default:
                return false;
        }
    }
    return true;
}
bool is_iupac_nucleotides(const std::string& s) {
    size_t sz = s.size();
    if (sz == 0) {
        return false;
    }
    for (size_t i = 0; i < sz; i++) {
        switch (s[i]) {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
            case 'U':
            case 'R':
            case 'Y':
            case 'S':
            case 'W':
            case 'K':
            case 'M':
            case 'B':
            case 'D':
            case 'H':
            case 'V':
            case 'N':
                continue;
            default:
                return false;
        }
    }
    return true;
}

// Add src alleles to dest alleles. Identical alleles alleles are merged,
// updating topAQ and combining zygosity_by_GQ
Status merge_discovered_alleles(const discovered_alleles& src, discovered_alleles& dest) {
    for (auto& dsal : src) {
        UNPAIR(dsal,allele,ai)
        auto p = dest.find(allele);
        if (p == dest.end()) {
            dest[allele] = ai;
        } else {
            if (ai.is_ref != p->second.is_ref) {
                return Status::Invalid("allele appears as both REF and ALT", allele.dna + "@" + allele.pos.str());
            }
            p->second.all_filtered = p->second.all_filtered && ai.all_filtered;
            p->second.topAQ += ai.topAQ;
            p->second.zGQ += ai.zGQ;
        }
    }

    return Status::OK();
}


Status range_yaml(const std::vector<std::pair<std::string,size_t> >& contigs,
                  const range& r, YAML::Emitter& yaml, bool omit_ref) {
    if (!omit_ref && (r.rid < 0 || r.rid >= contigs.size())) {
        return Status::NotFound("range_yaml: invalid rid", r.str());
    }
    yaml << YAML::Flow;
    yaml << YAML::BeginMap;
    if (!omit_ref) {
        yaml << YAML::Key << "ref" << YAML::Value << contigs[r.rid].first;
    }
    yaml << YAML::Key << "beg" << YAML::Value << r.beg+1;
    yaml << YAML::Key << "end" << YAML::Value << r.end;
    yaml << YAML::EndMap;
    return Status::OK();
}

Status range_of_yaml(const YAML::Node& yaml, const vector<pair<string,size_t> >& contigs,
                     range& ans, int default_rid) {
    #define V(pred,msg) if (!(pred)) return Status::Invalid("range_of_yaml: " msg)

    V(yaml.IsMap(), "not a map at top level");

    int rid;
    const auto n_ref = yaml["ref"];
    if (n_ref) {
        V(n_ref.IsScalar(), "invalid 'ref' field");
        const string& ref = n_ref.Scalar();
        for (rid=0; rid < contigs.size(); rid++) {
            if (ref == contigs[rid].first) {
                break;
            }
        }
        if (rid == contigs.size()) {
            return Status::Invalid("range_of_yaml: unknown contig", ref);
        }
    } else if (default_rid >= 0) {
        rid = default_rid;
    } else {
        return Status::Invalid("range_of_yaml: missing 'ref' field");
    }

    const auto n_beg = yaml["beg"];
    V(n_beg && n_beg.IsScalar(), "missing/invalid 'beg' field");
    int beg = n_beg.as<int>()-1;

    const auto n_end = yaml["end"];
    V(n_end && n_end.IsScalar(), "missing/invalid 'end' field");
    int end = n_end.as<int>();

    V(beg >= 0 && end >= beg, "invalid beg/end coordinates");

    ans = range(rid, beg, end);
    return Status::OK();
    #undef V
}

Status yaml_of_one_discovered_allele(const allele& allele,
                                     const discovered_allele_info& ainfo,
                                     const std::vector<std::pair<std::string,size_t> >& contigs,
                                     YAML::Emitter& out) {
    Status s;
    out << YAML::BeginMap;

    out << YAML::Key << "range" << YAML::Value;
    S(range_yaml(contigs, allele.pos, out));
    out << YAML::Key << "dna"
        << YAML::Value << allele.dna;
    out << YAML::Key << "is_ref"
        << YAML::Value << ainfo.is_ref;
    out << YAML::Key << "all_filtered"
        << YAML::Value << ainfo.all_filtered;
    out << YAML::Key << "top_AQ"
        << YAML::Value << YAML::Flow << YAML::BeginSeq;
    for (unsigned i = 0; i < top_AQ::COUNT; i++) {
        auto x = ainfo.topAQ.V[i];
        if (i > 0 && x <= 0) {
            break;
        }
        out << x;
    }
    out << YAML::EndSeq;

    out << YAML::Key << "zygosity_by_GQ"
        << YAML::Value << YAML::Flow << YAML::BeginSeq;
    for (unsigned i = 0; i < zygosity_by_GQ::GQ_BANDS; i++) {
        out << YAML::Flow << YAML::BeginSeq;
        for (unsigned j = 0; j < zygosity_by_GQ::PLOIDY; j++) {
            out << ainfo.zGQ.M[i][j];
        }
        out << YAML::EndSeq;
    }
    out << YAML::EndSeq;

    out << YAML::EndMap;
    return Status::OK();
}

Status yaml_of_discovered_alleles(const discovered_alleles& dals,
                                  const std::vector<std::pair<std::string,size_t> >& contigs,
                                  YAML::Emitter& yaml) {
    Status s;

    yaml << YAML::BeginSeq;
    for (const auto& p : dals) {
        yaml_of_one_discovered_allele(p.first, p.second, contigs, yaml);
    }
    yaml << YAML::EndSeq;

    return Status::OK();
}


Status one_discovered_allele_of_yaml(const YAML::Node& yaml,
                                     const std::vector<std::pair<std::string,size_t> >& contigs,
                                     allele& dsal,
                                     discovered_allele_info &ainfo) {
    Status s;
    #define V(pred,msg) if (!(pred)) return Status::Invalid("discovered_alleles_of_yaml: " msg)

    V(yaml.IsMap(), "invalid entry");

    range rng(-1,-1,-1);
    const auto n_range = yaml["range"];
    V(n_range, "missing 'range' field in entry");
    S(range_of_yaml(n_range, contigs, rng));
    #define VR(pred,msg) if (!(pred)) return Status::Invalid("discovered_alleles_of_yaml: " msg, rng.str(contigs))

    const auto n_dna = yaml["dna"];
    VR(n_dna && n_dna.IsScalar(), "missing/invalid 'dna' field in entry");
    const string& dna = n_dna.Scalar();
    VR(dna.size() > 0, "empty 'dna' in entry");
    VR(is_iupac_nucleotides(dna), "invalid allele DNA");
    allele al(rng, dna);
    dsal = std::move(al);

    discovered_allele_info ai;
    const auto n_is_ref = yaml["is_ref"];
    VR(n_is_ref && n_is_ref.IsScalar(), "missing/invalid 'is_ref' field in entry");
    ai.is_ref = n_is_ref.as<bool>();

    const auto n_all_filtered = yaml["all_filtered"];
    VR(n_all_filtered && n_all_filtered.IsScalar(), "missing/invalid 'all_filtered' field in entry");
    ai.all_filtered = n_all_filtered.as<bool>();

    const auto n_topAQ = yaml["top_AQ"];
    VR(n_topAQ && n_topAQ.IsSequence(), "missing/invalid 'top_AQ' field in entry");
    VR(n_topAQ.size() <= top_AQ::COUNT, "unexpected top_AQ size");
    unsigned i = 0;
    for (YAML::const_iterator aq = n_topAQ.begin(); aq != n_topAQ.end(); ++aq, ++i) {
        VR(aq->IsScalar() && aq->as<int>() >= 0, "invalid entry in top_AQ");
        ai.topAQ.V[i] = aq->as<int>();
    }

    const auto n_zygosity_by_GQ = yaml["zygosity_by_GQ"];
    VR(n_zygosity_by_GQ && n_zygosity_by_GQ.IsSequence(), "missing/invalid 'zygosity_by_GQ' field in entry");
    VR(n_zygosity_by_GQ.size() == zygosity_by_GQ::GQ_BANDS, "unexpected row count in zygosity_by_GQ");
    i = 0;
    for (YAML::const_iterator q = n_zygosity_by_GQ.begin(); q != n_zygosity_by_GQ.end(); ++q, ++i) {
        VR(q->IsSequence(), "invalid row in zygosity_by_GQ");
        VR(q->size() == zygosity_by_GQ::PLOIDY, "unexpected row size in zygosity_by_GQ");
        unsigned j = 0;
        for (YAML::const_iterator r = q->begin(); r != q->end(); ++r, ++j) {
            VR(r->IsScalar() && r->as<int>() >= 0, "invalid entry in zygosity_by_GQ");
            ai.zGQ.M[i][j] = (unsigned) r->as<int>();
        }
    }
    #undef VR
    #undef V

    ainfo = std::move(ai);
    return Status::OK();
}

Status discovered_alleles_of_yaml(const YAML::Node& yaml,
                                  const std::vector<std::pair<std::string,size_t> >& contigs,
                                  discovered_alleles& ans) {
    Status s;
    #define V(pred,msg) if (!(pred)) return Status::Invalid("discovered_alleles_of_yaml: " msg)

    V(yaml.IsSequence(), "not a sequence at top level");
    ans.clear();
    for (YAML::const_iterator p = yaml.begin(); p != yaml.end(); ++p) {
        V(p->IsMap(), "invalid entry");

        allele dsal(range(-1,-1,-1), "A");
        discovered_allele_info ainfo;
        S(one_discovered_allele_of_yaml((*p), contigs, dsal, ainfo));

        V(ans.find(dsal) == ans.end(), "duplicate alleles");
        ans[dsal] = ainfo;
    }

    #undef V
    return Status::OK();
}

// write discovered_alleles structure to a file-descriptor, with cap'n proto serialization
static Status capnp_write_discovered_alleles_fd(unsigned int sample_count,
                                                const std::vector<std::pair<std::string,size_t> >& contigs,
                                                const discovered_alleles &dsals,
                                                int fd) {
    ::capnp::MallocMessageBuilder message;
    capnp::DiscoveredAlleles::Builder dsals_pk = message.initRoot<capnp::DiscoveredAlleles>();
    dsals_pk.setSampleCount(sample_count);

    // serialize contigs
    ::capnp::List<capnp::Contig>::Builder contigs_pk = dsals_pk.initContigs(contigs.size());
    int cursor = 0;
    for (auto const &kv : contigs) {
        const string &name = kv.first;
        size_t size = kv.second;

        capnp::Contig::Builder ctg_pk = contigs_pk[cursor];
        ctg_pk.setName(name);
        ctg_pk.setSize(size);
        cursor++;
    }

    // serialize allele-info pairs
    ::capnp::List<capnp::AlleleInfoPair>::Builder aips_pk = dsals_pk.initAips(dsals.size());
    cursor = 0;
    for (auto const &kv : dsals) {
        auto &key = kv.first;
        auto &val = kv.second;

        capnp::AlleleInfoPair::Builder al_pk = aips_pk[cursor];
        capnp::Range::Builder range = al_pk.initRange();
        range.setRid(key.pos.rid);
        range.setBeg(key.pos.beg);
        range.setEnd(key.pos.end);
        al_pk.setDna(key.dna.c_str());

        capnp::DiscoveredAlleleInfo::Builder dai = al_pk.initDai();
        dai.setIsRef(val.is_ref);
        dai.setAllFiltered(val.all_filtered);

        ::capnp::List<int64_t>::Builder topAQ = dai.initTopAQ(top_AQ::COUNT);
        for (int k=0; k < top_AQ::COUNT; k++) {
            topAQ.set(k, val.topAQ.V[k]);
        }

        // skipping addbuf, it is used for scratch space

        ::capnp::List<uint64_t>::Builder zGQ0 = dai.initZGQ0(zygosity_by_GQ::GQ_BANDS);
        for (int k=0; k < zygosity_by_GQ::GQ_BANDS; k++)
            zGQ0.set(k, val.zGQ.M[k][0]);

        ::capnp::List<uint64_t>::Builder zGQ1 = dai.initZGQ1(zygosity_by_GQ::GQ_BANDS);
        for (int k=0; k < zygosity_by_GQ::GQ_BANDS; k++)
            zGQ1.set(k, val.zGQ.M[k][1]);

        cursor++;
    }

    // Capnp throws exceptions on errors, and only in extreme cases.
    writePackedMessageToFd(fd, message);
    return Status::OK();
}

// Variant of the above command, but write to a file
Status capnp_of_discovered_alleles_fd(unsigned int sample_count,
                                      const std::vector<std::pair<std::string,size_t> >& contigs,
                                      const discovered_alleles &dsals,
                                      int fd) {
    try {
        return capnp_write_discovered_alleles_fd(sample_count, contigs, dsals, fd);
    } catch (exception &e) {
        return Status::IOError("Capnproto error during de-serialization", e.what());
    }
}

// Variant of the above command, but write to a file
Status capnp_of_discovered_alleles(unsigned int sample_count,
                                   const std::vector<std::pair<std::string,size_t> >& contigs,
                                   const discovered_alleles &dsals,
                                   const std::string &filename) {
    std::remove(filename.c_str());
    int fd = open(filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRWXU);
    if (fd < 0)
        return Status::IOError("Could not open file for writing", filename);
    Status s = capnp_of_discovered_alleles_fd(sample_count, contigs, dsals, fd);
    int retval = close(fd);
    if (retval < 0)
        return Status::IOError("file close failed", filename);
    return s;
}

// read discovered_alleles structure from a file-descriptor, as serialized by cap'n proto
static Status _capnp_read_discovered_alleles_fd(int fd,
                                                unsigned int &sample_count,
                                                std::vector<std::pair<std::string,size_t> >& contigs,
                                                discovered_alleles &dsals) {
    ::capnp::ReaderOptions ropt;
    ropt.traversalLimitInWords = CAPNP_TRAVERSAL_LIMIT;
    ::capnp::PackedFdMessageReader message(fd, ropt);
    capnp::DiscoveredAlleles::Reader dsals_pk = message.getRoot<capnp::DiscoveredAlleles>();
    sample_count = dsals_pk.getSampleCount();

    // contigs
    contigs.clear();
    for (capnp::Contig::Reader ctg : dsals_pk.getContigs()) {
        auto p = make_pair<string,size_t>(ctg.getName(), ctg.getSize());
        contigs.push_back(p);
    }

    // Discovered Alleles
    dsals.clear();
    for (capnp::AlleleInfoPair::Reader aip : dsals_pk.getAips()) {
        // Unpack allele
        capnp::Range::Reader range_pk = aip.getRange();
        range r(range_pk.getRid(),
                range_pk.getBeg(),
                range_pk.getEnd());
        string dna = aip.getDna();
        allele alle(r, dna);

        // Unpack allele-info
        discovered_allele_info dai;
        capnp::DiscoveredAlleleInfo::Reader dai_pk = aip.getDai();
        dai.is_ref = dai_pk.getIsRef();
        dai.all_filtered = dai_pk.getAllFiltered();

        ::capnp::List<int64_t>::Reader topAQ_pk = dai_pk.getTopAQ();
        if (topAQ_pk.size() != top_AQ::COUNT) {
            return Status::Invalid("Wrong number of elements in topAQ", to_string(topAQ_pk.size()));
        }
        for (int k=0; k < topAQ_pk.size(); k++)  {
            dai.topAQ.V[k] =  topAQ_pk[k];
        }
        // skipping addbuf, it is used for scratch space

        // Unpack Zygosity information
        ::capnp::List<uint64_t>::Reader zGQ0_pk = dai_pk.getZGQ0();
        ::capnp::List<uint64_t>::Reader zGQ1_pk = dai_pk.getZGQ1();
        if (zGQ0_pk.size() != zygosity_by_GQ::GQ_BANDS)
            return Status::Invalid("Wrong number of GQ_BANDS", to_string(zGQ0_pk.size()));
        if (zGQ1_pk.size() != zygosity_by_GQ::GQ_BANDS)
            return Status::Invalid("Wrong number of GQ_BANDS", to_string(zGQ1_pk.size()));
        for (int k=0; k < zygosity_by_GQ::GQ_BANDS; k++) {
            dai.zGQ.M[k][0] = zGQ0_pk[k];
            dai.zGQ.M[k][1] = zGQ1_pk[k];
        }

        dsals[alle] = dai;
    }

    return Status::OK();
}

Status discovered_alleles_of_capnp(const std::string &filename,
                                   unsigned int &sample_count,
                                   std::vector<std::pair<std::string,size_t> >& contigs,
                                   discovered_alleles &dsals) {
    int fd = open(filename.c_str(), O_RDONLY);
    if (fd < 0)
        return Status::IOError("Could not open file for reading", filename);
    Status s;
    try {
        // Capnp throws exceptions on errors, and only in extreme cases.
        s = _capnp_read_discovered_alleles_fd(fd, sample_count, contigs, dsals);
    } catch (exception &e) {
        return Status::IOError("Capnproto error during de-serialization", e.what());
    }

    int retval = close(fd);
    if (retval < 0)
        return Status::IOError("file close failed", filename);
    return s;
}


Status capnp_discover_alleles_verify(unsigned int sample_count,
                                     const std::vector<std::pair<std::string,size_t> >& contigs,
                                     const discovered_alleles &dsals,
                                     const string &filename) {
    Status s;

    // Write to file
    S(capnp_of_discovered_alleles(sample_count, contigs, dsals, filename));

    // read from it
    unsigned int N;
    vector<pair<string,size_t>> contigs2;
    discovered_alleles dsals2;
    S(discovered_alleles_of_capnp(filename, N, contigs2, dsals2));

    // verify we get the same alleles back
    if (dsals == dsals2 &&
        contigs == contigs2 &&
        sample_count == N) {
        return Status::OK();
    }

    return Status::Invalid("capnp serialization/deserialization of discovered alleles does not return original value");
}


Status unified_site::yaml(const std::vector<std::pair<std::string,size_t> >& contigs,
                          YAML::Emitter& ans) const {
    Status s;

    ans << YAML::BeginMap;

    ans << YAML::Key << "range" << YAML::Value;
    S(range_yaml(contigs, pos, ans));

    if (in_target.rid >= 0) {
        ans << YAML::Key << "in_target" << YAML::Value;
        S(range_yaml(contigs, in_target, ans));
    }

    ans << YAML::Key << "alleles";
    ans << YAML::Value << YAML::Flow << YAML::BeginSeq;
    for (const auto& allele : alleles) {
        ans << allele;
    }
    ans << YAML::EndSeq;

    ans << YAML::Key << "allele_frequencies";
    ans << YAML::Value << YAML::Flow << YAML::BeginSeq;
    for (auto f : allele_frequencies) {
        if (f == f)
            ans << f;
        else
            // working around yaml-cpp issue: it emits "nan", not ".nan"
            ans << ".nan";
    }
    ans << YAML::EndSeq;

    if (lost_allele_frequency != 0.0f) {
        ans << YAML::Key << "lost_allele_frequency";
        ans << YAML::Value << lost_allele_frequency;
    }

    ans << YAML::Key << "quality" << YAML::Value << qual;

    if (monoallelic) {
        ans << YAML::Key << "monoallelic" << YAML::Value << true;
    }

    ans << YAML::Key << "unification";
    ans << YAML::Value << YAML::BeginSeq;
    // sort unification entries by to, then by range, then by dna
    vector<pair<allele,int>> u(unification.begin(), unification.end());
    sort(u.begin(), u.end(),
         [] (const pair<allele,int>& p1, const pair<allele,int>& p2) {
            if (p1.second != p2.second) return p1.second < p2.second;
            if (p1.first.pos != p2.first.pos) return p1.first.pos < p2.first.pos;
            return p1.first.dna < p2.first.dna;
         });
    for (const auto& p : u) {
        ans << YAML::BeginMap;
        ans << YAML::Key << "range" << YAML::Value;
        S(range_yaml(contigs, p.first.pos, ans, true));
        ans << YAML::Key << "dna";
        ans << YAML::Value << p.first.dna;
        ans << YAML::Key << "to";
        ans << YAML::Value << p.second;
        ans << YAML::EndMap;
    }
    ans << YAML::EndSeq;

    ans << YAML::EndMap;
    return Status::OK();
}

Status unified_site::of_yaml(const YAML::Node& yaml, const vector<pair<string,size_t> >& contigs,
                             unified_site& ans) {
    Status s;
    #define V(pred,msg) if (!(pred)) return Status::Invalid("unified_site_of_yaml: " msg)

    V(yaml.IsMap(), "not a map at top level");

    const auto n_range = yaml["range"];
    V(n_range, "missing 'range' field");
    S(range_of_yaml(n_range, contigs, ans.pos));
    #define VR(pred,msg) if (!(pred)) return Status::Invalid("unified_site_of_yaml: " msg, ans.pos.str(contigs))

    const auto n_in_target = yaml["in_target"];
    if (n_in_target) {
        S(range_of_yaml(n_in_target, contigs, ans.in_target));
    }

    ans.alleles.clear();
    const auto n_alleles = yaml["alleles"];
    VR(n_alleles && n_alleles.IsSequence(), "missing 'alleles' field");
    for (YAML::const_iterator al = n_alleles.begin(); al != n_alleles.end(); ++al) {
        V(al->IsScalar(), "invalid allele");
        ans.alleles.push_back(al->Scalar());
    }
    VR(ans.alleles.size() >= 2, "not enough alleles");

    ans.unification.clear();
    const auto n_unification = yaml["unification"];
    VR(n_unification && n_unification.IsSequence(), "missing 'unification' field");
    for (YAML::const_iterator p = n_unification.begin(); p != n_unification.end(); ++p) {
        VR(p->IsMap(), "invalid unification entry");

        range urange(-1,-1,-1);
        const auto n_urange = (*p)["range"];
        VR(n_urange, "missing 'range' field in unification entry");
        S(range_of_yaml(n_urange, contigs, urange, ans.pos.rid));
        VR(ans.pos.overlaps(urange), "unification entry range does not overlap site range");

        const auto n_udna = (*p)["dna"];
        VR(n_udna && n_udna.IsScalar(), "missing/invalid 'dna' field in unification entry");
        const string& dna = n_udna.Scalar();
        VR(dna.size() > 0, "empty 'dna' in unification entry");

        const auto n_uto = (*p)["to"];
        VR(n_uto && n_uto.IsScalar(), "missing/invalid 'to' field in unification entry");
        int to = n_uto.as<int>();
        VR(to >= 0 && to < ans.alleles.size(), "invalid 'to' field in unification entry");

        allele al(urange,dna);
        VR(ans.unification.find(al) == ans.unification.end(), "duplicate unification entries");
        ans.unification[al] = to;
    }
    VR(ans.unification.size() >= 1, "not enough unification entries");

    ans.allele_frequencies.clear();
    const auto n_obs = yaml["allele_frequencies"];
    VR(n_obs && n_obs.IsSequence(), "missing 'allele_frequencies' field");
    for (YAML::const_iterator freq = n_obs.begin(); freq != n_obs.end(); ++freq) {
        VR(freq->IsScalar(), "invalid allele_frequencies");
        float freqf = freq->as<float>();
        if (freqf == freqf) {
            VR(freqf >= 0.0 && freqf <= 1.0, "invalid allele_frequencies");
        }
        ans.allele_frequencies.push_back(freqf);
    }
    VR(ans.allele_frequencies.size() == ans.alleles.size(), "allele_frequencies list has wrong length");

    ans.lost_allele_frequency = 0.0f;
    const auto n_laf = yaml["lost_allele_frequency"];
    if (n_laf) {
        VR(n_laf.IsScalar(), "invalid 'lost_allele_frequency' field");
        ans.lost_allele_frequency = n_laf.as<float>();
        VR(ans.lost_allele_frequency == ans.lost_allele_frequency && ans.lost_allele_frequency >= 0.0f && ans.lost_allele_frequency <= 1.0f, "invalid 'lost_allele_frequency' field");
    }

    const auto n_qual = yaml["quality"];
    VR(n_qual && n_qual.IsScalar(), "missing 'quality' field");
    ans.qual = n_qual.as<int>();
    VR(ans.qual >= 0, "invalid 'qual' field");

    const auto n_monoallelic = yaml["monoallelic"];
    if (n_monoallelic) {
        VR(n_monoallelic.IsScalar(), "invalid 'monoallelic' field");
        ans.monoallelic = n_monoallelic.as<bool>();
    }

    #undef V
    #undef VR
    return Status::OK();
}

// write unified sites list to a file-descriptor, with cap'n proto serialization
static Status capnp_write_unified_sites_fd(const vector<const vector<unified_site>*>& sitess, int fd) {
    ::capnp::MallocMessageBuilder message;
    auto us_b = message.initRoot<capnp::UnifiedSites>();

    size_t total = 0;
    for (const auto& psites : sitess) {
        assert(psites != nullptr);
        total += psites->size();
    }

    auto sites_b = us_b.initSites(total);
    int i = 0;
    for (const auto& psites : sitess) {
        for (const auto& site : *psites) {
            assert(i<total);
            auto site_b = sites_b[i++];
            int j;

            auto pos_b = site_b.initPos();
            pos_b.setRid(site.pos.rid);
            pos_b.setBeg(site.pos.beg);
            pos_b.setEnd(site.pos.end);

            if (site.in_target.rid > -1) {
                auto ct_b = site_b.getContainingTargetOption().initContainingTarget();
                ct_b.setRid(site.in_target.rid);
                ct_b.setBeg(site.in_target.beg);
                ct_b.setEnd(site.in_target.end);
            } else {
                site_b.getContainingTargetOption().setNoContainingTarget(::capnp::VOID);
            }

            auto alleles_b = site_b.initAlleles(site.alleles.size());
            for (j = 0; j < site.alleles.size(); j++) {
                alleles_b.set(j, site.alleles[j]);
            }

            auto unification_b = site_b.initUnification(site.unification.size());
            j = 0;
            for (const auto& p : site.unification) {
                auto oa_b = unification_b[j++];
                auto oa_pos_b = oa_b.initPos();
                oa_pos_b.setRid(p.first.pos.rid);
                oa_pos_b.setBeg(p.first.pos.beg);
                oa_pos_b.setEnd(p.first.pos.end);
                oa_b.setDna(p.first.dna);
                oa_b.setUnifiedAllele(p.second);
            }

            auto allele_frequencies_b = site_b.initAlleleFrequencies(site.allele_frequencies.size());
            for (j = 0; j < site.allele_frequencies.size(); j++) {
                allele_frequencies_b.set(j, site.allele_frequencies[j]);
            }

            site_b.setLostAlleleFrequency(site.lost_allele_frequency);
            site_b.setQual(site.qual);
            site_b.setMonoallelic(site.monoallelic);
        }
    }

    // Capnp throws exceptions on errors, and only in extreme cases.
    writePackedMessageToFd(fd, message);
    return Status::OK();
}


Status capnp_of_unified_sites_fd(const vector<const vector<unified_site>*> &sites, int fd) {
    try {
        return capnp_write_unified_sites_fd(sites, fd);
    } catch (exception &e) {
        return Status::IOError("Capnproto error during de-serialization", e.what());
    }
}

Status capnp_of_unified_sites_fd(const vector<unified_site> &sites, int fd) {
    vector<const vector<unified_site>*> tmp {&sites};
    return capnp_write_unified_sites_fd(tmp, fd);
}

Status capnp_of_unified_sites(const vector<unified_site>& sites, const std::string &filename) {
    int fd = open(filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRWXU);
    if (fd < 0) {
        return Status::IOError("Could not open file for writing", filename);
    }
    Status s = capnp_of_unified_sites_fd(sites, fd);
    if (close(fd) != 0) {
        return Status::IOError("file close failed", filename);
    }
    return s;
}

// read unified_site list from a file-descriptor, as serialized by cap'n proto
static Status _capnp_read_unified_sites_fd(int fd, vector<unified_site>& sites) {
    ::capnp::ReaderOptions ropt;
    ropt.traversalLimitInWords = CAPNP_TRAVERSAL_LIMIT;
    ::capnp::PackedFdMessageReader message(fd, ropt);
    auto sites_r = message.getRoot<capnp::UnifiedSites>();

    sites.clear();
    for(const auto site_r : sites_r.getSites()) {
        range pos(-1,-1,-1);
        const auto pos_r = site_r.getPos();
        pos.rid = pos_r.getRid();
        pos.beg = pos_r.getBeg();
        pos.end = pos_r.getEnd();

        unified_site site(pos);

        const auto cto_r = site_r.getContainingTargetOption();
        if (cto_r.hasContainingTarget()) {
            const auto ct_r = cto_r.getContainingTarget();
            site.in_target.rid = ct_r.getRid();
            site.in_target.beg = ct_r.getBeg();
            site.in_target.end = ct_r.getEnd();
        }

        for (const auto& al : site_r.getAlleles()) {
            site.alleles.push_back(string(al));
        }

        for (const auto oa_r : site_r.getUnification()) {
            range oa_pos(-1,-1,-1);
            auto oa_pos_r = oa_r.getPos();
            oa_pos.rid = oa_pos_r.getRid();
            oa_pos.beg = oa_pos_r.getBeg();
            oa_pos.end = oa_pos_r.getEnd();

            site.unification[allele(oa_pos, oa_r.getDna())] = oa_r.getUnifiedAllele();
        }

        for (const auto x : site_r.getAlleleFrequencies()) {
            site.allele_frequencies.push_back(x);
        }

        site.lost_allele_frequency = site_r.getLostAlleleFrequency();
        site.qual = site_r.getQual();
        site.monoallelic = site_r.getMonoallelic();

        sites.push_back(std::move(site));
    }

    return Status::OK();
}

Status unified_sites_of_capnp(const std::string &filename, vector<unified_site>& sites) {
    int fd = open(filename.c_str(), O_RDONLY);
    if (fd < 0) {
        return Status::IOError("Could not open file for reading", filename);
    }
    Status s;
    try {
        // Capnp throws exceptions on errors, and only in extreme cases.
        s = _capnp_read_unified_sites_fd(fd, sites);
    } catch (exception &e) {
        return Status::IOError("Capnproto error during de-serialization", e.what());
    }
    if (close(fd) != 0) {
        return Status::IOError("file close failed", filename);
    }
    return s;
}

// verify roundtrip of sites
Status capnp_unified_sites_verify(const vector<unified_site>& sites, const string& filename) {
    Status s;
    vector<unified_site> sites2;
    S(capnp_of_unified_sites(sites, filename));
    S(unified_sites_of_capnp(filename, sites2));
    if (sites == sites2) {
        return Status::OK();
    }
    return Status::Invalid("capnp serialization/deserialization of unified sites does not return original value");
}

Status unifier_config::of_yaml(const YAML::Node& yaml, unifier_config& ans) {
    Status s;

    ans = unifier_config();

    #define V(pred,msg) if (!(pred)) return Status::Invalid("unifier_config::of_yaml: " msg)
    V(yaml.IsMap(), "not a map at top level");

    const auto n_drop_filtered = yaml["drop_filtered"];
    if (n_drop_filtered) {
        V(n_drop_filtered.IsScalar(), "invalid drop_filtered");
        ans.drop_filtered = n_drop_filtered.as<bool>();
    }

    const auto n_min_allele_copy_number = yaml["min_allele_copy_number"];
    if (n_min_allele_copy_number) {
        V(n_min_allele_copy_number.IsScalar(), "invalid min_allele_copy_number");
        ans.min_allele_copy_number = n_min_allele_copy_number.as<float>();
        V(ans.min_allele_copy_number >= 0, "invalid min_allele_copy_number");
    }

    const auto n_min_AQ1 = yaml["min_AQ1"];
    if (n_min_AQ1) {
        V(n_min_AQ1.IsScalar(), "invalid min_AQ1");
        ans.min_AQ1 = n_min_AQ1.as<int>();
        V(ans.min_AQ1 >= 0, "invalid min_AQ1");
    }

    const auto n_min_AQ2 = yaml["min_AQ2"];
    if (n_min_AQ2) {
        V(n_min_AQ2.IsScalar(), "invalid min_AQ2");
        ans.min_AQ2 = n_min_AQ2.as<int>();
        V(ans.min_AQ2 >= 0, "invalid min_AQ2");
    }

    const auto n_min_GQ = yaml["min_GQ"];
    if (n_min_GQ) {
        V(n_min_GQ.IsScalar(), "invalid min_GQ");
        ans.min_GQ = n_min_GQ.as<int>();
        V(ans.min_GQ >= 0, "invalid min_GQ");
    }

    const auto n_max_alleles_per_site = yaml["max_alleles_per_site"];
    if (n_max_alleles_per_site) {
        V(n_max_alleles_per_site.IsScalar(), "invalid max_alleles_per_site");
        int max_alleles_per_site = n_max_alleles_per_site.as<int>();
        V(max_alleles_per_site > 0, "invalid max_alleles_per_site");
        ans.max_alleles_per_site = (size_t) max_alleles_per_site;
    }

    const auto n_monoallelic_sites_for_lost_alleles = yaml["monoallelic_sites_for_lost_alleles"];
    if (n_monoallelic_sites_for_lost_alleles) {
        V(n_monoallelic_sites_for_lost_alleles.IsScalar(), "invalid monoallelic_sites_for_lost_alleles");
        ans.monoallelic_sites_for_lost_alleles = n_monoallelic_sites_for_lost_alleles.as<bool>();
    }

    const auto n_preference = yaml["preference"];
    if (n_preference) {
        V(n_preference.IsScalar(), "invalid preference");
        const string& preference = n_preference.Scalar();
        if (preference == "common") {
            ans.preference = UnifierPreference::Common;
        } else if (preference == "small") {
            ans.preference = UnifierPreference::Small;
        } else {
            V(false, "invalid preference");
        }
    }

    #undef V
    return Status::OK();
}

Status unifier_config::yaml(YAML::Emitter& ans) const {
    ans << YAML::BeginMap;
    ans << YAML::Key << "drop_filtered" << YAML::Value << drop_filtered;
    ans << YAML::Key << "min_allele_copy_number" << YAML::Value << min_allele_copy_number;
    ans << YAML::Key << "min_AQ1" << YAML::Value << min_AQ1;
    ans << YAML::Key << "min_AQ2" << YAML::Value << min_AQ2;
    ans << YAML::Key << "min_GQ" << YAML::Value << min_GQ;
    ans << YAML::Key << "max_alleles_per_site" << YAML::Value << max_alleles_per_site;
    ans << YAML::Key << "monoallelic_sites_for_lost_alleles" << YAML::Value << monoallelic_sites_for_lost_alleles;

    ans << YAML::Key << "preference" << YAML::Value;
    if (preference == UnifierPreference::Common) {
        ans << "common";
    } else if (preference == UnifierPreference::Small) {
        ans << "small";
    } else {
        return Status::Invalid("invalid preference");
    }

    ans << YAML::EndMap;

    return Status::OK();
}

Status retained_format_field::yaml(YAML::Emitter& ans) const {
    Status s;
    ans << YAML::Flow;
    ans << YAML::BeginMap;

    ans << YAML::Key << "orig_names";
    ans << YAML::Value << YAML::Flow << YAML::BeginSeq;
    for (auto& name : orig_names) {
        ans << name;
    }
    ans << YAML::EndSeq;

    ans << YAML::Key << "name" << YAML::Value << name;
    ans << YAML::Key << "description" << YAML::Value << description;

    ans << YAML::Key << "type" << YAML::Value;
    if (type == RetainedFieldType::INT) {
        ans << "int";
    } else if (type == RetainedFieldType::FLOAT) {
        ans << "float";
    } else if (type == RetainedFieldType::STRING) {
        ans << "string";
    } else {
        return Status::Invalid("retained_format_field::yaml: invalid type");
    }

    ans << YAML::Key << "number" << YAML::Value;
    if (number == RetainedFieldNumber::BASIC) {
        ans << "basic";
    } else if (number == RetainedFieldNumber::ALT) {
        ans << "alt";
    } else if (number == RetainedFieldNumber::GENOTYPE) {
        ans << "genotype";
    } else if (number == RetainedFieldNumber::ALLELES) {
        ans << "alleles";
    } else {
        return Status::Invalid("retained_format_field::yaml: invalid number");
    }

    ans << YAML::Key << "default_type" << YAML::Value;
    if (default_type == DefaultValueFiller::MISSING) {
        ans << "missing";
    } else if (default_type == DefaultValueFiller::ZERO) {
        ans << "zero";
    } else {
        return Status::Invalid("retained_format_field::yaml: invalid default_type");
    }

    ans << YAML::Key << "count" << YAML::Value << count;

    ans << YAML::Key << "combi_method" << YAML::Value;
    if (combi_method == FieldCombinationMethod::MIN) {
        ans << "min";
    } else if (combi_method == FieldCombinationMethod::MAX) {
        ans << "max";
    } else if (combi_method == FieldCombinationMethod::MISSING) {
        ans << "missing";
    } else if (combi_method == FieldCombinationMethod::SEMICOLON) {
        ans << "semicolon";
    } else {
        return Status::Invalid("retained_format_field::yaml: invalid combi_method");
    }

    ans << YAML::Key << "ignore_non_variants" << YAML::Value;
    if (ignore_non_variants) {
        ans << true;
    } else {
        ans << false;
    }

    ans << YAML::EndMap;

    return Status::OK();
}

Status retained_format_field::of_yaml(const YAML::Node& yaml, unique_ptr<retained_format_field>& ans) {
    Status s;
    #define V(pred,msg) if (!(pred)) return Status::Invalid("retained_format_field::of_yaml: " msg);

    V(yaml.IsMap(), "not a map at top level");

    const auto n_orig_names = yaml["orig_names"];
    V(n_orig_names && n_orig_names.IsSequence() && (n_orig_names.size() > 0), "missing orig_names");
    vector<string> orig_names;
    for (YAML::const_iterator it = n_orig_names.begin(); it != n_orig_names.end(); ++it) {
        V(it->IsScalar(), "invalid orig_names");
        orig_names.push_back(it->Scalar());
    }

    string name;
    const auto n_name = yaml["name"];
    V(n_name && (n_name.Scalar().size() > 0), "missing name");
    name = n_name.Scalar();

    string description;
    const auto n_description = yaml["description"];
    V(n_description && (n_description.Scalar().size() > 0), "missing description");
    description = n_description.Scalar();

    RetainedFieldType type;
    const auto n_type = yaml["type"];
    V(n_type && n_type.IsScalar(), "missing type");
    string s_type = n_type.Scalar();
    if (s_type == "int") {
        type = RetainedFieldType::INT;
    } else if (s_type == "float") {
        type = RetainedFieldType::FLOAT;
    } else if (s_type == "string") {
        type = RetainedFieldType::STRING;
    } else {
        V(false, "invalid type");
    }

    RetainedFieldNumber number;
    const auto n_number = yaml["number"];
    V(n_number && n_number.IsScalar(), "missing number");
    string s_number = n_number.Scalar();
    if (s_number == "basic") {
        number = RetainedFieldNumber::BASIC;
    } else if (s_number == "alt") {
        number = RetainedFieldNumber::ALT;
    } else if (s_number == "genotype") {
        number = RetainedFieldNumber::GENOTYPE;
    } else if (s_number == "alleles") {
        number = RetainedFieldNumber::ALLELES;
    } else {
        V(false, "invalid number");
    }

    int count = 0;
    if (number == RetainedFieldNumber::BASIC) {
        const auto n_count = yaml["count"];
        V(n_count && n_count.IsScalar(), "missing count");
        count = n_count.as<int>();
    }

    DefaultValueFiller default_type = DefaultValueFiller::MISSING;
    const auto n_default_type = yaml["default_type"];
    if (n_default_type) {
        V(n_default_type.IsScalar(), "invalid default_type value");
        string s_default_type = n_default_type.Scalar();

        if (s_default_type == "missing") {
            default_type = DefaultValueFiller::MISSING;
        } else if (s_default_type == "zero") {
            default_type = DefaultValueFiller::ZERO;
        } else {
            V(false, "invalid default_type value");
        }
    }

    FieldCombinationMethod combi_method;
    const auto n_combi_method = yaml["combi_method"];
    V(n_combi_method && n_combi_method.IsScalar(), "missing combi_method");
    string s_combi_method = n_combi_method.Scalar();
    if (s_combi_method == "min") {
        combi_method = FieldCombinationMethod::MIN;
    } else if (s_combi_method == "max") {
        combi_method = FieldCombinationMethod::MAX;
    } else if (s_combi_method == "missing") {
        combi_method = FieldCombinationMethod::MISSING;
    } else if (s_combi_method == "semicolon") {
        combi_method = FieldCombinationMethod::SEMICOLON;
    } else {
        V(false, "invalid combi_method");
    }

    bool ignore_non_variants = false;
    const auto n_ignore_non_variants = yaml["ignore_non_variants"];
    if (n_ignore_non_variants) {
        V(n_ignore_non_variants.IsScalar(), "invalid ignore_non_variants value");
        try {
            ignore_non_variants = n_ignore_non_variants.as<bool>();
        } catch (YAML::Exception e) {
            V(false, "invalid ignore_non_variants value");
        }
    }

    ans.reset(new retained_format_field(orig_names, name, type, combi_method, number, count, default_type, ignore_non_variants));
    ans->description = description;
    #undef V
    return Status::OK();
}

Status genotyper_config::yaml(YAML::Emitter& ans) const {
    Status s;
    ans << YAML::Block;
    ans << YAML::BeginMap;

    ans << YAML::Key << "revise_genotypes" << YAML::Value << revise_genotypes;
    ans << YAML::Key << "min_assumed_allele_frequency" << YAML::Value << min_assumed_allele_frequency;
    ans << YAML::Key << "required_dp" << YAML::Value << required_dp;
    ans << YAML::Key << "allele_dp_format" << YAML::Value << allele_dp_format;
    ans << YAML::Key << "ref_dp_format" << YAML::Value << ref_dp_format;
    ans << YAML::Key << "output_residuals" << YAML::Value << output_residuals;

    ans << YAML::Key << "output_format" << YAML::Value;
    if (output_format == GLnexusOutputFormat::BCF) {
        ans << "BCF";
    } else if (output_format == GLnexusOutputFormat::VCF) {
        ans << "VCF";
    } else {
        return Status::Invalid("genotyper_config::yaml: invalid output_format");
    }

    ans << YAML::Key <<  "liftover_fields";
    ans << YAML::Value << YAML::BeginSeq;
    for (const auto& lo_field : liftover_fields) {
        S(lo_field.yaml(ans));
    }
    ans << YAML::EndSeq;

    ans << YAML::EndMap;

    return Status::OK();
}

Status genotyper_config::of_yaml(const YAML::Node& yaml, genotyper_config& ans) {
    Status s;
    ans = genotyper_config();
    #define V(pred,msg) if (!(pred)) return Status::Invalid("genotyper_config::of_yaml: " msg);
    V(yaml.IsMap(), "not a map at top level");

    const auto n_revise_genotypes = yaml["revise_genotypes"];
    if (n_revise_genotypes) {
        V(n_revise_genotypes.IsScalar(), "invalid revise_genotypes");
        ans.revise_genotypes = n_revise_genotypes.as<bool>();
    }

    const auto n_min_assumed_allele_frequency = yaml["min_assumed_allele_frequency"];
    if (n_min_assumed_allele_frequency) {
        V(n_min_assumed_allele_frequency.IsScalar(), "invalid min_assumed_allele_frequency");
        ans.min_assumed_allele_frequency = n_min_assumed_allele_frequency.as<float>();
        V(ans.min_assumed_allele_frequency >= 0 && ans.min_assumed_allele_frequency <= 1.0, "invalid min_assumed_allele_frequency");
    }

    const auto n_required_dp = yaml["required_dp"];
    if (n_required_dp) {
        V(n_required_dp.IsScalar() && n_required_dp.as<int>() >= 0, "invalid required_dp");
        ans.required_dp = n_required_dp.as<int>();
    }

    const auto n_allele_dp_format = yaml["allele_dp_format"];
    if (n_allele_dp_format) {
        V(n_allele_dp_format.IsScalar(), "invalid allele_dp_format");
        ans.allele_dp_format = n_allele_dp_format.Scalar();
    }

    const auto n_ref_dp_format = yaml["ref_dp_format"];
    if (n_ref_dp_format) {
        V(n_ref_dp_format.IsScalar(), "invalid ref_dp_format");
        ans.ref_dp_format = n_ref_dp_format.Scalar();
    }

    const auto n_output_residuals = yaml["output_residuals"];
    if (n_output_residuals) {
        V(n_output_residuals.IsScalar(), "invalid output_residuals");
        ans.output_residuals = n_output_residuals.as<bool>();
    }

    const auto n_output_format = yaml["output_format"];
    if (n_output_format) {
        V(n_output_format.IsScalar(), "invalid output_format");
        string s_output_format = n_output_format.Scalar();
        if (s_output_format == "BCF") {
            ans.output_format = GLnexusOutputFormat::BCF;
        } else if (s_output_format == "VCF") {
            ans.output_format = GLnexusOutputFormat::VCF;
        } else {
            return Status::Invalid("genotyper_config::of_yaml: invalid output_format. Must be one of {BCF, VCF}.");
        }
    }

    const auto n_liftover_fields = yaml["liftover_fields"];
    if (n_liftover_fields) {
        V(n_liftover_fields.IsSequence(), "invalid liftover_fields");
        for (YAML::const_iterator it = n_liftover_fields.begin(); it != n_liftover_fields.end(); ++it) {
            unique_ptr<retained_format_field> rff;
            S(retained_format_field::of_yaml(*it, rff));
            ans.liftover_fields.push_back(*rff);
        }
    }

    #undef V
    return Status::OK();
}

// regex for a VCF symbolic allele
static std::regex regex_symbolic_allele("<.*>");
bool is_symbolic_allele(const char* allele) {
    return regex_match(allele, regex_symbolic_allele);
}

// gVCF reference confidence records recognized as having either:
// a) zero ALT alleles
// b) one symbolic ALT allele
bool is_gvcf_ref_record(const bcf1_t* record) {
    return record->n_allele == 1 || (record->n_allele == 2 && is_symbolic_allele(record->d.allele[1]));
}

} // namespace GLnexus

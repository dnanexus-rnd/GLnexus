#include <type_traits>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <capnp/message.h>
#include <capnp/serialize-packed.h>

#include <defs.capnp.h>
#include <capnp_serialize.h>

using namespace std;

namespace GLnexus {
namespace capnp {

// We assume the chromosome ploidy is two in our serialization format.
static_assert(zygosity_by_GQ::PLOIDY == 2, "PLOIDY needs to be two");

// write discovered_alleles structure to a file-descriptor, with cap'n proto serialization
Status write_discovered_alleles(int fd,
                                const std::vector<std::pair<std::string,size_t> >& contigs,
                                const discovered_alleles &dsals_i) {
    ::capnp::MallocMessageBuilder message;
    DiscoveredAlleles::Builder dsals = message.initRoot<DiscoveredAlleles>();

    // serialize allele-info pairs
    ::capnp::List<AlleleInfoPair>::Builder aips = dsals.initAips(dsals_i.size());

    int cursor = 0;
    for (auto const &kv : dsals_i) {
        auto &key = kv.first;
        auto &val = kv.second;

        AlleleInfoPair::Builder al = aips[cursor];
        Range::Builder range = al.initRange();
        range.setRid(key.pos.rid);
        range.setBeg(key.pos.beg);
        range.setEnd(key.pos.end);
        al.setDna(key.dna.c_str());

        DiscoveredAlleleInfo::Builder dai = al.initDai();
        dai.setIsRef(val.is_ref);

        TopAQ::Builder topAQ = dai.initTopAQ();
        ::capnp::List<int64_t>::Builder v = topAQ.initV(top_AQ::COUNT);
        for (int k=0; k < top_AQ::COUNT; k++) {
            v.set(k, val.topAQ.V[k]);
        }

        int nelem = val.topAQ.addbuf.size();
        ::capnp::List<int64_t>::Builder addbuf = topAQ.initAddbuf(nelem);
        for (int k=0; k < nelem; k++)
            addbuf.set(k, val.topAQ.addbuf[k]);
        topAQ.setAddbuf(addbuf);

        ::capnp::List<uint64_t>::Builder zGQ0 = dai.initZGQ0(zygosity_by_GQ::GQ_BANDS);
        for (int k=0; k < zygosity_by_GQ::GQ_BANDS; k++)
            zGQ0.set(k, val.zGQ.M[k][0]);

        ::capnp::List<uint64_t>::Builder zGQ1 = dai.initZGQ1(zygosity_by_GQ::GQ_BANDS);
        for (int k=0; k < zygosity_by_GQ::GQ_BANDS; k++)
            zGQ1.set(k, val.zGQ.M[k][1]);

        cursor++;
    }

    writePackedMessageToFd(fd, message);
    return Status::OK();
}

// read discovered_alleles structure from a file-descriptor, as serialized by cap'n proto
Status read_discovered_alleles(int fd,
                               const std::vector<std::pair<std::string,size_t> >& contigs,
                               discovered_alleles &dsals) {
    ::capnp::PackedFdMessageReader message(fd);
    DiscoveredAlleles::Reader dsals_pk = message.getRoot<DiscoveredAlleles>();

    for (AlleleInfoPair::Reader aip : dsals_pk.getAips()) {
        // Unpack allele
        Range::Reader range_pk = aip.getRange();
        range r(range_pk.getRid(),
                range_pk.getBeg(),
                range_pk.getEnd());
        string dna = aip.getDna();
        allele alle(r, dna);

        // Unpack allele-info
        discovered_allele_info dai;
        DiscoveredAlleleInfo::Reader dai_pk = aip.getDai();
        dai.is_ref = dai_pk.getIsRef();

        TopAQ::Reader topAQ_pk = dai_pk.getTopAQ();
        ::capnp::List<int64_t>::Reader v_pk = topAQ_pk.getV();
        if (v_pk.size() != top_AQ::COUNT) {
            return Status::Invalid("Wrong number of elements in topAQ", to_string(v_pk.size()));
        }
        for (int k=0; k < v_pk.size(); k++)  {
            dai.topAQ.V[k] =  v_pk[k];
        }
        ::capnp::List<int64_t>::Reader addbuf_pk = topAQ_pk.getAddbuf();
        for (int k=0; k < addbuf_pk.size(); k++) {
            dai.topAQ.addbuf.push_back(addbuf_pk[k]);
        }

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

Status discover_alleles_verify(const std::vector<std::pair<std::string,size_t> >& contigs,
                               const discovered_alleles &dsals,
                               const string &filename) {
    Status s;

    {
        // Write to file
        int fd = open(filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRWXU);
        S(write_discovered_alleles(fd, contigs, dsals));
        fsync(fd);
        close(fd);
    }

    {
        // read from it
        int fd = open(filename.c_str(), O_RDONLY);
        discovered_alleles dsals2;
        S(read_discovered_alleles(fd, contigs, dsals2));
        close(fd);

/*        cerr << "Reading from file" << endl;
        cerr << "dsals: " << dsals.size() << endl;
        for (auto const &kv : dsals) {
            auto &key = kv.first;
            auto &val = kv.second;
            cerr << key.str() << " " << val.str() << endl;
        }

        cerr << "dsals2: " << dsals2.size() << endl;
        for (auto const &kv : dsals2) {
            auto &key = kv.first;
            auto &val = kv.second;
            cerr << key.str() << " " << val.str() << endl;
            }*/


        // verify we get the same alleles back
        if (dsals == dsals2) {
            return Status::OK();
        }
    }

    return Status::Invalid("capnp serialization/deserialization of discovered alleles does not return original value");
}

}} // namespace GLnexus

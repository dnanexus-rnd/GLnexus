#include <type_traits>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <capnp/message.h>
#include <capnp/serialize-packed.h>

#include <defs.capnp.h>
#include <capnp_serialize.h>

/*
The capnp methods typically do not throw exceptions, or return error
codes that the caller needs to handle.

From the documentation: the implementation prefers to handle
errors using exceptions. Exceptions are only used in circumstances
that should never occur in normal operation. For example, exceptions
are thrown on assertion failures (indicating bugs in the code),
network failures, and invalid input. Exceptions thrown by Cap’n Proto
are never part of the interface and never need to be caught in correct
usage.
*/

using namespace std;

namespace GLnexus {
namespace capnp {

// We assume the chromosome ploidy is two in our serialization format.
static_assert(zygosity_by_GQ::PLOIDY == 2, "PLOIDY needs to be two");

// This limit should be much larger than the message size. It is a security,
// preventing the demarshaller from using too much memory, and getting into
// infinite loops.
const uint64_t TRAVERSAL_LIMIT = 20 * 1024L * 1024L * 1024L;

// write discovered_alleles structure to a file-descriptor, with cap'n proto serialization
static Status _write_discovered_alleles_fd(unsigned int sample_count,
                                           const std::vector<std::pair<std::string,size_t> >& contigs,
                                           const discovered_alleles &dsals,
                                           int fd) {
    ::capnp::MallocMessageBuilder message;
    DiscoveredAlleles::Builder dsals_pk = message.initRoot<DiscoveredAlleles>();
    dsals_pk.setSampleCount(sample_count);

    // serialize contigs
    ::capnp::List<Contig>::Builder contigs_pk = dsals_pk.initContigs(contigs.size());
    int cursor = 0;
    for (auto const &kv : contigs) {
        const string &name = kv.first;
        size_t size = kv.second;

        Contig::Builder ctg_pk = contigs_pk[cursor];
        ctg_pk.setName(name);
        ctg_pk.setSize(size);
        cursor++;
    }

    // serialize allele-info pairs
    ::capnp::List<AlleleInfoPair>::Builder aips_pk = dsals_pk.initAips(dsals.size());
    cursor = 0;
    for (auto const &kv : dsals) {
        auto &key = kv.first;
        auto &val = kv.second;

        AlleleInfoPair::Builder al_pk = aips_pk[cursor];
        Range::Builder range = al_pk.initRange();
        range.setRid(key.pos.rid);
        range.setBeg(key.pos.beg);
        range.setEnd(key.pos.end);
        al_pk.setDna(key.dna.c_str());

        DiscoveredAlleleInfo::Builder dai = al_pk.initDai();
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

    // Capnp throws exceptions on errors, and only in extreme cases.
    writePackedMessageToFd(fd, message);
    return Status::OK();
}

// Variant of the above command, but write to a file
Status write_discovered_alleles_fd(unsigned int sample_count,
                                   const std::vector<std::pair<std::string,size_t> >& contigs,
                                   const discovered_alleles &dsals,
                                   int fd) {
    try {
        return _write_discovered_alleles_fd(sample_count, contigs, dsals, fd);
    } catch (exception &e) {
        return Status::IOError("Capnproto error during de-serialization", e.what());
    }
}

// Variant of the above command, but write to a file
Status write_discovered_alleles(unsigned int sample_count,
                                const std::vector<std::pair<std::string,size_t> >& contigs,
                                const discovered_alleles &dsals,
                                const std::string &filename) {
    std::remove(filename.c_str());
    int fd = open(filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRWXU);
    if (fd < 0)
        return Status::IOError("Could not open file for writing", filename);
    Status s = write_discovered_alleles_fd(sample_count, contigs, dsals, fd);
    int retval = close(fd);
    if (retval < 0)
        return Status::IOError("file close failed", filename);
    return s;
}

// read discovered_alleles structure from a file-descriptor, as serialized by cap'n proto
static Status _read_discovered_alleles_fd(int fd,
                                          unsigned int &sample_count,
                                          std::vector<std::pair<std::string,size_t> >& contigs,
                                          discovered_alleles &dsals) {
    ::capnp::ReaderOptions ropt;
    ropt.traversalLimitInWords =  TRAVERSAL_LIMIT;
    ::capnp::PackedFdMessageReader message(fd, ropt);
    DiscoveredAlleles::Reader dsals_pk = message.getRoot<DiscoveredAlleles>();
    sample_count = dsals_pk.getSampleCount();

    // contigs
    contigs.clear();
    for (Contig::Reader ctg : dsals_pk.getContigs()) {
        auto p = make_pair<string,size_t>(ctg.getName(), ctg.getSize());
        contigs.push_back(p);
    }

    // Discovered Alleles
    dsals.clear();
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

Status read_discovered_alleles(const std::string &filename,
                               unsigned int &sample_count,
                               std::vector<std::pair<std::string,size_t> >& contigs,
                               discovered_alleles &dsals) {
    int fd = open(filename.c_str(), O_RDONLY);
    if (fd < 0)
        return Status::IOError("Could not open file for reading", filename);
    Status s;
    try {
        // Capnp throws exceptions on errors, and only in extreme cases.
        s = _read_discovered_alleles_fd(fd, sample_count, contigs, dsals);
    } catch (exception &e) {
        return Status::IOError("Capnproto error during de-serialization", e.what());
    }

    int retval = close(fd);
    if (retval < 0)
        return Status::IOError("file close failed", filename);
    return s;
}


Status discover_alleles_verify(unsigned int sample_count,
                               const std::vector<std::pair<std::string,size_t> >& contigs,
                               const discovered_alleles &dsals,
                               const string &filename) {
    Status s;

    // Write to file
    S(write_discovered_alleles(sample_count, contigs, dsals, filename));

    // read from it
    unsigned int N;
    vector<pair<string,size_t>> contigs2;
    discovered_alleles dsals2;
    S(read_discovered_alleles(filename, N, contigs2, dsals2));

    // verify we get the same alleles back
    if (dsals == dsals2 &&
        contigs == contigs2 &&
        sample_count == N) {
        return Status::OK();
    }

    return Status::Invalid("capnp serialization/deserialization of discovered alleles does not return original value");
}

}} // namespace GLnexus

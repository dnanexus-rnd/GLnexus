@0x926cdc189d296294;

using Cxx = import "/capnp/c++.capnp";
$Cxx.namespace("GLnexus::capnp");

struct DiscoveredAlleleInfo {
    isRef @0 :Bool = false;

    # top_AQ statistics are used to adjudicate allele existence
    topAQ @1 :List(Int64);

     # zygosity_by_GQ statistics are used to estimate allele copy number
    zGQ0 @2 :List(UInt64);
    zGQ1 @3 :List(UInt64);
}

struct Range {
    rid @0 :Int64;
    beg @1 :Int64;
    end @2 :Int64;
}

# This is really a mapping from allele to information.
struct AlleleInfoPair {
    range @0 :Range;
    dna @1 :Text;
    dai @2 :DiscoveredAlleleInfo;
}

struct Contig {
    name @0 :Text;
    size @1 :UInt64;
}

# aip = allele info pairs
struct DiscoveredAlleles {
    sampleCount @0 : UInt64;
    contigs @1 : List(Contig);
    aips @2 :List(AlleleInfoPair);
}

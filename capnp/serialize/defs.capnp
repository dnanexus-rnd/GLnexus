@0x926cdc189d296294;

using Cxx = import "/capnp/c++.capnp";
$Cxx.namespace("GLnexus::capnp");

struct TopAQ {
    v @0 :List(Int64);
    addbuf @1 :List(Int64);
}

struct DiscoveredAlleleInfo {
    isRef @0 :Bool = false;

    # top_AQ statistics are used to adjudicate allele existence
    topAQ @1 :TopAQ;

     # zygosity_by_GQ statistics are used to estimate allele copy number
    zGQ0 @2 :List(UInt64);
    zGQ1 @3 :List(UInt64);
}

struct Range {
    rid @0 :Int64;
    beg @1 :Int64;
    end @2 :Int64;
}

struct AlleleInfoPair {
    range @0 :Range;
    dna @1 :Text;
    dai @2 :DiscoveredAlleleInfo;
}

# This is really a mapping from allele to information.
struct DiscoveredAlleles {
    # aip = allele info pairs
    aips @0 :List(AlleleInfoPair);
}

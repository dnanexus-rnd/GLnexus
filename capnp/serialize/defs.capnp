@0x926cdc189d296294;

using Cxx = import "/capnp/c++.capnp";
$Cxx.namespace("GLnexus::capnp");

struct TopAQ {
    v @0 :List(UInt64);
    addbuf @1 :List(UInt64);
}

struct ZygosityByGQ {
    m @0 :List(List(UInt64));
}

struct DiscoveredAlleleInfo {
    isRef @0 :Bool = false;

    # top_AQ statistics are used to adjudicate allele existence
    topAQ @1 :TopAQ;

    # zygosity_by_GQ statsitics are used to estimate allele copy number
    zGQ @2 :ZygosityByGQ;
}

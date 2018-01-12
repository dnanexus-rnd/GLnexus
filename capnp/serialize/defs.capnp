@0x926cdc189d296294;

using Cxx = import "/capnp/c++.capnp";
$Cxx.namespace("GLnexus::capnp");

struct Range {
    rid @0 :Int64;
    beg @1 :Int64;
    end @2 :Int64;
}

### discovered alleles

struct DiscoveredAlleleInfo {
    isRef @0 :Bool = false;

    # indicates whether all observations of this allele failed some VCF FILTER
    allFiltered @4 :Bool = false;

    # top_AQ statistics are used to adjudicate allele existence
    topAQ @1 :List(Int64);

     # zygosity_by_GQ statistics are used to estimate allele copy number
    zGQ0 @2 :List(UInt64);
    zGQ1 @3 :List(UInt64);

    inTargetOption : union {
        noInTarget @5 :Void;
        inTarget @6 :Range;
    }
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

### unified sites

struct OriginalAllele {
    pos @0 :Range;
    dna @1 :Text;
    unifiedAllele @2 :Int64;
}

struct UnifiedSite {
    pos @0 :Range;
    inTargetOption :union {
        noInTarget @1 :Void;
        inTarget @2 :Range;
    }
    alleles @3 :List(Text);
    unification @4 :List(OriginalAllele);
    alleleFrequencies @5 :List(Float32);
    lostAlleleFrequency @6 :Float32;
    qual @7 :Int64;
    monoallelic @8 :Bool;
}

struct UnifiedSites {
    sites @0 :List(UnifiedSite);
}

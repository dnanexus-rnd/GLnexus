@0x926cdc189d296294;

using Cxx = import "/capnp/c++.capnp";
$Cxx.namespace("GLnexus::capnp");

struct Range {
    rid @0 :Int64;
    beg @1 :Int64;
    end @2 :Int64;
}

### discovered alleles

struct DiscoveredAllele {
    range @0 :Range;
    dna @1 :Text;
    isRef @2 :Bool = false;

    # optional: discovery target range in which this allele was found
    inTargetOption : union {
        noInTarget @3 :Void;
        inTarget @4 :Range;
    }

    # indicates whether all observations of this allele failed some VCF FILTER
    allFiltered @5 :Bool = false;

    # top_AQ statistics are used to adjudicate allele existence
    topAQ @6 :List(Int64);

     # zygosity_by_GQ statistics are used to estimate allele copy number
    zGQ0 @7 :List(UInt64);
    zGQ1 @8 :List(UInt64);
}

struct Contig {
    name @0 :Text;
    size @1 :UInt64;
}

# aip = allele info pairs
struct DiscoveredAlleles {
    sampleCount @0 : UInt64;
    contigs @1 : List(Contig);
    alleles @2 :List(DiscoveredAllele);
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

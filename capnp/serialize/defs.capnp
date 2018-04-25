@0x926cdc189d296294;

using Cxx = import "/capnp/c++.capnp";
$Cxx.namespace("GLnexus::capnp");

struct Range {
    rid @0 :Int64;
    beg @1 :Int64;
    end @2 :Int64;
}

struct Allele {
    pos @0 :Range;
    dna @1 :Text;
}

### discovered alleles

struct DiscoveredAllele {
    allele @0 :Allele;
    isRef @1 :Bool = false;

    # optional: discovery target range in which this allele was found
    inTargetOption : union {
        noInTarget @2 :Void;
        inTarget @3 :Range;
    }

    # indicates whether all observations of this allele failed some VCF FILTER
    allFiltered @4 :Bool = false;

    # top_AQ statistics are used to adjudicate allele existence
    topAQ @5 :List(Int64);

     # zygosity_by_GQ statistics are used to estimate allele copy number
    zGQ0 @6 :List(UInt64);
    zGQ1 @7 :List(UInt64);
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

struct UnifiedAllele {
    dna @0 :Text;
    normalized @1 :Allele;
    quality @2 :Int64;
    frequency @3 :Float32;
}

struct UnificationEntry {
    representation @0 :Allele;
    unifiedAllele @1 :Int64;
}

struct UnifiedSite {
    pos @0 :Range;
    inTargetOption :union {
        noInTarget @1 :Void;
        inTarget @2 :Range;
    }
    alleles @3 :List(UnifiedAllele);
    unification @4 :List(UnificationEntry);
    lostAlleleFrequency @5 :Float32;
    qual @6 :Int64;
    monoallelic @7 :Bool;
}

struct UnifiedSites {
    sites @0 :List(UnifiedSite);
}



### BCFBucket (for internal database use)
struct BCFBucketSkipEntry {
    recordIndex @0 : UInt32;
    posBeg @1 : Int64;
}
struct BCFBucket {
    records @0 : List(Data);
    skips @1 : List(BCFBucketSkipEntry);
}

#include <iostream>
#include "types.h"
#include "catch.hpp"
using namespace std;
using namespace GLnexus;

TEST_CASE("range::overlaps") {
    range r(0, 1000, 1100);
    REQUIRE(r.overlaps(range(0, 1000, 1100)));
    REQUIRE(r.overlaps(range(0, 1001, 1100)));
    REQUIRE(r.overlaps(range(0, 1000, 1099)));
    REQUIRE(r.overlaps(range(0, 1050, 1090)));
    REQUIRE(r.overlaps(range(0, 999, 1001)));
    REQUIRE_FALSE(r.overlaps(range(0, 666, 999)));
    REQUIRE_FALSE(r.overlaps(range(0, 1101, 1200)));
    REQUIRE_FALSE(r.overlaps(range(1, 1000, 1100)));

    REQUIRE_FALSE(r.overlaps(range(0, 999, 1000)));
    REQUIRE(r.overlaps(range(0, 999, 1001)));

    REQUIRE_FALSE(r.overlaps(range(0, 1100, 1101)));
    REQUIRE(r.overlaps(range(0, 1099, 1100)));
}

TEST_CASE("range sort order") {
    REQUIRE(range(0, 1000, 1100) < range(1, 1000, 1100));
    REQUIRE(range(0, 1000, 1100) < range(1, 999, 1100));
    REQUIRE(range(0, 1000, 1100) < range(0, 1001, 1100));
    REQUIRE(range(0, 1000, 1100) < range(0, 1000, 1101));
    REQUIRE_FALSE(range(0, 1000, 1100) < range(0, 1000, 1100));
    REQUIRE_FALSE(range(0, 1000, 1100) < range(0, 999, 1100));
    REQUIRE_FALSE(range(1, 1000, 1100) < range(0, 999, 1100));
    REQUIRE_FALSE(range(1, 999, 1100) < range(0, 1000, 1100));
    REQUIRE_FALSE(range(0, 1000, 1099) < range(0, 999, 1100));

    vector<range> v { range(1, 1000, 1100), range(0, 1000, 1100), range(0, 1000, 1100), range(1, 999, 1100), range(0, 1001, 1100), range(1, 999, 1100) };
    sort(v.begin(), v.end());

    for (unsigned i = 0; i < v.size()-1; i++) {
        REQUIRE(v[i] <= v[i+1]);
    }
}

TEST_CASE("range_of_bcf") {
    #define UPD(T,name,ini,del) std::unique_ptr<T, void(*)(T*)> up_##name((ini), (del)); auto name = up_##name.get();
    UPD(vcfFile, vcf, bcf_open("test/data/NA12878D_HiSeqX.21.10009462-10009469.gvcf", "r"), [](vcfFile* f) { bcf_close(f); });
    UPD(bcf_hdr_t, hdr, bcf_hdr_read(vcf), &bcf_hdr_destroy);
    shared_ptr<bcf1_t> vt;
    vector<shared_ptr<bcf1_t>> records;

    do {
        if (vt) {
            REQUIRE(bcf_unpack(vt.get(), BCF_UN_ALL) == 0);
            records.push_back(vt);
        }
        vt = shared_ptr<bcf1_t>(bcf_init(), &bcf_destroy);
    } while (bcf_read(vcf, hdr, vt.get()) == 0);

    REQUIRE(records.size() == 5);

    range rng(records[0]);
    REQUIRE(rng.rid == 0);
    REQUIRE(rng.beg == 10009461);
    REQUIRE(rng.end == 10009463);
    REQUIRE(rng.size() == 2);

    rng = range(records[1]);
    REQUIRE(rng.rid == 0);
    REQUIRE(rng.beg == 10009463);
    REQUIRE(rng.end == 10009465);
    REQUIRE(rng.size() == 2);

    rng = range(records[2]);
    REQUIRE(rng.rid == 0);
    REQUIRE(rng.beg == 10009465);
    REQUIRE(rng.end == 10009466);
    REQUIRE(rng.size() == 1);
}

TEST_CASE("discovered_alleles_of_yaml") {
    vector<pair<string,size_t>> contigs;
    contigs.push_back(make_pair("16",12345));
    contigs.push_back(make_pair("17",23456));

    const char* da_yaml = 1 + R"(
- range: {ref: '17', beg: 100, end: 100}
  dna: A
  is_ref: true
  copy_number: 100
- range: {ref: '17', beg: 100, end: 100}
  dna: G
  is_ref: false
  copy_number: 10.5
)";

    #define VERIFY_DA(dal) \
        REQUIRE((dal).size() == 2); \
        REQUIRE((dal).find(allele(range(1, 99, 100),"A")) != (dal).end()); \
        REQUIRE((dal)[allele(range(1, 99, 100),"A")].is_ref == true); \
        REQUIRE((dal)[allele(range(1, 99, 100),"A")].copy_number == 100); \
        REQUIRE((dal).find(allele(range(1, 99, 100),"G")) != (dal).end()); \
        REQUIRE((dal)[allele(range(1, 99, 100),"G")].is_ref == false); \
        REQUIRE((dal)[allele(range(1, 99, 100),"G")].copy_number == 10.5);

    SECTION("basic") {
        YAML::Node n = YAML::Load(da_yaml);

        discovered_alleles dal;
        Status s = discovered_alleles_of_yaml(n, contigs, dal);
        REQUIRE(s.ok());
        VERIFY_DA(dal);
    }

    SECTION("bogus range") {
        YAML::Node n = YAML::Load(1 + R"(
- range: {ref: '17', beg: 100, end: 100}
  dna: A
  is_ref: true
  copy_number: 100
- range: {ref: 'bogus', beg: 100, end: 100}
  dna: G
  is_ref: false
  copy_number: 10.5
)");

        discovered_alleles dal;
        Status s = discovered_alleles_of_yaml(n, contigs, dal);
        REQUIRE(s.bad());

        n = YAML::Load(1 + R"(
- dna: A
  is_ref: true
  copy_number: 100
)");

        s = discovered_alleles_of_yaml(n, contigs, dal);
        REQUIRE(s.bad());

        n = YAML::Load(1 + R"(
- range: 123
  dna: A
  is_ref: true
  copy_number: 100
)");

        s = discovered_alleles_of_yaml(n, contigs, dal);
        REQUIRE(s.bad());
    }

    SECTION("bogus DNA") {
        YAML::Node n = YAML::Load(1 + R"(
- range: {ref: '17', beg: 100, end: 100}
  dna: A
  is_ref: true
  copy_number: 100
- range: {ref: '17', beg: 100, end: 100}
  dna: BOGUS
  is_ref: false
  copy_number: 10.5
)");

        discovered_alleles dal;
        Status s = discovered_alleles_of_yaml(n, contigs, dal);
        REQUIRE(s.bad());
    }

    SECTION("bogus copy_number") {
        YAML::Node n = YAML::Load(1 + R"(
- range: {ref: '17', beg: 100, end: 100}
  dna: A
  is_ref: true
  copy_number: 100
- range: {ref: '17', beg: 100, end: 100}
  dna: G
  is_ref: false
  copy_number: [x]
)");

        discovered_alleles dal;
        Status s = discovered_alleles_of_yaml(n, contigs, dal);
        REQUIRE(s.bad());
    }

    SECTION("duplicates") {
        string dup = string(da_yaml) + string(da_yaml);
        YAML::Node n = YAML::Load(dup.c_str());

        discovered_alleles dal;
        Status s = discovered_alleles_of_yaml(n, contigs, dal);
        REQUIRE(s.bad());
    }
}

TEST_CASE("yaml_of_discovered_alleles") {
    vector<pair<string,size_t>> contigs;
    contigs.push_back(make_pair("16",12345));
    contigs.push_back(make_pair("17",23456));

    SECTION("roundtrip") {
        const char* da_yaml = 1 + R"(
- range: {ref: '17', beg: 100, end: 100}
  dna: A
  is_ref: true
  copy_number: 100
- range: {ref: '17', beg: 100, end: 100}
  dna: G
  is_ref: false
  copy_number: 10.5
)";

        YAML::Node n = YAML::Load(da_yaml);
        discovered_alleles dal;
        REQUIRE(discovered_alleles_of_yaml(n, contigs, dal).ok());

        YAML::Emitter yaml;
        REQUIRE(yaml_of_discovered_alleles(dal, contigs, yaml).ok());
        n = YAML::Load(yaml.c_str());
        discovered_alleles dal2;
        REQUIRE(discovered_alleles_of_yaml(n, contigs, dal2).ok());
        REQUIRE(dal == dal2);
    }
}

TEST_CASE("unified_site::of_yaml") {
    vector<pair<string,size_t>> contigs;
    contigs.push_back(make_pair("16",12345));
    contigs.push_back(make_pair("17",23456));

    const char* snp = 1 + R"(
range: {ref: '17', beg: 100, end: 100}
alleles: [A, G]
copy_number: [100, 51]
unification:
  - range: {ref: '17', beg: 100, end: 100}
    alt: A
    to: 0
  - range: {ref: '17', beg: 100, end: 100}
    alt: G
    to: 1
)";


    #define VERIFY_SNP(us) \
        REQUIRE((us).pos == range(1, 99, 100)); \
        REQUIRE((us).alleles.size() == 2);      \
        REQUIRE((us).alleles[0] == "A");        \
        REQUIRE((us).alleles[1] == "G");        \
        REQUIRE((us).unification.size() == 2);  \
        REQUIRE((us).unification[allele(range(1, 99, 100), "A")] == 0); \
        REQUIRE((us).unification[allele(range(1, 99, 100), "G")] == 1); \
        REQUIRE((us).copy_number.size() == 2); \
        REQUIRE((us).copy_number[0] == 100.0); \
        REQUIRE((us).copy_number[1] == 51.0);


    const char* del = 1 + R"(
range: {ref: '17', beg: 1000, end: 1001}
alleles: [AG, AC, C]
copy_number: [100, 50, 1]
unification:
  - range: {ref: '17', beg: 1000, end: 1001}
    alt: AG
    to: 0
  - range: {ref: '17', beg: 1000, end: 1001}
    alt: AC
    to: 1
  - range: {ref: '17', beg: 1000, end: 1001}
    alt: C
    to: 2
  - range: {ref: '17', beg: 1001, end: 1001}
    alt: C
    to: 1
)";

    #define VERIFY_DEL(us) \
        REQUIRE((us).pos == range(1, 999, 1001)); \
        REQUIRE((us).alleles.size() == 3);      \
        REQUIRE((us).alleles[0] == "AG");       \
        REQUIRE((us).alleles[1] == "AC");       \
        REQUIRE((us).alleles[2] == "C");        \
        REQUIRE((us).unification.size() == 4);  \
        REQUIRE((us).unification[allele(range(1, 999, 1001), "AG")] == 0); \
        REQUIRE((us).unification[allele(range(1, 999, 1001), "AC")] == 1); \
        REQUIRE((us).unification[allele(range(1, 999, 1001), "C")] == 2);  \
        REQUIRE((us).unification[allele(range(1, 1000, 1001), "C")] == 1);
    
    SECTION("snp") {
        YAML::Node n = YAML::Load(snp);

        unified_site us(range(-1,-1,-1));
        Status s = unified_site::of_yaml(n, contigs, us);
        REQUIRE(s.ok());
        VERIFY_SNP(us);
    }

    SECTION("del") {
        YAML::Node n = YAML::Load(del);

        unified_site us(range(-1,-1,-1));
        Status s = unified_site::of_yaml(n, contigs, us);
        REQUIRE(s.ok());
        VERIFY_DEL(us);
    }
    
    SECTION("snp+del") {
        vector<YAML::Node> ns = YAML::LoadAll("---\n" + string(snp) + "\n---\n" + string(del) + "\n...");
        REQUIRE(ns.size() == 2);
        unified_site us(range(-1,-1,-1));
        Status s = unified_site::of_yaml(ns[0], contigs, us);
        REQUIRE(s.ok());
        VERIFY_SNP(us);
        s = unified_site::of_yaml(ns[1], contigs, us);
        REQUIRE(s.ok());
        VERIFY_DEL(us);
    }

    SECTION("optional ref") {
        const char* snp_opt = 1 + R"(
range: {ref: '17', beg: 100, end: 100}
alleles: [A, G]
copy_number: [100, 51]
unification:
  - range: {beg: 100, end: 100}
    alt: A
    to: 0
  - range: {beg: 100, end: 100}
    alt: G
    to: 1
)";
        YAML::Node n = YAML::Load(snp_opt);

        unified_site us(range(-1,-1,-1));
        Status s = unified_site::of_yaml(n, contigs, us);
        REQUIRE(s.ok());
        VERIFY_SNP(us);
    }

    SECTION("bogus range") {
        const char* snp_bogus = 1 + R"(
range: {ref: 'bogus', beg: 100, end: 100}
alleles: [A, G]
copy_number: [100, 51]
unification:
  - range: {beg: 100, end: 100}
    alt: A
    to: 0
  - range: {beg: 100, end: 100}
    alt: G
    to: 1
)";
        YAML::Node n = YAML::Load(snp_bogus);

        unified_site us(range(-1,-1,-1));
        REQUIRE(unified_site::of_yaml(n, contigs, us).bad());

        snp_bogus = 1 + R"(
range: 12345
alleles: [A, G]
copy_number: [100, 51]
unification:
  - range: {beg: 100, end: 100}
    alt: A
    to: 0
  - range: {beg: 100, end: 100}
    alt: G
    to: 1
)";
        n = YAML::Load(snp_bogus);
        REQUIRE(unified_site::of_yaml(n, contigs, us).bad());
    }

    SECTION("bogus alleles") {
        const char* snp_bogus = 1 + R"(
range: {ref: '17', beg: 100, end: 100}
alleles: [A]
copy_number: [100, 51]
unification:
  - range: {beg: 100, end: 100}
    alt: A
    to: 0
)";
        YAML::Node n = YAML::Load(snp_bogus);

        unified_site us(range(-1,-1,-1));
        REQUIRE(unified_site::of_yaml(n, contigs, us).bad());
    }

    SECTION("bogus unification") {
        const char* snp_bogus = 1 + R"(
range: {ref: '17', beg: 100, end: 100}
alleles: [A, G]
copy_number: [100, 51]
unification:
  - range: {ref: 'bogus', beg: 100, end: 100}
    alt: A
    to: 0
  - range: {beg: 100, end: 100}
    alt: G
    to: 1
)";
        YAML::Node n = YAML::Load(snp_bogus);

        unified_site us(range(-1,-1,-1));
        REQUIRE(unified_site::of_yaml(n, contigs, us).bad());

        snp_bogus = 1 + R"(
range: {ref: '17', beg: 100, end: 100}
alleles: [A, G]
copy_number: [100, 51]
unification:
  - range: {beg: 100, end: 100}
    alt: A
    to: 0
)";
        n = YAML::Load(snp_bogus);
        REQUIRE(unified_site::of_yaml(n, contigs, us).bad());
    }
}

TEST_CASE("unified_site::yaml") {
    vector<pair<string,size_t>> contigs;
    contigs.push_back(make_pair("16",12345));
    contigs.push_back(make_pair("17",23456));

    SECTION("roundtrip") {
        const char* del = 1 + R"(
range: {ref: '17', beg: 1000, end: 1001}
containing_target: {ref: '17', beg: 1, end: 10000}
alleles: [AG, AC, C]
copy_number: [100, 50, 1]
unification:
  - range: {ref: '17', beg: 1000, end: 1001}
    alt: AG
    to: 0
  - range: {ref: '17', beg: 1000, end: 1001}
    alt: AC
    to: 1
  - range: {ref: '17', beg: 1000, end: 1001}
    alt: C
    to: 2
  - range: {ref: '17', beg: 1001, end: 1001}
    alt: C
    to: 1
)";
        YAML::Node n = YAML::Load(del);
        unified_site us(range(-1,-1,-1));
        REQUIRE(unified_site::of_yaml(n, contigs, us).ok());

        YAML::Emitter yaml;
        REQUIRE(us.yaml(contigs,yaml).ok());
        n = YAML::Load(yaml.c_str());
        unified_site us2(range(-1,-1,-1));
        REQUIRE(unified_site::of_yaml(n, contigs, us2).ok());
        REQUIRE(us == us2);
    }
}

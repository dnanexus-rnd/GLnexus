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

TEST_CASE("unifier_config") {
    const char* buf1 = 1 + R"(
        {min_allele_copy_number: 5.0,
         max_alleles_per_site:  44,
         preference: common}
)";
    const char* buf2 = 1 + R"(
        {min_allele_copy_number: 5.0,
         max_alleles_per_site:  1,
         preference: small}
)";

    const char* good_examples[] = {buf1, buf2};

    SECTION("good examples") {
        Status s;

        for (const char* buf : good_examples) {
            YAML::Node node = YAML::Load(buf);
            unifier_config uc;
            s = unifier_config::of_yaml(node, uc);
            REQUIRE(s.ok());

            YAML::Emitter ans;
            s = uc.yaml(ans);
            REQUIRE(s.ok());
            node = YAML::Load(ans.c_str());
            unifier_config uc2;
            s = unifier_config::of_yaml(node, uc2);
            REQUIRE(uc == uc2);
        }
    }

    const char* buf3 = 1 + R"(
        {min_allele_copy_number: 5.0,
         max_alleles_per_site:  -1,
         preference: small}
)";
    const char* buf4 = 1 + R"(
        {min_allele_copy_number: 2.0,
         max_alleles_per_site:  1,
         preference: foobar}
)";
    const char* bad_examples[] = {buf3, buf4};

    SECTION("bad examples") {
        Status s;

        for (const char* buf : bad_examples) {
            YAML::Node node = YAML::Load(buf);
            unifier_config uc;
            s = unifier_config::of_yaml(node, uc);
            REQUIRE(s.bad());
        }
    }
}

TEST_CASE("retained_format_field") {
    const char* buf1 = 1 + R"(
        {orig_names: [XXX, YYY],
         name: AAA,
         description: foobar,
         type: int,
         number: genotype,
         default_to_zero: false,
         count: 5,
         combi_method: max}
)";

    const char* good_examples[] = {buf1};

    SECTION("good examples") {
        Status s;

        for (const char* buf : good_examples) {
            YAML::Node node = YAML::Load(buf);
            unique_ptr<retained_format_field> rff;
            s = retained_format_field::of_yaml(node, rff);
            REQUIRE(s.ok());

            YAML::Emitter ans;
            s = rff->yaml(ans);
            REQUIRE(s.ok());
            const char *x = ans.c_str();
            YAML::Node node2 = YAML::Load(x);
            unique_ptr<retained_format_field> rff2;
            s = retained_format_field::of_yaml(node, rff2);

            // We would like to compare the two retained-format-fields,
            // like this:
            //   REQUIRE(*rff == *rff2);
            // but we don't have an equality operator, so, we compare c-strings instead.
            YAML::Emitter ans3;
            s = rff2->yaml(ans3);
            REQUIRE(s.ok());
            const char *y = ans3.c_str();
            REQUIRE(*x == *y);
        }
    }

    const char* buf4 = 1 + R"(
        {orig_names: [XXX, YYY],
         name: AAA,
         description: foobar,
         type: int,
         number: GENOTYPE,
         default_to_zero: false,
         count: 5,
         combi_method: zzzzz}
)";
    const char* bad_examples[] = {buf4};

    SECTION("bad examples") {
        Status s;

        for (const char* buf : bad_examples) {
            YAML::Node node = YAML::Load(buf);
            unique_ptr<retained_format_field> rff;
            s = retained_format_field::of_yaml(node, rff);
            REQUIRE(s.bad());
        }
    }
}


TEST_CASE("genotyper_config") {
    const char* buf1 = 1 + R"(
         required_dp: 3
         allele_dp_format: AD
         ref_symbolic_allele: <NON_REF>
         ref_dp_format: MIN_DP
         output_residuals: true
         output_format: VCF
)";

    const char* buf2 = 1 + R"(
         required_dp: 0
         allele_dp_format: AD
         ref_symbolic_allele: <NON_REF>
         ref_dp_format: MIN_DP
         output_residuals: false
         output_format: BCF
)";

     const char* buf3 = 1 + R"(
required_dp: 0
allele_dp_format: AD
ref_symbolic_allele: <NON_REF>
ref_dp_format: MIN_DP
output_residuals: false
output_format: BCF
liftover_fields:

  - name: AAA
    description: foobar
    type: int
    number: genotype
    default_to_zero: false
    count: 5
    combi_method: max
    orig_names: [XXX, YYY]

  - name: BBB
    description: xxx
    type: int
    number: genotype
    default_to_zero: false
    count: 5
    combi_method: max
    orig_names: [uuu, zzz]

 )";

    const char* good_examples[] = {buf1, buf2, buf3};

    SECTION("good examples") {
        for (const char* buf : good_examples) {
            Status s;
            genotyper_config gc, gc2;

            YAML::Node node = YAML::Load(buf);
            s = genotyper_config::of_yaml(node, gc);
            REQUIRE(s.ok());

            YAML::Emitter ans;
            s = gc.yaml(ans);
            REQUIRE(ans.good());
            REQUIRE(s.ok());
            char const *res1 = ans.c_str();

            node = YAML::Load(res1);
            cout << res1 << endl;
            s = genotyper_config::of_yaml(node, gc2);
            REQUIRE(s.ok());

            // We would like to compare the two structures
            // like this:
            //   REQUIRE(gc == gc2);
            // but we don't have an equality operator, so, we compare c-strings instead.
            YAML::Emitter ans2;
            s = gc2.yaml(ans2);
            REQUIRE(ans2.good());
            REQUIRE(s.ok());
            const char *res2 = ans2.c_str();
            REQUIRE(*res1 == *res2);
        }
    }

    const char* bad_buf1 = 1 + R"(
         {required_dp: -1,
         allele_dp_format: AD,
         ref_symbolic_allele: <NON_REF>,
         ref_dp_format: MIN_DP,
         output_residuals: xxx,
         output_format: BCF}
)";

    const char* bad_buf2 = 1 + R"(
         required_dp: 2
         allele_dp_format: AD
         ref_symbolic_allele: <NON_REF>
         ref_dp_format: MIN_DP
         output_residuals: false
         output_format: xxxxxxxxx
)";

     const char* bad_buf3 = 1 + R"(
required_dp: 0
allele_dp_format: AD
ref_symbolic_allele: <NON_REF>
ref_dp_format: MIN_DP
output_residuals: false
output_format: BCF
liftover_fields:

  - name: CCC
    description: foobar
    type: int
    number: cat
    default_to_zero: false
    count: 5
    combi_method: max
    orig_names: [XXX, YYY]

  - name: DDD
    description: xxx
    type: int
    number: cow
    default_to_zero: false
    count: 5
    combi_method: max
    orig_names: [uuu, zzz]

 )";

    const char* bad_examples[] = {bad_buf1, bad_buf2, bad_buf3};

    SECTION("bad examples") {
        Status s;

        for (const char* buf : bad_examples) {
            YAML::Node node = YAML::Load(buf);
            genotyper_config gc;
            s = genotyper_config::of_yaml(node, gc);
            REQUIRE(s.bad());
        }
    }
}

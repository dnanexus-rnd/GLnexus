// GLnexus crude command-line interface. This is a temporary thing to get us
// bootstrapped with the core algorithms and storage engine before engineering
// an always-on "server"

#include <iostream>
#include <getopt.h>
#include "vcf.h"
#include "hfile.h"
#include "service.h"
#include "alleles.h"
#include "BCFKeyValueData.h"
#include "RocksKeyValue.h"

using namespace std;

GLnexus::Status s;
#define H(desc,expr) \
    s = expr; \
    if (s.bad()) { \
        cerr << "Failed to " << desc << endl << s.str() << endl; \
        return 1; \
    }

void help_init(const char* prog) {
    cerr << "usage: " << prog << " init [options] /desired/db/path exemplar.gvcf[.gz]" << endl
         << "Initializes a new GLnexus database in the specified directory; the parent directory" << endl
         << "must exist but the directory itself must not. The exemplar gVCF file is used to set" << endl
         << "the contigs in the database configuration; it is NOT loaded into the database."
         << endl;
}

int main_init(int argc, char *argv[]) {
    if (argc == 2) {
        help_init(argv[0]);
        return 1;
    }

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    optind = 2; // force optind past command positional argument
    while (-1 != (c = getopt_long(argc, argv, "h",
                                  long_options, nullptr))) {
        switch (c) {
            case 'h':
            case '?':
                help_init(argv[0]);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    if (optind != argc-2) {
        help_init(argv[0]);
        return 1;
    }
    string dbpath(argv[optind]);
    string exemplar_gvcf(argv[optind+1]);

    // load exemplar contigs
    unique_ptr<vcfFile, void(*)(vcfFile*)> vcf(bcf_open(exemplar_gvcf.c_str(), "r"),
                                               [](vcfFile* f) { bcf_close(f); });
    if (!vcf) {
        cerr << "Failed to open exemplar gVCF file at " << exemplar_gvcf << endl;
        return 1;
    }
    unique_ptr<bcf_hdr_t, void(*)(bcf_hdr_t*)> hdr(bcf_hdr_read(vcf.get()), &bcf_hdr_destroy);
    if (!hdr) {
        cerr << "Failed to read gVCF file header from " << exemplar_gvcf << endl;
        return 1;
    }
    int ncontigs = 0;
    const char **contignames = bcf_hdr_seqnames(hdr.get(), &ncontigs);
    vector<pair<string,size_t>> contigs;
    for (int i = 0; i < ncontigs; i++) {
        if (hdr->id[BCF_DT_CTG][i].val == nullptr) {
            cerr << "Invalid gVCF header in " << exemplar_gvcf << endl;
            return 1;
        }
        contigs.push_back(make_pair(string(contignames[i]),
                                    hdr->id[BCF_DT_CTG][i].val->info[0]));
    }

    // create and initialize the database
    unique_ptr<GLnexus::KeyValue::DB> db;
    H("create database", GLnexus::RocksKeyValue::Initialize(dbpath, db));
    H("initialize database", GLnexus::BCFKeyValueData::InitializeDB(db.get(), contigs));

    // report success
    cout << "Initialized GLnexus database in " << dbpath << endl
         << "contigs:";
    for (const auto& contig : contigs) {
        cout << " " << contig.first;
    }
    cout << endl;

    return 0;
}

void help_load(const char* prog) {
    cerr << "usage: " << prog << " load [options] /db/path sample.gvcf[.gz]" << endl
         << "Loads a gVCF file into an existing database. The data set name will be derived from" << endl
         << "the gVCF filename unless overridden."
         << endl;
}

int main_load(int argc, char *argv[]) {
    if (argc == 2) {
        help_load(argv[0]);
        return 1;
    }

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"dataset", required_argument, 0, 'd'},
        {0, 0, 0, 0}
    };

    string dataset;

    int c;
    optind = 2; // force optind past command positional argument
    while (-1 != (c = getopt_long(argc, argv, "hd:",
                                  long_options, nullptr))) {
        switch (c) {
            case 'd':
                dataset = string(optarg);
                if (dataset.size() == 0) {
                    cerr <<  "invalid --dataset" << endl;
                    return 1;
                }
                break;

            case 'h':
            case '?':
                help_load(argv[0]);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    if (optind != argc-2) {
        help_load(argv[0]);
        return 1;
    }
    string dbpath(argv[optind]);
    string gvcf(argv[optind+1]);

    // default dataset name (the gVCF filename)
    if (dataset.size() == 0) {
        size_t p = gvcf.find_last_of('/');
        if (p != string::npos && p < gvcf.size()-1) {
            dataset = gvcf.substr(p+1);
        } else {
            dataset = gvcf;
        }
    }

    // open the database
    unique_ptr<GLnexus::KeyValue::DB> db;
    H("open database", GLnexus::RocksKeyValue::Open(dbpath, db));

    {
        unique_ptr<GLnexus::BCFKeyValueData> data;
        H("open database", GLnexus::BCFKeyValueData::Open(db.get(), data));

        {
            unique_ptr<GLnexus::MetadataCache> metadata;
            H("instantiate metadata cache", GLnexus::MetadataCache::Start(*data, metadata));

            // load the gVCF
            set<string> samples;
            s = data->import_gvcf(*metadata, dataset, gvcf, samples);
            if (s.bad()) {
                cerr << "Failed to load gVCF " << gvcf << endl
                     << s.str() << endl;
                return 1;
            }

            // report success
            cout << "Loaded gVCF " << gvcf << " into database " << dbpath << endl
                 << "dataset: " << dataset << endl
                 << "sample(s):";
            for (const auto& sample : samples) {
                cout << " " << sample;
            }
            cout << endl;
        }
    }

    return 0;
}

void help_dump(const char* prog) {
    cerr << "usage: " << prog << " dump [options] /db/path chrom 1234 2345" << endl
         << "Dump all gVCF records in the database overlapping the given range. The positions"
         << "are one-based, inclusive."
         << endl;
}

int main_dump(int argc, char *argv[]) {
    if (argc == 2) {
        help_dump(argv[0]);
        return 1;
    }

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    optind = 2; // force optind past command positional argument
    while (-1 != (c = getopt_long(argc, argv, "h",
                                  long_options, nullptr))) {
        switch (c) {
            case 'h':
            case '?':
                help_dump(argv[0]);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    if (optind != argc-4) {
        help_dump(argv[0]);
        return 1;
    }
    string dbpath(argv[optind]);
    string rname(argv[optind+1]);
    string beg_txt(argv[optind+2]);
    string end_txt(argv[optind+3]);

    // open the database
    unique_ptr<GLnexus::KeyValue::DB> db;
    H("open database", GLnexus::RocksKeyValue::Open(dbpath, db));

    {
        unique_ptr<GLnexus::BCFKeyValueData> data;
        H("open database", GLnexus::BCFKeyValueData::Open(db.get(), data));

        {
            unique_ptr<GLnexus::MetadataCache> metadata;
            H("instantiate metadata cache", GLnexus::MetadataCache::Start(*data, metadata));

            // resolve the user-supplied contig name to rid
            const auto& contigs = metadata->contigs();
            int rid = 0;
            for(; rid<contigs.size(); rid++) {
                if (contigs[rid].first == rname) {
                    break;
                }
            }
            if (rid == contigs.size()) {
                cerr << "Unknown contig " << rname << endl
                     << "Known contigs:";
                for (const auto& p : contigs) {
                    cerr << " " << p.first;
                }
                cerr << endl;
                return 1;
            }
            // parse positions
            GLnexus::range query(rid, strtol(beg_txt.c_str(), nullptr, 10)-1,
                                      strtol(end_txt.c_str(), nullptr, 10));
            if (query.beg < 0 || query.end < 1 || query.end <= query.beg) {
                cerr << "Invalid query range" << endl;
                return 1;
            }

            // query and output records
            vcfFile *vcfout = vcf_open("-", "w");
            shared_ptr<const set<string>> samples, datasets;
            H("sampleset_datasets", metadata->sampleset_datasets(string("*"), samples, datasets));
            for (const auto& dataset : *datasets) {
                shared_ptr<const bcf_hdr_t> hdr;
                vector<shared_ptr<bcf1_t>> records;
                H("dataset_range", data->dataset_range_and_header(dataset, query, hdr, records));

                int bcf_nsamples = bcf_hdr_nsamples(hdr.get());
                for (const auto& record : records) {
                    cout << dataset << "\t";
                    for (int i = 0; i < bcf_nsamples; i++) {
                        if (i) cout << ",";
                        cout << bcf_hdr_int2id(hdr.get(), BCF_DT_SAMPLE, i);
                    }
                    cout << "\t" << flush;
                    vcf_write(vcfout, hdr.get(), record.get());
                    ignore_retval(hflush(vcfout->fp.hfile));
                }
            }
            // intentionally leaking vcfout because closing it closes stdout
        }
    }

    return 0;
}

void help_genotype(const char* prog) {
    cerr << "usage: " << prog << " genotype [options] /db/path chrom 1234 2345" << endl
         << "Genotype all samples in the database in the given interval. The positions are"
         << "one-based, inclusive."
         << endl;
}

int main_genotype(int argc, char *argv[]) {
    if (argc == 2) {
        help_genotype(argv[0]);
        return 1;
    }

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    optind = 2; // force optind past command positional argument
    while (-1 != (c = getopt_long(argc, argv, "h",
                                  long_options, nullptr))) {
        switch (c) {
            case 'h':
            case '?':
                help_genotype(argv[0]);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    if (optind != argc-4) {
        help_genotype(argv[0]);
        return 1;
    }
    string dbpath(argv[optind]);
    string rname(argv[optind+1]);
    string beg_txt(argv[optind+2]);
    string end_txt(argv[optind+3]);

    // open the database
    unique_ptr<GLnexus::KeyValue::DB> db;
    H("open database", GLnexus::RocksKeyValue::Open(dbpath, db));

    {
        unique_ptr<GLnexus::BCFKeyValueData> data;
        H("open database", GLnexus::BCFKeyValueData::Open(db.get(), data));

        // resolve the user-supplied contig name to rid
        std::vector<std::pair<std::string,size_t> > contigs;
        H("read contig metadata", data->contigs(contigs));
        int rid = 0;
        for(; rid<contigs.size(); rid++) {
            if (contigs[rid].first == rname) {
                break;
            }
        }
        if (rid == contigs.size()) {
            cerr << "Unknown contig " << rname << endl
                 << "Known contigs:";
            for (const auto& p : contigs) {
                cerr << " " << p.first;
            }
            cerr << endl;
            return 1;
        }
        // parse positions
        GLnexus::range query(rid, strtol(beg_txt.c_str(), nullptr, 10)-1,
                                  strtol(end_txt.c_str(), nullptr, 10));
        if (query.beg < 0 || query.end < 1 || query.end <= query.beg) {
            cerr << "Invalid query range" << endl;
            return 1;
        }

        {
            // start service, discover alleles, unify sites, genotype sites
            unique_ptr<GLnexus::Service> svc;
            H("start GLnexus service", GLnexus::Service::Start(*data, *data, svc));

            const string sampleset("*");
            GLnexus::discovered_alleles alleles;
            H("discover alleles", svc->discover_alleles(sampleset, query, alleles));

            vector<GLnexus::unified_site> sites;
            H("unify sites", GLnexus::unified_sites(alleles, sites));

            GLnexus::consolidated_loss losses;
            H("genotype sites",
              svc->genotype_sites(GLnexus::genotyper_config(), sampleset, sites, string("-"), losses));

            cerr << "\nReporting loss for " << sites.size() << " site(s) genotyped for "<< losses.size() << " sample(s)." << endl;

            cerr << "============" << endl;
            for (auto& loss : losses) {
                cerr << "Sample " << loss.first << ": ";
                cerr << loss.second.str() << endl;
            }
        }
    }

    return 0;
}

void help(const char* prog) {
    cerr << "usage: " << prog << " <command> [options]" << endl
             << endl
             << "commands:" << endl 
             << "  init     initialize new database" << endl
             << "  load     load a gVCF file into an existing database" << endl
             << "  genotype genotype samples in the database" << endl
             << endl;
}

int main(int argc, char *argv[]) {
    
    if (argc == 1) {
        help(argv[0]);
        return 1;
    }

    string command = argv[1];
    if (command == "init") {
        return main_init(argc, argv);
    } else if (command == "load") {
        return main_load(argc, argv);
    } else if (command == "dump") {
        return main_dump(argc, argv);
    } else if (command == "genotype") {
        return main_genotype(argc, argv);
    } else {
        cerr << "unknown command " << command << endl;
        help(argv[0]);
        return 1;
    }

    return 0;
}

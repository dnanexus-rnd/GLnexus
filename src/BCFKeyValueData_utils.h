// Helper code for BCFKeyValueData.cc

namespace GLnexus {

// Memory efficient representation of a bucket range. This could
// be turned into a standard C++ iterator, although, that might be
// a bit of an overkill.
class BucketExtent {
private:
    int rid_ = 0;
    int step_ = 0;
    int bgn_ = 0;
    int end_ = 0;
    int current_ = 0;

    // disable copy and assignment constructors
    BucketExtent(const BucketExtent&);
    BucketExtent& operator=(const BucketExtent&);

public:
    BucketExtent(const range &query, int interval_len) {
        rid_ = query.rid;
        step_ = interval_len;
        bgn_ = (query.beg / interval_len) * interval_len;
        end_ = std::max(bgn_, ((query.end-1) / interval_len) * interval_len);
        current_ = 0;
    }

    range begin() {
        range r = range(rid_, bgn_, bgn_ + step_);
        current_ = bgn_;
        return r;
    }

    range next() {
        current_ += step_;
        range r = range(rid_, current_, current_ + step_);
        return r;
    }

    range end() {
        range r = range(rid_, end_, end_ + step_);
        return r;
    }
};

// Map a range into a set of buckets. The records in the range
// can be found by scanning all the buckets. Note that a BCF record
// is placed in a bucket based on its start position.
// It could start in one bucket, and extend into an adjacent bucket(s).
//
// This class separates out the logic for answering the following questions:
//   1) Which buckets should I scan for this query range?
//   2) Which bucket does a bcf1_t with this range go into?
class BCFBucketRange {
private:
    // disable copy and assignment constructors
    BCFBucketRange(const BCFBucketRange&);
    BCFBucketRange& operator=(const BCFBucketRange&);

public:
    static const size_t PREFIX_LENGTH = 8;
    int interval_len;

    // constructor
    BCFBucketRange(int interval_len) : interval_len(interval_len) {};

    // Given the range of a bucket, produce the key prefix for the bucket.
    // Important: the range must be exactly that of the bucket.
    // BCFBucketRange::bucket below translates an arbitrary range into a
    // bucket's range.
    std::string bucket_prefix(const range& rng) {
        uint64_t rid_be = htobe64(rng.rid);
        uint64_t beg_be = htobe64(rng.beg);
        assert(be64toh(rid_be) < MAX_NUM_CONTIGS_PER_GVCF);
        assert(be64toh(beg_be) < MAX_CONTIG_LEN);
        char buf[8]; static_assert(PREFIX_LENGTH == 8, "assumption");
        memcpy(buf, ((char*)&rid_be)+5, 3);
        memcpy(buf+3, ((char*)&beg_be)+3, 5);
        return string(buf, 8);
    }

    // Produce the complete key for a bucket (given the prefix) in a dataset
    std::string bucket_key(const std::string& prefix, const std::string& dataset) {
        assert(prefix.size() == PREFIX_LENGTH);
        return prefix+dataset;
    }

    // Same as bucket_key(bucket_prefix(rng), dataset)
    // Important: the range must be exactly that of the bucket.
    // BCFBucketRange::bucket below translates an arbitrary range into a
    // bucket's range.
    std::string bucket_key(const range& rng, const std::string& dataset) {
        return bucket_key(bucket_prefix(rng), dataset);
    }

    // Decompose the key into bucket prefix and dataset
    Status parse_key(const string& key, string& bucket, string& dataset) {
        if (key.size() < PREFIX_LENGTH) {
            return Status::Invalid("BCFBucketRange::parse_key: key too small", key);
        }
        bucket = key.substr(0, PREFIX_LENGTH);
        dataset = key.substr(PREFIX_LENGTH);
        return Status::OK();
    }

    // Given arbitrary [query] range, return a structure describing one or
    // more buckets to search through in order to find all records overlapping
    // query. This may be multiple buckets, even for small query ranges, to
    // account for the possibility of records spanning multiple buckets.
    std::shared_ptr<BucketExtent> scan(const range& query) {
        return make_shared<BucketExtent>(query, interval_len);
    }

    // Which bucket does this BCF record start in?
    range bucket(bcf1_t *rec) {
        int bgn = (rec->pos / interval_len) * interval_len;
        return range(rec->rid, bgn, bgn + interval_len);
    }
    // The bucket after [rng], assuming [rng] is a bucket.
    range inc_bucket(range &rng) {
        assert((rng.end - rng.beg) == interval_len);
        return range(rng.rid,
                     rng.beg + interval_len,
                     rng.end + interval_len);
    }

    // Create a ficticious bucket marking the end of a chromosome.
    range bucket_at_end_of_chrom(int rid,
                                 const std::vector<std::pair<std::string,size_t> >&contigs) {
        //const string &contig_name = contigs[rid].first;
        size_t contig_len = contigs[rid].second;
        int bgn = ((contig_len / interval_len) + 2) * interval_len;
        return range(rid, bgn, bgn + interval_len);
    }
};

// A "BCF Bucket" is the value serialized into the database containing some
// number of BCF records. The records are ordered by position and must all
// lie on the same contig. They may overlap.
//
// Range queries within the bucket may be performed by a linear scan of the
// records. Such a scan can be truncated upon seeing a record whose beg is
// greater than the end of the query range. For further optimization, we
// save a "skip index" along with the list of records. An entry in the skip
// index has an index into the list of records, and the position.beg of the
// corresponding record. A record may be represented in the skip index
// only if no preceding records in the bucket overlap it. Given this, the
// scan can begin at the index indicated by the last skip index entry whose
// beg is <= the query beg, thus 'skipping' the preceding records.

class BCFBucketWriter {
    vector<vector<uint8_t>> records_;
    int rid_, last_beg_, end_;
    vector<pair<int,int>> skips_;

public:
    BCFBucketWriter()
        : rid_(-1), last_beg_(-1), end_(-1)
        {
    }

    void clear() {
        records_.clear();
        skips_.clear();
        rid_ = last_beg_ = end_ = -1;
    }

    Status add(bcf1_t* rec) {
        range rng(rec);
        if (rid_ == -1) {
            rid_ = rng.rid;
        } else if (rid_ != rng.rid) {
            return Status::Invalid("BCFBucketWriter: contig mismatch (BUG)");
        }
        if (rng.beg < last_beg_) {
            return Status::Invalid("BCFBucketWriter: records not sorted (BUG)");
        }
        last_beg_ = rng.beg;
        if (end_ <= rng.beg && records_.size()/10 > skips_.size()) {
            // no preceding records in the bucket overlap this one, so record
            // this skip index entry (for every 10th record, 10 arbitrary)
            skips_.push_back(make_pair((int)records_.size(), rng.beg));
        }
        end_ = max(end_, rng.end);

        size_t reclen = bcf_raw_calc_packed_len(rec);
        assert(reclen > 0);
        vector<uint8_t> buf(reclen);
        bcf_raw_write_to_mem(rec, reclen, buf.data());
        records_.push_back(move(buf));
        return Status::OK();
    }

    size_t get_num_entries() const {
        return records_.size();
    }

    Status contents(string& ans) const {
        try {
            ::capnp::MallocMessageBuilder b;
            auto msg_b = b.initRoot<capnp::BCFBucket>();
            auto records_b = msg_b.initRecords(records_.size());
            for (int i = 0; i < records_.size(); i++) {
                records_b.set(i, kj::arrayPtr((kj::byte*) records_[i].data(), records_[i].size()));
                assert(records_b[i].begin() != nullptr); assert(records_b[i].size() == records_[i].size());
            }

            if (skips_.size()) {
                auto skips_b = msg_b.initSkips(skips_.size());
                for (int i = 0; i < skips_.size(); i++) {
                    skips_b[i].setRecordIndex(skips_[i].first);
                    skips_b[i].setPosBeg(skips_[i].second);
                }
            }

            auto msg_words = ::capnp::messageToFlatArray(b);
            auto msg_bytes = msg_words.asBytes();
            ans.assign((char*)msg_bytes.begin(), msg_bytes.size());

            #ifndef NDEBUG
            {
                // verify correct deserialization, with different memory alignments
                // more info about why memory alignment is of interest below, in ScanBCFBucket.
                int alignments_to_test = sizeof(::capnp::word);
                #ifndef __x86_64__
                alignments_to_test = 1;
                #endif
                for (int ofs = 0; ofs < alignments_to_test; ofs++) {
                    auto buf = (uint8_t*) calloc(ans.size() + ofs, 1);
                    memcpy(buf+ofs, ans.data(), ans.size());
                    ::capnp::FlatArrayMessageReader message(kj::ArrayPtr<const ::capnp::word>((::capnp::word*)(buf+ofs), ans.size() / sizeof(::capnp::word)));
                    capnp::BCFBucket::Reader bucket_reader = message.getRoot<capnp::BCFBucket>();
                    auto records = bucket_reader.getRecords();
                    assert(records.size() == records_.size());
                    for (int i = 0; i < records.size(); i++) {
                        assert(records[i].size() == records_[i].size());
                        assert(memcmp(records[i].begin(), records_[i].data(), records[i].size()) == 0);
                    }
                    auto skips = bucket_reader.getSkips();
                    assert(skips.size() == skips_.size());
                    for (int i = 0; i < skips.size(); i++) {
                        assert(skips[i].getRecordIndex() == skips_[i].first);
                        assert(skips[i].getPosBeg() == skips_[i].second);
                    }
                    free(buf);
                }
            }
            #endif

            return Status::OK();
        } catch (exception &e) {
            return Status::IOError("exception serializing BCF bucket", e.what());
        }
    }
};

// Search the bucket's 'skip index' to find the index of a record in the bucket
// from which it is safe to begin a scan for records overlapping query.
static int SearchBCFBucketSkipIndex(const capnp::BCFBucket::Reader& bucket, const range& query) {
    auto skips = bucket.getSkips();
    if (skips.size() == 0 || skips[0].getPosBeg() > query.beg) {
        return 0;
    }
    int i = 0, hi = skips.size(); // invariant: skips[i].getPosBeg() <= query.beg
    while (i < hi-1) {
        int mid = (i + hi) / 2;
        if (skips[mid].getPosBeg() <= query.beg) {
            assert(mid > i);
            i = mid;
        } else {
            assert(mid < hi);
            hi = mid;
        }
    }
    #ifndef NDEBUG
    int naive;
    for (naive = skips.size()-1; naive >= 0; naive--) {
        if (skips[naive].getPosBeg() <= query.beg) {
            break;
        }
    }
    assert(i == naive);
    #endif
    return skips[i].getRecordIndex();
}

// test whether a gVCF file is compatible for deposition into the database
static bool gvcf_compatible(const MetadataCache& metadata, const bcf_hdr_t *hdr) {
    Status s;
    const auto& contigs = metadata.contigs();

    // verify contigs match exactly. even the order matters

    int ncontigs = 0;
    const char **contignames = bcf_hdr_seqnames(hdr, &ncontigs);
    bool ans = true;

    if (((uint)ncontigs) != contigs.size()) {
        ans = false;
    } else {
        for (int i = 0; i < ncontigs; i++) {
            if (string(contignames[i]) != contigs[i].first ||
                hdr->id[BCF_DT_CTG][i].val->info[0] != contigs[i].second) {
                ans = false;
                break;
            }
        }
    }

    free(contignames);
    return ans;
}

// Sanity-check an individual bcf1_t record before ingestion.
static Status validate_bcf(BCFBucketRange& rangeHelper,
                           const std::vector<std::pair<std::string,size_t> >&contigs,
                           const std::string& filename,
                           const bcf_hdr_t *hdr,
                           bcf1_t *bcf,
                           int prev_rid, int prev_pos) {
    if (bcf_unpack(bcf, BCF_UN_ALL) != 0 || bcf->errcode != 0) {
        return Status::Invalid("invalid VCF record (corrupt format, or undeclared info/format fields)",
                               filename + " " + range(bcf).str(contigs));
    }

    // Check that bcf->rlen is calculated correctly based on POS,END if
    // available or POS,strlen(REF) otherwise
    bcf_info_t *info = bcf_get_info(hdr, bcf, "END");
    if (info) {
        if (info->type != BCF_BT_INT8 && info->type != BCF_BT_INT16 && info->type != BCF_BT_INT32) {
            return Status::Invalid("gVCF record's END field has unexpected type", filename + " " + range(bcf).str(contigs));
        }
        if (info->len != 1) {
            return Status::Invalid("gVCF record has multiple END fields", filename + " " + range(bcf).str(contigs));
        }
        if (info->v1.i < bcf->pos) {
            return Status::Invalid("gVCF record has END < POS", filename + " " + range(bcf).str(contigs));
        }
        if (info->v1.i - bcf->pos != bcf->rlen) {
            return Status::Invalid("gVCF record END-POS doesn't match rlen", filename + " " + range(bcf).str(contigs));
        }
    } else {
        if (bcf->d.allele == nullptr) {
            return Status::Invalid("gVCF allele is null", filename + " " + range(bcf).str(contigs));
        }
        if (bcf->rlen != (int) strlen(bcf->d.allele[0])) {
            return Status::Invalid("gVCF rlen doesn't match strlen(REF) (and no END field)", filename + " " + range(bcf).str(contigs));
        }
    }

    // verify record ordering is non-decreasing within a contig
    if (prev_rid == bcf->rid && prev_pos > bcf->pos) {
        return Status::Invalid("gVCF records are out-of-order ",
                               filename + " " + std::to_string(prev_pos+1) + " >= " + range(bcf).str(contigs));
    }

    // verify record does not go over the length of the contig
    const string &contig_name = contigs[bcf->rid].first;
    size_t contig_len = contigs[bcf->rid].second;
    if (bcf->pos + bcf->rlen > contig_len) {
        return Status::Invalid("gVCF record is longer than contig ",
                               filename + " " + range(bcf).str(contigs) + " " + to_string(contig_len) + " " + contig_name);
    }

    // check that alleles are all distinct, and that all alleles are valid strings
    // of IUPAC nucleotides, except the last ALT allele which may be symbolic.
    set<string> alleles;
    for (int i=0; i < bcf->n_allele; i++) {
        const string allele_i(bcf->d.allele[i]);
        if (!(is_iupac_nucleotides(allele_i) ||
              (i == bcf->n_allele-1 && is_symbolic_allele(allele_i.c_str())))) {
            return Status::Invalid("allele is not a DNA sequence ",
                                filename + " " + allele_i +  " " + range(bcf).str(contigs));
        }
        alleles.insert(allele_i);
    }
    if (bcf->n_allele<1 || alleles.size() != bcf->n_allele) {
        return Status::Invalid("alleles are not distinct ", filename + " " + range(bcf).str(contigs));
    }

    // validate genotypes (all entries either missing or in [0, n_allele))
    if (bcf->n_sample != bcf_hdr_nsamples(hdr)) {
        return Status::Invalid("gVCF record doesn't have expected # of samples", filename + " " + range(bcf).str(contigs));
    }
    htsvecbox<int> gt;
    int nGT = bcf_get_genotypes(hdr, bcf, &gt.v, &gt.capacity);
    if (nGT != 2*bcf->n_sample) {
        return Status::Invalid("gVCF record doesn't have expected # of GT entries", filename + " " + range(bcf).str(contigs));
    }
    for (int i = 0; i < nGT; i++) {
        if (!bcf_gt_is_missing(gt[i]) && (bcf_gt_allele(gt[i]) < 0 || bcf_gt_allele(gt[i]) >= bcf->n_allele)) {
            return Status::Invalid("invalid GT entry in gVCF record", filename + " " + range(bcf).str(contigs));
        }
    }

    // validate genotype likelihoods (n_samples*nGT entries; PL: all entries nonnegative; GL: all entries nonpositive)
    htsvecbox<int32_t> pl;
    int nPL = bcf_get_format_int32(hdr, bcf, "PL", &pl.v, &pl.capacity);
    if (nPL >= 0) {
        if (nPL != bcf->n_sample * diploid::genotypes(bcf->n_allele) && bcf->n_allele > 1) {
            // the exception when bcf->n_allele == 1 accommodates xAtlas
            return Status::Invalid("gVCF record doesn't have expected # of PL entries", filename + " " + range(bcf).str(contigs));
        }
        for (int i = 0; i < nPL; i++) {
            if (pl[i] < 0) {
                return Status::Invalid("negative PL entry in gVCF record", filename + " " + range(bcf).str(contigs));
            }
        }
    }
/*
    htsvecbox<float> gl;
    int nGL = bcf_get_format_float(hdr, bcf, "GL", &gl.v, &gl.capacity);
    if (nGL >= 0) {
        if (nGL != bcf->n_sample * diploid::genotypes(bcf->n_allele)) {
            return Status::Invalid("gVCF record doesn't have expected # of GL entries", filename + " " + range(bcf).str(contigs));
        }
        for (int i = 0; i < nGL; i++) {
            if (gl[i] > 0.0) {
                return Status::Invalid("positive GL entry in gVCF record", filename + " " + range(bcf).str(contigs));
            }
        }
    }
*/
    return Status::OK();
}


// Verify that a VCF file is well formed.
// AND, fill in the [samples_out]
static Status vcf_validate_basic_facts(MetadataCache& metadata,
                                       const string& dataset,
                                       const string& filename,
                                       bcf_hdr_t *hdr,
                                       vcfFile *vcf,
                                       set<string>& samples_out)
{
    if (!hdr) return Status::IOError("reading gVCF header", filename);
    if (!gvcf_compatible(metadata, hdr)) {
        return Status::Invalid("Incompatible gVCF. The reference contigs must match the database configuration exactly.", filename);
    }

    vector<string> samples;
    unsigned n = bcf_hdr_nsamples(hdr);
    if (n == 0) {
        return Status::Invalid("gVCF contains no samples", dataset + " (" + filename + ")");
    }
    for (unsigned i = 0; i < n; i++) {
        string sample(bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i));
        if (!regex_match(sample, regex_id)) {
            return Status::Invalid("gVCF contains invalid sample name", dataset + " (" + filename + ")" + " " + sample);
        }
        samples.push_back(move(sample));
    }
    samples_out.clear();
    samples_out.insert(samples.begin(), samples.end());
    if (samples.size() != samples_out.size()) {
        return Status::Invalid("gVCF sample names are not unique", dataset + " (" + filename + ")");
    }

    return Status::OK();
}

// helper class for bulk_insert_gvcf_key_values: accumulate sizable batches of
// key/value pairs before insertion into the KeyValue database.
// This is to reduce database write lock contention during intense multi-
// threaded bulk loads, as each thread makes fewer larger inserts instead
// of many smaller inserts.
class BulkInsertBuffer {
    const size_t LIMIT = 16777216;
    KeyValue::DB& db_;
    std::unique_ptr<KeyValue::WriteBatch> buf_;
    size_t bufsz_ = 0;

public:
    BulkInsertBuffer(KeyValue::DB& db) : db_(db) {}
    ~BulkInsertBuffer() {
        assert(!buf_);
    }

    Status put(KeyValue::CollectionHandle coll, const std::string& key, const std::string& value) {
        Status s;
        size_t delta = key.size() + value.size() + 32;
        if (bufsz_ + delta >= LIMIT) {
            S(flush());
        }
        if (!buf_) {
            S(db_.begin_writes(buf_));
        }
        S(buf_->put(coll, key, value));
        bufsz_ += delta;
        return Status::OK();
    }

    // make sure to call when finished
    Status flush() {
        if (buf_ && bufsz_) {
            Status s;
            S(buf_->commit());
        }
        buf_.reset();
        bufsz_ = 0;
        return Status::OK();
    }
};

} // namespace GLnexus
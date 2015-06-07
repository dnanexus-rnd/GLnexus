#ifndef GLNEXUS_KEYVALUE_H
#define GLNEXUS_KEYVALUE_H
#include "data.h"

namespace GLnexus {
namespace KeyValue {

/// Abstract interface to a database underlying BCFKeyValueData. The database
/// has one or more collections of key-value records. Each collection is
/// ordered by key.

/// In-order iterator over records in a collection
class Iterator {
public:
    /// Update key and value with the next record and return OK. Return
    /// NotFound if there are no remaining records, or any error code
    virtual Status next(std::string& key, std::string& value) = 0;
};

/// A DB snapshot providing consistent multiple reads if possible
template<typename CollectionHandle, class IteratorImpl>
class Reader {
public:
    using iterator_type = IteratorImpl;
    static_assert(std::is_base_of<Iterator, IteratorImpl>::value, "IteratorImpl must implement Iterator interface");

    /// Get the value corresponding to the key and return OK. Return NotFound
    /// if no corresponding record exists in the collection, or any error code
    virtual Status get(const CollectionHandle& coll, const std::string& key, std::string& value) const = 0;

    /// Create an iterator over the whole collection.
    virtual Status iterator(const CollectionHandle& coll, std::unique_ptr<IteratorImpl>& it) const = 0;

    /// Create an iterator beginning at the key if a corresponding record
    /// exists, or the first subsequent record otherwise
    virtual Status iterator(const CollectionHandle& coll, const std::string& key, std::unique_ptr<IteratorImpl>& it) const = 0;
};

/// A batch of writes to apply atomically if possible
template<typename CollectionHandle>
class WriteBatch {
public:
    virtual Status put(const CollectionHandle& coll, const std::string& key, const std::string& value) = 0;
    //virtual Status delete(Collection* coll, const std::string& key) = 0;
};

/// Main database interface for retrieving collection handles, generating
/// snapshopts to read from, and creating and applying write batches. The DB
/// object itself implements the Reader interface (with no consistency
/// guarantees between multiple calls) and the WriteBatch interface (which
/// applies one write immediately, no atomicity guarantees between multiple
/// calls). Caller must ensure that the parent DB object still exists when any
/// Reader or WriteBatch object is used.
template<typename CollectionHandle, class ReaderImpl, class IteratorImpl, class WriteBatchImpl>
class DB : public ReaderImpl, public WriteBatchImpl {
public:
    using collection_handle_type = CollectionHandle;
    using iterator_type = IteratorImpl;
    using write_batch_type = WriteBatchImpl;
    static_assert(std::is_base_of<Reader<CollectionHandle,IteratorImpl>, ReaderImpl>::value, "ReaderImpl must implement Reader interface");
    static_assert(std::is_base_of<WriteBatch<CollectionHandle>, WriteBatchImpl>::value, "WriterImpl must implement Writer interface");
    static_assert(std::is_base_of<Iterator, IteratorImpl>::value, "IteratorImpl must implement Iterator interface");
    
    /// Get the handle to a collection, or return NotFound.
    virtual Status collection(const std::string& name, CollectionHandle& coll) const = 0;

    /// Create a new collection, or return Exists.
    virtual Status create_collection(const std::string& name) = 0;

    /// Get an up-to-date snapshot.
    virtual Status current(std::unique_ptr<ReaderImpl>& snapshot) const = 0;

    /// Begin preparing a batch of writes.
    virtual Status begin_writes(std::unique_ptr<WriteBatchImpl>& writes) = 0;

    /// Apply a batch of writes.
    virtual Status commit_writes(WriteBatchImpl* writes) = 0;

    // Base implementations of Reader and WriteBatch interfaces. They simply
    // create a snapshot just to read one record (or begin one iterator), or
    // apply a "batch" of one write. Derived classes may want to provide more
    // efficient overrides.
    Status get(const CollectionHandle& coll, const std::string& key, std::string& value) const override;
    Status iterator(const CollectionHandle& coll, std::unique_ptr<IteratorImpl>& it) const override;
    Status iterator(const CollectionHandle& coll, const std::string& key, std::unique_ptr<IteratorImpl>& it) const override;
    Status put(const CollectionHandle& coll, const std::string& key, const std::string& value) override;
};

// Trivial in-memory KeyValue::DB implementation used in unit tests. I'd
// rather define this in the unit test code, but, C++ templates aren't
// accommodating...
namespace Mem {
    class Iterator : public KeyValue::Iterator {
        std::map<std::string,std::string> data_;
        std::map<std::string,std::string>::const_iterator it_;
        friend class Reader;

    public:
        Status next(std::string& key, std::string& value) override {
            if (it_ == data_.end()) return Status::NotFound();
            key = it_->first;
            value = it_->second;
            it_++;
            return Status::OK();
        }
    };

    class Reader : public KeyValue::Reader<uint64_t,Iterator> {
        std::vector<std::map<std::string,std::string>> data_;
        friend class DB;

    public:
        Status get(const uint64_t& coll, const std::string& key, std::string& value) const override {
            assert(coll < data_.size());
            const auto& m = data_[coll];
            auto p = m.find(key);
            if (p == m.end()) return Status::NotFound("key", key);
            value = p->second;
            return Status::OK();
        }

        Status iterator(const uint64_t& coll, std::unique_ptr<Iterator>& it) const override {
            assert(coll < data_.size());
            auto it2 = std::make_unique<Iterator>();
            it2->data_ = data_[coll];
            it2->it_ = it2->data_.begin();
            it.reset(it2.release());
            return Status::OK();
        }

        Status iterator(const uint64_t& coll, const std::string& key, std::unique_ptr<Iterator>& it) const override {
            assert(coll < data_.size());
            auto it2 = std::make_unique<Iterator>();
            it2->data_ = data_[coll];
            it2->it_ = it2->data_.lower_bound(key);
            it.reset(it2.release());
            return Status::OK();
        }
    };

    class WriteBatch : public KeyValue::WriteBatch<uint64_t> {
        std::vector<std::map<std::string,std::string>> data_;
        friend class DB;

    public:
        Status put(const uint64_t& coll, const std::string& key, const std::string& value) override {
            assert(coll < data_.size());
            data_[coll][key] = value;
            return Status::OK();
        };
    };

    class DB : public KeyValue::DB<uint64_t,Mem::Reader,Mem::Iterator,Mem::WriteBatch> {
        std::map<std::string,uint64_t> collections_;
        std::vector<std::map<std::string,std::string>> data_;

    public:
        DB(const std::vector<std::string>& collections) {
            for (uint64_t i = 0; i < collections.size(); i++) {
                assert(collections_.find(collections[i]) == collections_.end());
                collections_[collections[i]] = i;
            }
            data_ = std::vector<std::map<std::string,std::string>>(collections_.size());
        }

        Status collection(const std::string& name, uint64_t& coll) const override {
            auto p = collections_.find(name);
            if (p != collections_.end()) {
                coll = p->second;
                return Status::OK();
            }
            return Status::NotFound("KeyValue::Mem::collection", name);
        }

        Status create_collection(const std::string& name) override {
            if (collections_.find(name) != collections_.end()) {
                return Status::Exists("key-value collection already exists", name);
            }
            collections_[name] = data_.size();
            data_.emplace_back();
            return Status::OK();
        }

        Status current(std::unique_ptr<Mem::Reader>& reader) const override {
            auto p = std::make_unique<KeyValue::Mem::Reader>();
            p->data_ = data_;
            reader = std::move(p);
            return Status::OK();
        }

        Status begin_writes(std::unique_ptr<Mem::WriteBatch>& writes) override {
            auto p = std::make_unique<Mem::WriteBatch>();
            p->data_ = std::vector<std::map<std::string,std::string>>(data_.size());
            writes = std::move(p);
            return Status::OK();
        }

        Status commit_writes(Mem::WriteBatch* writes) override {
            assert(writes != nullptr);
            assert(writes->data_.size() <= data_.size());
            for (size_t i = 0; i < writes->data_.size(); i++) {
                for (const auto& p : writes->data_[i]) {
                    data_[i][p.first] = p.second;
                }
            }
            return Status::OK();
        }

        void wipe() {
            collections_.clear();
            data_.clear();
        }
    };
}

}}

#endif

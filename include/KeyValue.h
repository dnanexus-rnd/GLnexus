#ifndef GLNEXUS_KEYVALUE_H
#define GLNEXUS_KEYVALUE_H
#include "data.h"

namespace GLnexus {
namespace KeyValue {

/// Abstract interface to a database underlying BCFKeyValueData. The database
/// has one or more collections of key-value records. Each collection is
/// ordered by key.

using CollectionHandle = void*;

/// In-order iterator over records in a collection
class Iterator {
public:
    /// Update key and value with the next record and return OK. Return
    /// NotFound if there are no remaining records, or any error code
    virtual Status next(std::string& key, std::string& value) = 0;
};

/// A DB snapshot providing consistent multiple reads if possible
class Reader {
public:
    /// Get the value corresponding to the key and return OK. Return NotFound
    /// if no corresponding record exists in the collection, or any error code
    virtual Status get(CollectionHandle coll, const std::string& key, std::string& value) const = 0;

    /// Create an iterator over the whole collection.
    virtual Status iterator(CollectionHandle coll, std::unique_ptr<Iterator>& it) const = 0;

    /// Create an iterator beginning at the key if a corresponding record
    /// exists, or the first subsequent record otherwise
    virtual Status iterator(CollectionHandle coll, const std::string& key, std::unique_ptr<Iterator>& it) const = 0;
};

/// A batch of writes to apply atomically if possible
class WriteBatch {
public:
    virtual Status put(CollectionHandle coll, const std::string& key, const std::string& value) = 0;
    //virtual Status delete(Collection* coll, const std::string& key) = 0;

    /// Apply a batch of writes.
    virtual Status commit() = 0;
};

/// Main database interface for retrieving collection handles, generating
/// snapshopts to read from, and creating and applying write batches. The DB
/// object itself implements the Reader interface (with no consistency
/// guarantees between multiple calls) and the WriteBatch interface (which
/// applies one write immediately, no atomicity guarantees between multiple
/// calls). Caller must ensure that the parent DB object still exists when any
/// Reader or WriteBatch object is used.
class DB : public Reader {
public:    
    /// Get the handle to a collection, or return NotFound.
    virtual Status collection(const std::string& name, CollectionHandle& coll) const = 0;

    /// Create a new collection, or return Exists.
    virtual Status create_collection(const std::string& name) = 0;

    /// Get an up-to-date snapshot.
    virtual Status current(std::unique_ptr<Reader>& snapshot) const = 0;

    /// Begin preparing a batch of writes.
    virtual Status begin_writes(std::unique_ptr<WriteBatch>& writes) = 0;

    // Base implementations of Reader and WriteBatch interfaces. They simply
    // create a snapshot just to read one record (or begin one iterator), or
    // apply a "batch" of one write. Derived classes may want to provide more
    // efficient overrides.
    Status get(CollectionHandle coll, const std::string& key, std::string& value) const override;
    Status iterator(CollectionHandle coll, std::unique_ptr<Iterator>& it) const override;
    Status iterator(CollectionHandle coll, const std::string& key, std::unique_ptr<Iterator>& it) const override;
    Status put(CollectionHandle coll, const std::string& key, const std::string& value);
};

// Trivial in-memory KeyValue::DB implementation used in unit tests.
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

    class Reader : public KeyValue::Reader {
        std::vector<std::map<std::string,std::string>> data_;
        friend class DB;

    public:
        Status get(CollectionHandle _coll, const std::string& key, std::string& value) const override {
            auto coll = reinterpret_cast<uint64_t>(_coll);
            assert(coll < data_.size());
            const auto& m = data_[coll];
            auto p = m.find(key);
            if (p == m.end()) return Status::NotFound("key", key);
            value = p->second;
            return Status::OK();
        }

        Status iterator(CollectionHandle _coll, std::unique_ptr<KeyValue::Iterator>& it) const override {
            auto coll = reinterpret_cast<uint64_t>(_coll);
            assert(coll < data_.size());
            auto it2 = std::make_unique<Iterator>();
            it2->data_ = data_[coll];
            it2->it_ = it2->data_.begin();
            it.reset(it2.release());
            return Status::OK();
        }

        Status iterator(CollectionHandle _coll, const std::string& key, std::unique_ptr<KeyValue::Iterator>& it) const override {
            auto coll = reinterpret_cast<uint64_t>(_coll);
            assert(coll < data_.size());
            auto it2 = std::make_unique<Iterator>();
            it2->data_ = data_[coll];
            it2->it_ = it2->data_.lower_bound(key);
            it.reset(it2.release());
            return Status::OK();
        }
    };

    class DB;

    class WriteBatch : public KeyValue::WriteBatch {
        std::vector<std::map<std::string,std::string>> data_;
        DB* db_;
        friend class DB;

    public:
        Status put(CollectionHandle _coll, const std::string& key, const std::string& value) override;
        Status commit() override;
    };

    class DB : public KeyValue::DB {
        std::map<std::string,uint64_t> collections_;
        std::vector<std::map<std::string,std::string>> data_;
        friend class WriteBatch;

    public:
        DB(const std::vector<std::string>& collections) {
            for (uint64_t i = 0; i < collections.size(); i++) {
                assert(collections_.find(collections[i]) == collections_.end());
                collections_[collections[i]] = i;
            }
            data_ = std::vector<std::map<std::string,std::string>>(collections_.size());
        }

        Status collection(const std::string& name, CollectionHandle& coll) const override {
            auto p = collections_.find(name);
            if (p != collections_.end()) {
                coll = reinterpret_cast<CollectionHandle>(p->second);
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

        Status current(std::unique_ptr<KeyValue::Reader>& reader) const override {
            auto p = std::make_unique<KeyValue::Mem::Reader>();
            p->data_ = data_;
            reader = std::move(p);
            return Status::OK();
        }

        Status begin_writes(std::unique_ptr<KeyValue::WriteBatch>& writes) override {
            auto p = std::make_unique<Mem::WriteBatch>();
            p->db_ = this;
            p->data_ = std::vector<std::map<std::string,std::string>>(data_.size());
            writes = std::move(p);
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

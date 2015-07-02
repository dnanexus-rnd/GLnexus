#include "KeyValue.h"
#include <assert.h>
using namespace std;

namespace GLnexus { namespace KeyValue {

    // KeyValue::DB base implementations of Reader and WriteBatch interfaces.
    // They simply create a snapshot just to read one record (or begin one
    // iterator), or apply a "batch" of one write. Derived classes may want to
    // provide more efficient overrides.

    template<typename CollectionHandle, class ReaderImpl, class IteratorImpl, class WriteBatchImpl>
    Status DB<CollectionHandle,ReaderImpl,IteratorImpl,WriteBatchImpl>::get(CollectionHandle& coll, const std::string& key, std::string& value) const {
        Status s;
        unique_ptr<ReaderImpl> curr;
        S(current(curr));
        return curr->get(coll, key, value);
    }

    template<typename CollectionHandle, class ReaderImpl, class IteratorImpl, class WriteBatchImpl>
    Status DB<CollectionHandle,ReaderImpl,IteratorImpl,WriteBatchImpl>::iterator(CollectionHandle& coll, unique_ptr<IteratorImpl>& it) const {
        Status s;
        unique_ptr<ReaderImpl> curr;
        S(current(curr));
        return curr->iterator(coll, it);
    }

    template<typename CollectionHandle, class ReaderImpl, class IteratorImpl, class WriteBatchImpl>
    Status DB<CollectionHandle,ReaderImpl,IteratorImpl,WriteBatchImpl>::iterator(CollectionHandle& coll, const string& key, unique_ptr<IteratorImpl>& it) const {
        Status s;
        unique_ptr<ReaderImpl> curr;
        S(current(curr));
        return curr->iterator(coll, key, it);
    }

    template<typename CollectionHandle, class ReaderImpl, class IteratorImpl, class WriteBatchImpl>
    Status DB<CollectionHandle,ReaderImpl,IteratorImpl,WriteBatchImpl>::put(CollectionHandle& coll, const std::string& key, const std::string& value) {
        Status s;
        unique_ptr<WriteBatchImpl> batch;
        S(begin_writes(batch));
        assert(batch);
        S(batch->put(coll, key, value));
        return commit_writes(batch.get());
    }

    // instantiate these template member functions for the test in-memory implementation
    template Status DB<uint64_t,Mem::Reader,Mem::Iterator,Mem::WriteBatch>::get(uint64_t& coll, const std::string& key, std::string& value) const;
    template Status DB<uint64_t,Mem::Reader,Mem::Iterator,Mem::WriteBatch>::iterator(uint64_t& coll, unique_ptr<Mem::Iterator>& it) const;
    template Status DB<uint64_t,Mem::Reader,Mem::Iterator,Mem::WriteBatch>::iterator(uint64_t& coll, const std::string& key, unique_ptr<Mem::Iterator>& it) const;
    template Status DB<uint64_t,Mem::Reader,Mem::Iterator,Mem::WriteBatch>::put(uint64_t& coll, const std::string& key, const std::string& value);


}}

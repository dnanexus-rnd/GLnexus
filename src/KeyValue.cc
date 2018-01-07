#include "KeyValue.h"
#include <assert.h>
using namespace std;

namespace GLnexus { namespace KeyValue {

    // KeyValue::DB base implementations of Reader and WriteBatch interfaces.
    // They simply create a snapshot just to read one record (or begin one
    // iterator), or apply a "batch" of one write. Derived classes may want to
    // provide more efficient overrides.

    Status DB::get0(CollectionHandle coll, const std::string& key, std::shared_ptr<Data>& value) const {
        Status s;
        unique_ptr<Reader> curr;
        S(current(curr));
        return curr->get0(coll, key, value);
    }

    Status DB::iterator(CollectionHandle coll, const string& key, unique_ptr<Iterator>& it) const {
        Status s;
        unique_ptr<Reader> curr;
        S(current(curr));
        return curr->iterator(coll, key, it);
    }

    Status DB::put(CollectionHandle coll, const std::string& key, const KeyValue::Data& value) {
        Status s;
        unique_ptr<WriteBatch> batch;
        S(begin_writes(batch));
        assert(batch);
        S(batch->put(coll, key, value));
        return batch->commit();
    }
}}

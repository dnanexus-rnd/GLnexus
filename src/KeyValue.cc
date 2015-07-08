#include "KeyValue.h"
#include <assert.h>
using namespace std;

namespace GLnexus { namespace KeyValue {

    // KeyValue::DB base implementations of Reader and WriteBatch interfaces.
    // They simply create a snapshot just to read one record (or begin one
    // iterator), or apply a "batch" of one write. Derived classes may want to
    // provide more efficient overrides.

    Status DB::get(CollectionHandle coll, const std::string& key, std::string& value) const {
        Status s;
        unique_ptr<Reader> curr;
        S(current(curr));
        return curr->get(coll, key, value);
    }

    Status DB::iterator(CollectionHandle coll, unique_ptr<Iterator>& it) const {
        Status s;
        unique_ptr<Reader> curr;
        S(current(curr));
        return curr->iterator(coll, it);
    }

    Status DB::iterator(CollectionHandle coll, const string& key, unique_ptr<Iterator>& it) const {
        Status s;
        unique_ptr<Reader> curr;
        S(current(curr));
        return curr->iterator(coll, key, it);
    }

    Status DB::put(CollectionHandle coll, const std::string& key, const std::string& value) {
        Status s;
        unique_ptr<WriteBatch> batch;
        S(begin_writes(batch));
        assert(batch);
        S(batch->put(coll, key, value));
        return batch->commit();
    }

    namespace Mem {
         Status WriteBatch::put(CollectionHandle _coll, const std::string& key, const std::string& value) {
            auto coll = reinterpret_cast<uint64_t>(_coll);
            assert(coll < data_.size());
            data_[coll][key] = value;
            return Status::OK();
        };

        Status WriteBatch::commit() {
            assert(data_.size() <= db_->data_.size());
            for (size_t i = 0; i < data_.size(); i++) {
                for (const auto& p : data_[i]) {
                    db_->data_[i][p.first] = p.second;
                }
            }
            return Status::OK();
        }
    }
}}

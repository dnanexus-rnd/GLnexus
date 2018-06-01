@0x926cdc189d296294;

using Cxx = import "/capnp/c++.capnp";
$Cxx.namespace("GLnexus::capnp");

### BCFBucket (for internal database use)
struct BCFBucketSkipEntry {
    recordIndex @0 : UInt32;
    posBeg @1 : Int64;
}
struct BCFBucket {
    records @0 : List(Data);
    skips @1 : List(BCFBucketSkipEntry);
}

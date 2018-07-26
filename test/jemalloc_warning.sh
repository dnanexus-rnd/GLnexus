#!/bin/bash
# test that we can dynamically link jemalloc into glnexus_cli, and its
# built-in warning about whether jemalloc is present
#find /usr -name "libjemalloc.so*"
set -ex -o pipefail

HERE=$(dirname "$0")
LIBJEMALLOC_SO=$(ldconfig -p | grep "libjemalloc.so " | cut -f2 -d '>' | tr -d ' ')
echo "path to libjemalloc.so: $LIBJEMALLOC_SO"
unset LD_PRELOAD

exit_code=0
MALLOC_CONF=bogus:true LD_PRELOAD=$LIBJEMALLOC_SO $HERE/../glnexus_cli -h 2>&1 | grep 'Invalid conf pair' || exit_code=$?
if [ "$exit_code" -ne 0 ]; then
    echo "couldn't confirm jemalloc dynamic linking"
    exit 1
fi
$HERE/../glnexus_cli -h 2>&1 | grep 'jemalloc absent' || exit_code=$?
if [ "$exit_code" -ne 0 ]; then
    echo "false-negative jemalloc warning"
    exit 1
fi
LD_PRELOAD=$LIBJEMALLOC_SO $HERE/../glnexus_cli -h 2>&1 | grep 'jemalloc absent' || exit_code=$?
if [ "$exit_code" -eq 0 ]; then
    echo "false-positive jemalloc warning"
    exit 1
fi

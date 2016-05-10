#!/bin/bash

set -ex -o pipefail

perf_pid=""
dstat_pid=""
iostat_pid=""

function clean_shutdown {
    if [ -n $perf_pid ]; then
        sudo kill -s SIGINT $perf_pid
    fi
    if [ -n $dstat_pid ]; then
        kill -9 $dstat_pid
    fi
    if [ -n $iostat_pid ]; then
        kill -9 $iostat_pid
    fi
}
trap clean_shutdown EXIT
trap clean_shutdown SIGINT

#disable_huge_pages

# Allow kernel tracing
sudo sysctl -w kernel.kptr_restrict=0
sudo sh -c 'echo -1 > /proc/sys/kernel/perf_event_paranoid'

# parameters for performance recording
recordFreq=99

# log detailed utilization
dstat -cmdn 10 &
dstat_pid=$!
iostat -x 60 &
iostat_pid=$!

mkdir -p out/perf
sudo perf record -F $recordFreq -a -g &
perf_pid=$!

time numactl --interleave=all glnexus_cli genotype GLnexus.db --config test -t $(nproc) --bed ranges.bed | bcftools view - | bgzip -c > out/matches.vcf.gz


# Try to kill the perf process nicely; this does not always work
sleep 2
clean_shutdown
sudo chmod 644 perf.data
perf script > out/perf/perf_genotype
stackcollapse-perf.pl < out/perf/perf_genotype > out/perf/genotype.stacks
flamegraph.pl out/perf/genotype.stacks > out/perf/genotype_perf.svg





#!/bin/bash

# Constant for tracing, 99Hz
recordFreq=99

main() {
    set -ex -o pipefail

    compressor="pigz"
    compress_ext="gz"

    # log detailed utilization
    dstat -cmdn 60 &

    # Replace malloc with jemalloc
    export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.1

    # Make sure jemalloc is in use
    check_use_of_jemalloc_lib

    # download inputs.
    dx-download-all-inputs --parallel

    if [ -n "$gvcf_manifest" ]; then
        dx download "$gvcf_manifest" -o gvcf_manifest.txt
        mkdir -p in/gvcf/manifest
        cat gvcf_manifest.txt | xargs -n 100 -P 5 dx download --no-progress -o in/gvcf/manifest
    fi

    # Make a list of all gvcf files
    find in/gvcf -type f > /tmp/gvcf_list
    wc -l /tmp/gvcf_list

    debug_flags=""
    if [[ $bed_ranges_to_genotype ]]; then
        bed_ranges=$(find in/bed_ranges_to_genotype -type f)
    fi
    if [[ $debug == "true" ]]; then
        debug_flags+=" --debug"
    fi
    if [[ $iter_compare == "true" ]]; then
        debug_flags+=" --iter_compare"
    fi
    if [[ "$generate_perf_kernel" == "true" ]]; then
        # Kernel tracing implies user-space space tracing
        generate_perf="true"
    fi
    if [[ "$generate_perf" == "true" ]]; then
        setup_system_for_tracing

        # Run a perf command that will record everything the
        # machine is doing, this allows tracking down all activity,
        mkdir -p out/perf
        sudo perf record -F $recordFreq -a --call-graph dwarf &
        perf_pid=$!
    fi
    bucket_size_arg=""
    if [ -n "$bucket_size" ]; then
        bucket_size_arg="--bucket_size $bucket_size"
    fi
    squeeze_arg=""
    squeeze_cmd="cat"
    if [[ "$squeeze" == "true" ]]; then
        squeeze_arg="--squeeze"
        squeeze_cmd="spvcf squeeze -q -"
    fi
    config_arg="--config $config"
    if compgen -G "in/config_yml/*.yml" > /dev/null; then
        config_arg="--config in/config_yml/*.yml"
    fi

    mkdir -p out/vcf
    time numactl --interleave=all glnexus_cli $config_arg $squeeze_arg --list --bed $bed_ranges $bucket_size_arg $debug_flags /tmp/gvcf_list \
        | bcftools view - | $squeeze_cmd | bgzip --threads $(nproc) -c > "out/vcf/${output_name}.vcf.${compress_ext}"

    if [[ "$perf" == "true" ]]; then
        # Try to kill the perf process nicely; this does not always work
        sleep 2
        sudo kill -s SIGINT $perf_pid
        sudo chmod 644 perf.data
        perf script > perf_genotype
        FlameGraph/stackcollapse-perf.pl < perf_genotype > out/perf/genotype.stacks
        /bin/rm -f perf.data
    fi

    # we are writing the generated VCF to stdout, so the residuals will
    # be placed in a default location.
    if [ -f /tmp/residuals.yml ]; then
        mkdir -p out/residuals/
        $compressor -c /tmp/residuals.yml > "out/residuals/${output_name}.residuals.yml.${compress_ext}"
    fi

    # discovered alleles file
    if [ -f /tmp/dsals.yml ]; then
        mkdir -p out/discovered_alleles
        $compressor -c /tmp/dsals.yml > out/discovered_alleles/"${output_name}.discovered_alleles.yml.${compress_ext}"
    fi
    if [ -f /tmp/dsals.cflat ]; then
        mkdir -p out/discovered_alleles
        $compressor -c /tmp/dsals.cflat > out/discovered_alleles/"${output_name}.discovered_alleles.cflat.${compress_ext}"
    fi

    # unified alleles file
    if [ -f /tmp/sites.yml ]; then
        mkdir -p out/unified_sites
        $compressor -c /tmp/sites.yml > "out/unified_sites/${output_name}.sites.yml.${compress_ext}"
    fi

    ls -Rlh GLnexus.DB
    ls -Rlh out

    # upload
    dx-upload-all-outputs --parallel

    # sleep if requested
    if [ "$sleep" -gt "0" ]; then
        sleep "$sleep"
    fi
}

# Make sure the jemalloc library is used.
function check_use_of_jemalloc_lib {
    MALLOC_CONF=invalid_flag:foo glnexus_cli >& /tmp/malloc_dbg || true
    num_appear=$(grep jemalloc /tmp/malloc_dbg | wc -l)
    if [[ $num_appear == "0" ]]; then
        echo "Error, jemalloc is supposed to be in use"
        exit 1
    fi
}

function setup_system_for_tracing {
    # We want to have visibility into kernel symbols. This is under a
    # special flag, because normally, we do not have sufficient
    # permissions in a platform container. We do this first, so, if there
    # are any permission problems, we will fail early.
    if [[ "$perf_kernel" == "true" ]]; then
        apt_get_add_debug_repos
        install_kernel_debug_symbols
    fi

    # install the perf utility
    # install flame-graph package
    git clone https://github.com/brendangregg/FlameGraph

    # install linux-tools for the current kernel version
    linux_version=`uname -r`
    sudo apt-get -y -qq install "linux-tools-${linux_version}"

    # test that perf is working correctly
    perf record -F $recordFreq -g /bin/ls
    rm -f perf.data
}


# Add special debian repositories holding kernel debugging symbols, and
# libraries with debug information.
function apt_get_add_debug_repos {
    echo "Adding apt-get repositories with debug symbols"
    lsb_rel=$(lsb_release -cs)
    echo "deb http://ddebs.ubuntu.com ${lsb_rel} main restricted universe multiverse" > ddebs.list
    echo "deb http://ddebs.ubuntu.com ${lsb_rel}-updates main restricted universe multiverse" >> ddebs.list
    echo "deb http://ddebs.ubuntu.com ${lsb_rel}-proposed main restricted universe multiverse" >> ddebs.list
    sudo cp ddebs.list /etc/apt/sources.list.d/

    sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys C8CAB6595FDFF622
    sudo apt-get update
}

function install_kernel_debug_symbols {
    echo "Setting Linux permissions to allow seeing kernel symbols"
    sudo sysctl -w kernel.kptr_restrict=0
    sudo sh -c 'echo -1 > /proc/sys/kernel/perf_event_paranoid'

    echo "Installing kernel debug symbols"
    # Check for an exact match to the kernel version
    kernel_version=$(uname -r)
    retval=$(apt-cache search linux-image-$kernel_version-dbgsym)
    if [[ -n "$retval" ]]; then
        sudo apt-get install linux-image-$kernel_version-dbgsym
        return
    fi

    echo "An exact match for kernel [$kernel_version] debug symbols is not available"
    kernel_version_short=$(uname -r | cut --delimiter='-' --fields=1)
    echo "Searching for close match for kernel $kernel_version_short"

    pkg=$(apt-cache search linux-image-$kernel_version_short | grep generic-dbgsym | head -n 1 | cut --delimiter=' ' --fields=1)
    echo "Found package $pkg, installing ..."
    sudo apt-get install $pkg
}

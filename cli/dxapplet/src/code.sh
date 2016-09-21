#!/bin/bash

# http://stackoverflow.com/a/17841619
function join { local IFS="$1"; shift; echo "$*"; }

main() {
    set -ex -o pipefail

    # parameters for performance recording
    recordFreq=99

    # log detailed utilization
    dstat -cmdn 60 &
    iostat -x 600 &

    # Replace malloc with jemalloc
    export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.1

    # Make sure jemalloc is in use
    check_use_of_jemalloc_lib

    # Kernel tracing implies user-space space tracing
    if [ "$enable_kernel_perf" == "true" ]; then
        enable_perf=true
    fi

    # We want to have visibility into kernel symbols. This is under a
    # special flag, because normally, we do not have the right
    # permissions in a platform container. We do this first, so, if there
    # are any permission problems, we will fail early.
    if [ "$enable_kernel_perf" == "true" ]; then
        apt_get_add_debug_repos
        install_kernel_debug_symbols
    fi

    # install the perf utility
    if [ "$enable_perf" == "true" ]; then
        # install flame-graph package
        git clone https://github.com/brendangregg/FlameGraph

        # install linux-tools for the current kernel version
        linux_version=`uname -r`
        sudo apt-get -y -qq install "linux-tools-${linux_version}"

        # test that perf is working correctly
        perf record -F $recordFreq -g /bin/ls
        rm -f perf.data
    fi

    # download inputs.
    # TODO: it would be nice to overlap the staging of the input gVCFs with
    # the bulk load, by streaming the filenames into glnexus_cli as they're
    # staged.
    mkdir -p "in/gvcf" "in/gvcf_tar"
    dx-download-all-inputs --parallel --except gvcf_tar --except existing_db
    for tar_dxlink in "${gvcf_tar[@]}"
    do
        tar_dxid=$(dx-jobutil-parse-link --no-project "$tar_dxlink")
        dn="in/gvcf/${tar_dxid}"
        mkdir -p "$dn"
        dx cat "$tar_dxid" | tar x -C "$dn" --strip-components=1
    done
    find in/gvcf -type f > all_gvcfs.txt
    wc -l all_gvcfs.txt

    # initialize database if none provided
    if [ -n "$existing_db" ]; then
        dx cat "$existing_db" | tar x
    else
        bucket_size_arg=""
        if [ -n "$bucket_size" ]; then
            bucket_size_arg="--bucket-size $bucket_size"
        fi
        glnexus_cli init $bucket_size_arg GLnexus.db $(find in/gvcf -type f | head -n 1)
    fi

    # load gVCFs
    ranges_to_load_arg=""
    if [ -n "$ranges_to_load" ]; then
        ranges_to_load_arg=$(join , "${ranges_to_load[@]}")
        ranges_to_load_arg="--range $ranges_to_load_arg"
    fi
    cat all_gvcfs.txt | time glnexus_cli load $ranges_to_load_arg --and-delete GLnexus.db -
    ls -lh GLnexus.db
    mkdir -p out/db_load_log
    cp GLnexus.db/LOG "out/db_load_log/${output_name}.LOG"


    # Test that the iterators work correctly
    if [ "$iter_compare" == "true" ]; then
        echo "Comparing iterator implementations"
        glnexus_cli iter_compare GLnexus.db
    fi

    # assemble BED file of ranges to genotype
    mkdir -p in/bed_ranges_to_genotype
    touch in/bed_ranges_to_genotype/DUMMY
    cat in/bed_ranges_to_genotype/* > ranges.bed
    for range in "${ranges_to_genotype[@]}"; do
        range_sp=$(echo "$range" | tr -d "," | tr ":-" "\t")
        echo "$range_sp" >> ranges.bed
    done
    wc -l ranges.bed

    # genotype the ranges
    if [ "$(cat ranges.bed | wc -l)" -gt "0" ]; then
        if [ "$enable_perf" == "true" ]; then
            # Run a perf command that will record everything the
            # machine is doing, this allows tracking down all activity,
            mkdir -p out/perf
            sudo perf record -F $recordFreq -a --call-graph dwarf &
            perf_pid=$!
        fi

        residuals_flag=""
        if [[ $residuals == "true" ]]; then
            residuals_flag="--residuals"
        fi

        config_flag=""
        if [ -n "$config" ]; then
            config_flag="--config $config"
        fi

        # set up for gzip or lz4 output compression
        compressor="pigz"
        vcf_compressor="bgzip --threads $(nproc)"
        compress_ext="gz"
        if [ "$lz4_output" == "true" ]; then
            compressor="lz4"
            vcf_compressor="lz4"
            compress_ext="lz4"
        fi

        # numactl explanation: https://blog.jcole.us/2010/09/28/mysql-swap-insanity-and-the-numa-architecture/
        time numactl --interleave=all glnexus_cli discover_alleles GLnexus.db --bed ranges.bed $config_flag --dsals discovered_alleles.cflat
        ls -l discovered_alleles.cflat

        time numactl --interleave=all glnexus_cli unify_sites discovered_alleles.cflat --bed ranges.bed $config_flag > unified_sites.yml

        mkdir -p out/vcf
        time numactl --interleave=all glnexus_cli genotype GLnexus.db unified_sites.yml $residuals_flag $config_flag | bcftools view - | $vcf_compressor -c > "out/vcf/${output_name}.vcf.${compress_ext}"

        # we are writing the generated VCF to stdout, so the residuals will
        # be placed in a default location.
        if [ -f /tmp/residuals.yml ]; then
            mkdir -p out/residuals/
            $compressor -c /tmp/residuals.yml > "out/residuals/${output_name}.residuals.yml.${compress_ext}"
        fi

        if [ "$enable_perf" == "true" ]; then
            # Try to kill the perf process nicely; this does not always work
            sleep 2
            sudo kill -s SIGINT $perf_pid
            sudo chmod 644 perf.data
            perf script > perf_genotype
            FlameGraph/stackcollapse-perf.pl < perf_genotype > out/perf/genotype.stacks
            /bin/rm -f perf.data
        fi

        # compress the discovered alleles file, and place it in the output directory
        mkdir -p out/discovered_alleles
        $compressor -c discovered_alleles.cflat > out/discovered_alleles/"${output_name}.discovered_alleles.cflat.${compress_ext}"

        # compress unified sites file, place in output directory
        mkdir -p out/unified_sites
        $compressor -c unified_sites.yml > out/unified_sites/"${output_name}.unified_sites.yml.${compress_ext}"
    fi

    # upload
    dx-upload-all-outputs --parallel

    # upload database (unless we used an existing one that we didn't change)
    if [ -z "$existing_db" ] || [ "$(cat all_gvcfs.txt | wc -l)" -gt "0" ]; then
        dxdb=$(tar c GLnexus.db | dx upload --destination "${output_name}.db.tar" --type GLnexus_db --brief -)
        dx-jobutil-add-output db "$dxdb" --class=file
    fi

    # sleep if requested
    if [ "$sleep" -gt "0" ]; then
        sleep "$sleep"
    fi
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

# Make sure the jemalloc library is used.
function check_use_of_jemalloc_lib {
    MALLOC_CONF=invalid_flag:foo glnexus_cli >& /tmp/malloc_dbg || true
    num_appear=$(grep jemalloc /tmp/malloc_dbg | wc -l)
    if [[ $num_appear == "0" ]]; then
        echo "Error, jemalloc is supposed to be in use"
        exit 1
    fi
}

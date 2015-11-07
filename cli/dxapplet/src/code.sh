#!/bin/bash

main() {
    set -ex -o pipefail

    # parameters for performance recording
    recordFreq=99

    # log detailed utilization
    dstat -cmdn 10 &

    # install the perf utility
    if [ "$enable_perf" == "1" ]; then
        ubuntu_release=$(lsb_release -r | cut -f 2)

        case "$ubuntu_release" in
            "12.04")
                linux_version=`uname -r`
                sudo apt-get -y install linux-tools-${linux_version}
                sudo apt-get -y install linux-tools
                ;;
            "14.04")
                sudo apt-get -y install linux-tools-common linux-tools-generic
                ;;
            *)
                echo "Error, ubuntu release is not 12.04 or 14.04, can't install perf"
                exit 1
        esac

        # test that perf is working correctly
        perf record -F $recordFreq -g /bin/ls
        rm -f perf.data

        # install flame-graph package
        git clone https://github.com/brendangregg/FlameGraph
    fi

    # install dependencies
    sudo rm -f /etc/apt/apt.conf.d/99dnanexus
    sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
    sudo apt-get -qq update
    sudo apt-get -qq install -y gcc-4.9 g++-4.9 liblz4-dev
    sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 100 \
                             --slave /usr/bin/g++ g++ /usr/bin/g++-4.9

    # download inputs
    dx-download-all-inputs --parallel --except gvcf_tar --except existing_db
    mkdir -p in/gvcf
    if [ -n "$gvcf_tar" ]; then
        dx cat "$gvcf_tar" | tar x -C in/gvcf --strip-components=1
    fi
    find in/gvcf -type f > all_gvcfs.txt
    wc -l all_gvcfs.txt

    # initialize database if none provided
    if [ -n "$existing_db" ]; then
        dx cat "$existing_db" | tar x
    else
        glnexus_cli init GLnexus.db $(find in/gvcf -type f | head -n 1)
    fi

    # load gVCFs
    cat all_gvcfs.txt | time glnexus_cli load --and-delete GLnexus.db -
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
        if [ "$enable_perf" == "1" ]; then
            mkdir -p out/perf
            sudo perf record -F $recordFreq -a -g &
        fi

        mkdir -p out/vcf
        time glnexus_cli genotype GLnexus.db --bed ranges.bed \
            | bcftools view - | bgzip -c > "out/vcf/${output_name}.vcf.gz"

        if [ "$enable_perf" == "1" ]; then
            sudo killall perf
            sudo chmod 644 perf.data
            perf script > perf_genotype
            FlameGraph/stackcollapse-perf.pl < perf_genotype > out/perf/genotype.stacks
            mv out/perf/genotype.stacks out/perf/genotype.stacks.${range}
            /bin/rm -f perf.data
        fi
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

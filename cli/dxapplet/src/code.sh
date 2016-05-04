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

#    sudo bash -c 'echo "0" > /proc/sys/kernel/kptr_restrict'

    # install the perf utility
    if [ "$enable_perf" == "true" ]; then
        # install linux-tools for the current kernel version
        linux_version=`uname -r`
        sudo apt-get -y -qq install "linux-tools-${linux_version}"

        kernel_dbgsym_pkg=""
        case $linux_version in
            "3.13.0-85-generic")
                kernel_dbgsym_pkg=linux-image-3.13.0-85-generic-dbgsym_3.13.0-85.129_amd64.ddeb
                ;;
            "3.13.0-74-generic")
                kernel_dbgsym_pkg=linux-image-3.13.0-74-generic-dbgsym_3.13.0-74.118_i386.ddeb
                ;;
            *)

                echo "We don't have the necessary debugging symbols for $linux_version"
                exit 1
                ;;
        esac

        echo "installing kernel debug symbols from project resources"
        sudo dpkg -i $kernel_dbgsym_pkg
        sudo apt-get install -f

        echo "print all the loadable kernel modules"
        lsmod

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
            mkdir -p out/perf
            sudo perf record -F $recordFreq -a -g &
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

        mkdir -p out/vcf
        time glnexus_cli genotype GLnexus.db $residuals_flag $config_flag -t $(nproc) --bed ranges.bed \
            | bcftools view - | bgzip -c > "out/vcf/${output_name}.vcf.gz"

        # we are writing the generated VCF to stdout, so the residuals will
        # be placed in a default location.
        if [ -f /tmp/residuals.yml ]; then
            mkdir -p out/residuals/
            mv /tmp/residuals.yml out/residuals/
        fi

        if [ "$enable_perf" == "true" ]; then
            sleep 5
            sudo kill -s SIGINT $perf_pid
            sudo chmod 644 perf.data
            perf script > perf_genotype
            FlameGraph/stackcollapse-perf.pl < perf_genotype > out/perf/genotype.stacks
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

#!/bin/bash

main() {
    set -ex -o pipefail

    # parameters for performance recording
    recordFreq=99

    # log detailed utilization
    dstat -cmdn 10 &

    # install the perf utility
    if [ "$enable_perf" == "1" ]; then
        linux_version=`uname -r`
        sudo apt-get -y install linux-tools-${linux_version}
        sudo apt-get -y install linux-tools

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
    sudo apt-get -qq install -y gcc-4.9 g++-4.9
    sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 100 \
                             --slave /usr/bin/g++ g++ /usr/bin/g++-4.9
    LZ4_REVISION="r131"
    tar zxf lz4-${LZ4_REVISION}.tar.gz
    make -C lz4-${LZ4_REVISION} -j$(nproc)
    sudo make -C lz4-${LZ4_REVISION} install
    sudo ldconfig

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
    cat all_gvcfs.txt | time glnexus_cli load GLnexus.db -
    ls -lh GLnexus.db
    mkdir -p out/db_load_log
    cp GLnexus.db/LOG "out/db_load_log/${output_name}.LOG"

    # genotype specified ranges
    if [ "${#ranges_to_genotype[@]}" -gt "0" ]; then
        mkdir bcf
        i=0
        for range in "${ranges_to_genotype[@]}"; do
            if [ "$enable_perf" == "1" ]; then
                mkdir -p out/perf
                sudo perf record -F $recordFreq -a -g &
            fi

            range_sp=$(echo "$range" | tr -d "," | tr ":-" " ")
            outfn=$(printf "bcf/%05d.bcf" "$i")
            time glnexus_cli genotype GLnexus.db $range_sp > "$outfn"
            i=$(expr $i + 1)

            if [ "$enable_perf" == "1" ]; then
                sudo killall perf
                sudo chmod 644 perf.data
                perf script > perf_genotype
                FlameGraph/stackcollapse-perf.pl < perf_genotype > out/perf/genotype.stacks
                mv out/perf/genotype.stacks out/perf/genotype.stacks.${range}
                /bin/rm -f perf.data
            fi
        done
        ls -lh bcf

        # generate merged VCF
        mkdir -p out/vcf
        if [ "${#ranges_to_genotype[@]}" -gt "1" ]; then
            find bcf/ -type f -name "*.bcf" | xargs -n 1 -t bcftools index
            bcftools concat bcf/*.bcf | bgzip -c > "out/vcf/${output_name}.vcf.gz"
        else
            bcftools view bcf/*.bcf | bgzip -c > "out/vcf/${output_name}.vcf.gz"
        fi
        ls -lh out/vcf
    fi

    # upload
    dx-upload-all-outputs --parallel

    # upload database
    dxdb=$(tar c GLnexus.db | dx upload --destination "${output_name}.db.tar" --type GLnexus_db --brief -)
    dx-jobutil-add-output db "$dxdb" --class=file

    # sleep if requested
    if [ "$sleep" -gt "0" ]; then
        sleep "$sleep"
    fi
}

#!/bin/bash

main() {
    set -ex -o pipefail

    # log detailed utilization
    dstat -cmdn 10 &

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
    dx-download-all-inputs --parallel --except gvcf_tar
    if [ -n "$gvcf_tar" ]; then
        mkdir -p in/gvcf
        dx cat "$gvcf_tar" | tar x -C in/gvcf --strip-components=1
    fi
    find in/gvcf -type f > all_gvcfs.txt
    wc -l all_gvcfs.txt

    # initialize and load database
    glnexus_cli init GLnexus.db $(find in/gvcf -type f | head -n 1)
    cat all_gvcfs.txt | time glnexus_cli load GLnexus.db -
    ls -lh GLnexus.db
    mkdir -p out/db_load_log
    cp GLnexus.db/LOG "out/db_load_log/${output_name}.LOG"

    # genotype specified ranges
    if [ "${#ranges_to_genotype[@]}" -gt "0" ]; then
        mkdir bcf
        i=0
        for range in "${ranges_to_genotype[@]}"; do
            range_sp=$(echo "$range" | tr -d "," | tr ":-" " ")
            outfn=$(printf "bcf/%05d.bcf" "$i")
            time glnexus_cli genotype --cache-bytes "$cache_bytes" GLnexus.db $range_sp > "$outfn"
            i=$(expr $i + 1)
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

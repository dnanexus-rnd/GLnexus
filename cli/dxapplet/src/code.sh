#!/bin/bash

main() {
    set -ex -o pipefail

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
    dx-download-all-inputs --parallel

    # initialize and load database
    glnexus_cli init GLnexus.db $(find in/gvcf -type f | head -n 1)
    date
    find in/gvcf -type f | xargs -n 1 -t time glnexus_cli load GLnexus.db
    ls -lh GLnexus.db
    date

    # genotype specified ranges
    if [ "${#ranges_to_genotype[@]}" -gt "0" ]; then
        mkdir bcf
        i=0
        for range in "${ranges_to_genotype[@]}"; do
            range_sp=$(echo "$range" | tr -d "," | tr ":-" " ")
            outfn=$(printf "bcf/%05d.bcf" "$i")
            time glnexus_cli genotype GLnexus.db $range_sp > "$outfn"
            i=$(expr $i + 1)
        done
        ls -lh bcf

        # generate merged VCF
        mkdir -p out/vcf
        if [ "${#ranges_to_genotype[@]}" -gt "1" ]; then
            bcftools merge bcf/*.bcf | bgzip -c > "out/vcf/${output_name}.vcf.gz"
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

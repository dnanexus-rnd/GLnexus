#!/bin/bash

main() {
    set -ex -o pipefail

    compressor="pigz"
    vcf_compressor="bgzip --threads $(nproc)"
    compress_ext="gz"

    # log detailed utilization
    dstat -cmdn 20 &

    # Replace malloc with jemalloc
    export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.1

    # Make sure jemalloc is in use
    check_use_of_jemalloc_lib

    # download inputs.
    # TODO: it would be nice to overlap the staging of the input gVCFs with
    # the bulk load, by streaming the filenames into glnexus_cli as they're
    # staged.
    mkdir -p "in/gvcf" "in/gvcf_tar"
    dx-download-all-inputs --parallel --except gvcf_tar
    for tar_dxlink in "${gvcf_tar[@]}"
    do
        tar_dxid=$(dx-jobutil-parse-link --no-project "$tar_dxlink")
        dn="in/gvcf/${tar_dxid}"
        mkdir -p "$dn"
        dx cat "$tar_dxid" | tar x -C "$dn" --strip-components=1
    done

    # Make a list of all gvcf files
    gvcfs=$(find in/gvcf -type f)

    if [[ $bed_ranges_to_genotype ]]; then
        bed_ranges=$(find in/bed_ranges_to_genotype -type f)
    fi
    if [[ $residuals == "true" ]]; then
        residuals_flag="--residuals"
    fi
    if [[ $debug == "true" ]]; then
        debug_flag="--debug"
    fi

    mkdir -p out/vcf
    time numactl --interleave=all glnexus_cli --bed $bed_ranges $residuals_flag $debug_flag $gvcfs | bcftools view - | $vcf_compressor -c > "out/vcf/${output_name}.vcf.${compress_ext}"

    # we are writing the generated VCF to stdout, so the residuals will
    # be placed in a default location.
    if [ -f /tmp/residuals.yml ]; then
        mkdir -p out/residuals/
        $compressor -c /tmp/residuals.yml > "out/residuals/${output_name}.residuals.yml.${compress_ext}"
    fi

    # discovered alleles file
    if [ -f /tmp/dsals.yml ]; then
        mkdir -p out/discovered_alleles
        $compressor -c /tmp/dsals.yml > "out/discovered_alleles/${output_name}.dsals.yml.${compress_ext}"
    fi

    # unified alleles file
    if [ -f /tmp/sites.yml ]; then
        mkdir -p out/unified_sites
        $compressor -c /tmp/sites.yml > "out/unified_sites/${output_name}.sites.yml.${compress_ext}"
    fi

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

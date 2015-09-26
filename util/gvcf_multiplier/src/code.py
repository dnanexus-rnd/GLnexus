#!/usr/bin/env python
import dxpy
import subprocess
import multiprocessing
from multiprocessing.pool import ThreadPool as Pool
import os
import sys

@dxpy.entry_point("main")
def main(gvcf,N,sample_name_prefix,output_name):
    K = len(gvcf)

    # download all the source gVCFs
    sh("dx-download-all-inputs --parallel")

    # create output directory
    os.mkdir("gvcf")

    # parallel generate gVCF files
    pool = Pool(multiprocessing.cpu_count())
    inputs = [{
        "source_index": i%K,
        "sample_name_prefix": sample_name_prefix,
        "dest_index": i
    } for i in xrange(N)]
    pool.map(generate_gvcf_kwargs, inputs)

    # tar and upload
    dxid = subprocess.check_output(["/bin/bash", "-e", "-o", "pipefail", "-c", 'tar cv gvcf | dx upload --brief --destination "{}.tar" -'.format(output_name)]).strip()
    return {
        "tar": dxpy.dxlink(dxid)
    }

def generate_gvcf(source_index,sample_name_prefix,dest_index):
    dest_index_str = "{0:06d}".format(dest_index)
    sample_name = "{}{}".format(sample_name_prefix, dest_index_str)
    sh('echo "{}" > "{}.reheader"'.format(sample_name, sample_name))
    sh('bcftools reheader -s "{}.reheader" in/gvcf/{}/* > "gvcf/{}.gvcf.gz"'.format(sample_name, source_index, sample_name))
    os.unlink("{}.reheader".format(sample_name))

def generate_gvcf_kwargs(kwargs):
    return generate_gvcf(**kwargs)

def sh(cmd):
    print cmd
    subprocess.check_call(["/bin/bash", "-e", "-o", "pipefail", "-c", cmd])

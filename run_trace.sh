#!/bin/bash

dx build -f --destination :/Trace/GLnexus cli/dxapplet

dx run :/Trace/GLnexus -i gvcf_tar=:/Fake/gvcfs1024.tar -i bed_ranges_to_genotype=:/vcr_b37.bed -i config=test -i enable_perf=true --folder :/Trace --instance-type mem3_ssd1_x16 -y --watch

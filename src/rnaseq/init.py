#!/usr/bin/env python
# -*- coding: utf-8 -*-

from rnaseq import utility

MANIFEST = """####################################################################################################
# Manifest file for running rnaseq pipeline
#
# This CSV manifest file should contain at least 2 columns (SampleName and FASTQ1)
#  1. SampleName: in the format of group_replicate, e.g., p53OE_Veh_1, p53OE_Veh_2, p53OE_Veh_3
#  2. FASTQ1: absolute path to the forward strand (read 1) FASTQ file
#  3. FASTQ2: absolute path to the reverse strand (read 2) FASTQ file, optional for pair-end data
#
# Note:
#  a. The names in SampleName column must have at least 1 underscore (_) to be used for
#     right splitting, thus we can get group and replicate of each sample
#  b. For single-end data, you should either delete ",FASTQ2" in the header line or
#     add a comma (,) after each path of FASTQ1 file
# Any line (include this line) starts with '#' will be ignore and thus can be deleted
######################################################################################################
SampleName,FASTQ1,FASTQ2
"""


def init(args):
    outdir = utility.mkdir(args)
    manifest = outdir.joinpath('manifest.csv')
    
    if manifest.exists():
        utility.logger.info(f'Manifest file already exists, skip recreating it')
    else:
        if args.dry:
            utility.logger.info(f'Manifest file not found, will create {manifest}')
        else:
            utility.logger.info(f'Manifest file not found, creating {manifest}')
            manifest.parent.mkdir(exist_ok=True, parents=True)
            with manifest.open('w') as o:
                o.write(MANIFEST)
            utility.logger.info(f'Manifest file created: {manifest}')

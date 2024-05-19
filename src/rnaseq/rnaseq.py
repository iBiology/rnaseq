#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A toolkit for processing RNA-Seq data to identify gene differential expression.
"""

import argparse
from pathlib import Path

import cmder
from rnaseq import utility, init, pipeline, dega, fastp


def main():
    parser = argparse.ArgumentParser(description=__doc__, prog='rnaseq')
    
    parser.add_argument('-c', '--comparison', nargs='?',
                        help='Comparisons (in the format of a.vs.c, b.vs.c, ...) needed to '
                             'run gene differential expression analysis')
    parser.add_argument('-o', '--outdir',
                        help=f'Path to a directory for writing output files, default: %(default)s',
                        default=Path.cwd(), type=Path)
    parser.add_argument('-p', '--process', help=f'Number of CPUs to use, default: %(default)s',
                        default=16, type=int)
    parser.add_argument('-d', '--dry', action='store_true',
                        help='Print out steps and files involved in each step without actually processing data.')
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {utility.VERSION}')
    
    subparsers = parser.add_subparsers(help='', title='subcommands')
    
    description = 'Submit RNA-Seq job to run RNA-Seq pipeline'
    submit_parser = subparsers.add_parser('submit', help=description, description=description)
    
    description = 'Run RNA-Seq pipeline'
    pipeline_parser = subparsers.add_parser('run', help=description, description=description)
    pipeline_parser.add_argument('manifest', help=f'Path to a CSV manifest file', type=cmder.File)
    pipeline_parser.add_argument('--genome', help='Path to a directory contains STAR index files',
                                 type=cmder.Dir, required=True)
    pipeline_parser.add_argument('--gtf', help='Path to a GTF file contains genomic annotation',
                                 type=cmder.File, required=True)
    pipeline_parser.add_argument('--bed', help='Path to a BED file contains genomic annotation',
                                 type=cmder.File, required=True)
    pipeline_parser.add_argument('--assembly',
                                 help='Name of the genome assembly, e.g., hg19, hg38, mm10, default: %(default)s',
                                 choices=['hg19', 'hg38', 'mm10'], default='mm10')
    pipeline_parser.set_defaults(func=pipeline.run)
    
    description = 'Initialize RNA-Seq pipeline'
    init_parser = subparsers.add_parser('init', help=description, description=description)
    init_parser.set_defaults(func=init.init)
    
    description = 'Process FASTQ files using fastp program'
    fastp_parser = subparsers.add_parser('fastp', help=description, description=description)
    fastp_parser.add_argument('fastq', help='Forward read (read 1) FASTQ file(s)', nargs='+', type=cmder.File)
    fastp_parser.set_defaults(func=fastp.main)
    
    description = 'Map FASTQ files using STAR program'
    map_parser = subparsers.add_parser('map', help=description, description=description)
    map_parser.add_argument('--genome', help='Path to a directory contains STAR index files', type=cmder.Dir)
    map_parser.add_argument('fastq', help='Forward read (read 1) FASTQ file(s)', nargs='+', type=cmder.File)
    
    description = 'Perform feature count using featureCounts program'
    count_parser = subparsers.add_parser('count', help=description, description=description)
    count_parser.add_argument('--gtf', help='Path to a GTF file contains genomic annotation', type=cmder.File)
    count_parser.add_argument('--bed', help='Path to a BED file contains genomic annotation', type=cmder.File)
    count_parser.add_argument('bam', help='Path to BAM file(s)', nargs='+', type=cmder.Files)
    
    description = 'Determine gene differential expression'
    de_parser = subparsers.add_parser('dega', help=description, description=description)
    de_parser.add_argument('table', help='Path to a TSV file contains feature count table', type=cmder.File)
    de_parser.add_argument('--assembly',
                           help='Name of the genome assembly, e.g., hg19, hg38, mm10, default: %(default)s',
                           choices=['hg19', 'hg38', 'mm10'], default='mm10')
    de_parser.set_defaults(func=dega.dega)
    
    args = parser.parse_args()
    print(args)
    args.func(args)


if __name__ == '__main__':
    main()

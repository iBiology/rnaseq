#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Run the pipeline
"""

from pathlib import Path
from cmder import File, Dir
from rnaseq import utility, fastp, star, count, tools


class RNASeq(tools.Pipeline):
    manifest: File = None  # Path to a manifest file
    genome: Dir = None  # Path to STAR the directory contains STAR genome index files
    gtf: File = None  # Path to a GTF file contains genomic annotations
    bed: File = None  # Path to a BED file contains genomic annotations
    
    def pre_process(self):
        self.outdir = self.mkdir(self.outdir or Path.cwd())
        
    def command(self):
        cmd = (f'rnaseq --manifest {self.manifest.path} --genome {self.genome.path} '
               f'--gtf {self.gtf.path} --bed {self.bed.path} --out {self.outdir} '
               f'--cpu {self.cpu}{self.qvd}')
        return cmd
    
    def process(self):
        df, bams = tools.load_manifest(self.manifest.path), []
        for row in df.itertuples():
            sample_name = row.SampleName
            outdir = self.outdir / 'fastp'
            f1, f2 = fastp.fastq(row.FASTQ1, row.FASTQ2, row.SampleName, outdir, dryrun=self.dryrun)
            
            outdir = self.outdir / 'bam'
            bam = star.star(f1, f2, sample_name, outdir, self.genome.path, self.cpu, dryrun=self.dryrun)
            bams.append(bam)
            
        


def run(args):
    outdir = utility.mkdir(args)
    
    setattr(args, 'outdir', outdir / 'fastp')
    files = fastp.fastp(args)
    
    setattr(args, 'files', files)
    setattr(args, 'outdir', outdir / 'bam')
    bam = star.star(args)
    
    setattr(args, 'bam', bam)
    setattr(args, 'outdir', outdir / 'count')
    table = count.count(args)
    
    setattr(args, 'table', table)
    setattr(args, 'outdir', outdir)

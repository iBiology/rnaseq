#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Run the pipeline
"""

from pathlib import Path
from cmder import File, Dir
from rnaseq import fastp, star, count, tools


class RNASeq(tools.Pipeline):
    manifest: File  # Path to a manifest file
    genome: Dir  # Path to STAR the directory contains STAR genome index files
    gtf: File  # Path to a GTF file contains genomic annotations
    bed: File  # Path to a BED file contains genomic annotations
    
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
            
            
def main():
    RNASeq().parse_args().fire()
    
    
if __name__ == '__main__':
    main()

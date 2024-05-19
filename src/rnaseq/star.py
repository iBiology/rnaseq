#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Perform feature counts using featureCount program
"""

from pathlib import Path

import cmder
from cmder import File, Dir
from rnaseq import tools


def star(f1, f2, sample_name, outdir, genome, cpu, dryrun=False):
    bam, log, wd = outdir / f'{sample_name}.bam', outdir / f'{sample_name}.star.log', outdir / sample_name
    if bam.exists():
        cmder.logger.info(f'BAM file {bam} exists, skip re-mapping')
    else:
        if dryrun:
            cmder.logger.info(f'Will create directory: {wd}')
        else:
            wd.mkdir(exist_ok=True)
    
        cmd = (f'STAR \\\n  --genomeDir {genome} \\\n  --runThreadN {cpu} \\\n  --outFileNamePrefix {wd}/'
               f' \\\n  --outSAMtype BAM Unsorted \\\n  --outSAMunmapped None')
        reads = f'{f1} \\\n                {f2}' if f2.exists() else f1
        cmd = f"{cmd} \\\n  --readFilesCommand 'zcat <' \\\n  --readFilesIn {reads}"
    
        cmder.logger.info(f'Mapping reads of sample {sample_name} using STAR')
        if dryrun:
            cmder.logger.info(cmd)
        else:
            p = cmder.run(cmd, fmt_cmd=False)
            if p.returncode:
                cmder.logger.error(f'STAR failed with return code {p.returncode}')
            else:
                cmder.run(f'mv {wd}/Log.final.out {log}')
                memory = int(32 / cpu) or 1
                cmder.run(f'samtools sort -@ {cpu} -m {memory}G -o {bam} {wd}/Aligned.out.bam',
                          exit_on_error=True)
                cmder.run(f'samtools index -@ {cpu} {bam}', exit_on_error=True)
                cmder.run(f'rm -r {wd}')
    return bam


class STAR(tools.Pipeline):
    fastq1: File  # Forward read (read 1) FASTQ file(s)
    fastq2: File = None  # Reverse read (read 2) FASTQ file(s)
    sample_name: str = ''  # The corresponding sample name of each FASTQ file (or pair)
    genome: Dir  # Path to STAR the directory contains STAR genome index files
    
    def pre_process(self):
        name = self.sample_name or tools.fastq_basename(self.fastq1.path)
        bam = self.outdir / f'{name}.bam'
        if bam.exists():
            self.logger.info(f'BAM file {bam} exists, skip re-mapping')
    
    def process(self):
        f2 = self.fastq2.path if self.fastq2 else None
        name = self.sample_name or tools.fastq_basename(self.fastq1.path)
        star(self.fastq1.path, f2, name, self.outdir, self.genome.path, self.cpu, dryrun=self.dryrun)
    
    def command(self):
        sample_name = f' --sample_name {self.sample_name or tools.fastq_basename(self.fastq1.path)}'
        f2 = f' --fastq2 {self.fastq2.path}' if self.fastq2 else ''
        cmd = (f'rnaseq-star --fastq1 {self.fastq1.path}{f2}{sample_name} --outdir {self.outdir} '
               f'--genome {self.genome.path} --cpu {self.cpu}{self.qvd}')
        return cmd


def main():
    STAR().parse_args().fire()
    
    
if __name__ == '__main__':
    main()

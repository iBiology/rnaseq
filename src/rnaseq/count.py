#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Perform feature count using featurCounts program
"""

import sys
import multiprocessing
from pathlib import Path

import cmder
from cmder import File
import pandas as pd
from rnaseq import utility, tools


def find_strand(bam, bed, dry):
    cmd = f'infer_experiment.py -r {bed} -i {bam}'
    data, strand, x, y = 'paired', '0', 0, 0
    
    if dry:
        utility.logger.info(f'{cmd} [DRYRUN (with fake strandness)]')
    else:
        p = cmder.run(cmd)
        for line in p.stdout.readlines():
            utility.logger.info(line.strip())
            if 'SingleEnd Data' in line:
                data = 'single'
            if '1++,1--,2+-,2-+' in line or '++,--' in line:
                x = float(line.strip().split(': ')[1])
            elif '1+-,1-+,2++,2--' in line or '+-,-+' in line:
                y = float(line.strip().split(': ')[1])
        if y:
            if x / y > 1.25:
                strand = '1'
            elif x / y < 0.75:
                strand = '2'
            else:
                strand = '0'
        else:
            strand = '1'
    return bam, (data, strand)


def worker(item):
    return find_strand(*item)


def counting(bams, gtf, bed, outdir, cpu, dryrun=False):
    processes = min(len(bams), cpu)
    with multiprocessing.Pool(processes=processes) as pool:
        meta = dict(pool.map(worker, [(bam, bed, dryrun) for bam in bams]))
    
    data = set(v[0] for v in meta.values())
    if len(data) == 1:
        data = list(data)[0]
    else:
        cmder.logger.error(f'BAM files are mixed with both paired-end and single-end data, cannot proceed')
    
    out, output, log = outdir / 'count', outdir / 'count.tsv', outdir / 'count.log'
    meta = {k: v[1] for k, v in meta.items()}
    exe = f'featureCounts -p' if data == 'paired' else 'featureCounts'
    files, strands = [], []
    for k, v in meta.items():
        files.append(f'{k}')
        strands.append(v)
    
    files, strand = ' \\\n  '.join(files), ','.join(strands)
    cmd = f'{exe} \\\n  -a {gtf} \\\n  -s {strand} \\\n  -o {out} \\\n  {files} \\\n  &> {log}'
    _ = utility.logger.info(cmd) if dryrun else cmder.run(cmd, fmt_cmd=False, exit_on_error=True)
    
    df = pd.read_csv(out, sep='\t', skiprows=1)
    df = df.drop(columns=['Chr', 'Start', 'End', 'Strand', 'Length'])
    df.columns = [Path(c).name.removesuffix('.bam') for c in df.columns]
    df.to_csv(output, index=False, sep='\t')
    cmder.logger.info(f'Feature count saved to {output}')


class Count(tools.Pipeline):
    bam: list[File, ]  # One or multiple BAM files need to count features
    gtf: File  # Path to a GTF file contains genomic annotations
    bed: File  # Path to a BED file contains genomic annotations
    
    def pre_process(self):
        output = self.outdir / 'feature.count.tsv'
        if output.exists():
            self.logger.info(f'Output file already exists: {output}', terminate=True)
            
    def command(self):
        bam = ' '.join([bam.path for bam in self.bam])
        return (f'rnaseq-count -bam {bam} -gtf {self.gtf.path} --bed {self.bed.path} --outdir {self.outdir} '
                f'--cpu {self.cpu}{self.qvd}')
    
    def process(self):
        counting([bam.path for bam in self.bam], self.gtf.path, self.bed.path, self.outdir,
                 self.cpu, dryrun=self.dryrun)


if __name__ == '__main__':
    Count().parse_args().fire()

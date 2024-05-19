#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Process FASTQ files using fastp program
"""

from pathlib import Path

import cmder
from cmder import File
from rnaseq import tools


def fastq(f1, f2, sample_name, outdir, dryrun=False):
    o1, o2 = outdir / f'{sample_name}.r1.fastq.gz', outdir / f'{sample_name}.r2.fastq.gz'
    data, html = outdir / f'{sample_name}.fastp.json', outdir / f'{sample_name}.fastp.html'

    cmd = f'fastp -5 -3 -W 4 -M 20 -l 15 -x -n 5 -z 9 -w 8 \\\n  -j {data} \\\n  -h {html} '
    if o1.exists():
        if f2:
            if o2.exists():
                cmd = ''
                cmder.logger.info(f'Output files for paired-end sample {sample_name} already exists')
            else:
                cmder.logger.info(f'Processing paired-end FASTQ files for sample {sample_name} (due to missing {o2})')
                cmd = f'{cmd} \\\n  -i {f1} \\\n  -o {o1} \\\n  -I {f2} \\\n  -O {o2}'
        else:
            cmd = ''
            cmder.logger.info(f'Output file for single-end sample {sample_name} already exists')
    else:
        if f2:
            cmder.logger.info(f'Processing paired-end FASTQ files for sample {sample_name}')
            cmd = f'{cmd} \\\n  -i {f1} \\\n  -o {o1} \\\n  -I {f2} \\\n  -O {o2}'
        else:
            cmder.logger.info(f'Processing single-end FASTQ files for sample {sample_name}')
            cmd = f'{cmd} \\\n  -i {f1} \\\n  -o {o1}'
    
    if cmd:
        _ = cmder.logger.info(cmd) if dryrun else cmder.run(cmd, fmt_cmd=False)
    return o1, o2


class FASTP(tools.Pipeline):
    fastq1: File  # Forward read (read 1) FASTQ file(s)
    fastq2: File = None  # Reverse read (read 2) FASTQ file(s)
    sample_name: str = None  # The corresponding sample name of each FASTQ file (or pair)
    
    def pre_process(self):
        o1, o2 = self.outdir / f'{self.sample_name}.r1.fastq.gz', self.outdir / f'{self.sample_name}.r2.fastq.gz'
        if o1.exists():
            if self.fastq2:
                if o2.exists():
                    self.logger.info(f'Output files for paired-end sample {self.sample_name} already exist',
                                     terminate=True)
                else:
                    self.logger.info(f'Output file for single-end sample {self.sample_name} already exist',
                                     terminate=True)
    
    def process(self):
        sample_name = self.sample_name or tools.fastq_basename(self.fastq1.path)
        f2 = self.fastq2.path if self.fastq2 else None
        fastq(self.fastq1.path, f2, sample_name, outdir=self.outdir, dryrun=self.dryrun)
    
    def command(self):
        sample_name = f' --sample_name {self.sample_name or tools.fastq_basename(self.fastq1.path)}'
        f2 = f' --fastq2 {self.fastq2.path}' if self.fastq2 else ''
        cmd = (f'rnaseq-fastp --fastq1 {self.fastq1.path}{f2}{sample_name} --outdir {self.outdir} '
               f'--cpu {self.cpu}{self.qvd}')
        return cmd
    
    
def main():
    FASTP().parse_args().fire()


if __name__ == '__main__':
    main()

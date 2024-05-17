#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Perform feature counts using featureCount program
"""

import cmder
from rnaseq import utility


def star(args):
    files, outdir = [], utility.mkdir(args)
    for f1 in args.files:
        name = f1.name.removesuffix('.r1.fastq.gz')
        f2 = f1.parent / f1.name.replace('.r1.fastq.gz', '.r2.fastq.gz')
        bam, log, wd = outdir / f'{name}.bam', outdir / f'{name}.log', outdir / name
        
        if args.dry:
            utility.logger.info(f'Will create directory: {wd}')
        else:
            wd.mkdir(exist_ok=True)
        
        cmd = (f'STAR \\\n  --genomeDir {args.genome} \\\n  --runThreadN {args.process} \\\n  --outFileNamePrefix {wd}/'
               f' \\\n  --outSAMtype BAM Unsorted \\\n  --outSAMunmapped None')
        if bam.exists():
            utility.logger.info(f'BAM file {bam} exists, skip re-mapping')
            cmd = ''
        else:
            reads = f'{f1} \\\n                {f2}' if f2.exists() else f1
            cmd = f"{cmd} \\\n  --readFilesCommand 'zcat <' \\\n  --readFilesIn {reads}"

        if cmd:
            utility.logger.info(f'Mapping reads of sample {name} using STAR')
            if args.dry:
                utility.logger.info(cmd)
            else:
                p = cmder.run(cmd, fmt_cmd=False)
                if p.returncode:
                    utility.logger.error(f'STAR failed with return code {p.returncode}')
                else:
                    cmder.run(f'mv {wd}/Log.final.out {log}')
                    memory = int(32 / args.process) or 1
                    cmder.run(f'samtools sort -@ {args.process} -m {memory}G -o {bam} {wd}/Aligned.out.bam',
                              exit_on_error=True)
                    cmder.run(f'samtools index -@ {args.process} {bam}', exit_on_error=True)
                    cmder.run(f'rm -r {wd}')
        
        files.append(bam)
    return files
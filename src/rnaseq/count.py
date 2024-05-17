#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Perform feature count using featurCounts program
"""

import sys
import multiprocessing
from pathlib import Path

import cmder
import pandas as pd
from rnaseq import utility


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


def count(args):
    outdir = utility.mkdir(args)
    out, output, log = outdir / 'count', outdir / 'count.tsv', outdir / 'count.log'
    if output.exists():
        utility.logger.info(f'Feature count output file already exists: {output}')
    else:
        processes = min(len(args.bam), args.process)
        with multiprocessing.Pool(processes=processes) as pool:
            meta = dict(pool.map(worker, [(bam, args.bed, args.dry) for bam in args.bam]))
            
        data = set(v[0] for v in meta.values())
        if len(data) == 1:
            data = list(data)[0]
        else:
            utility.logger.error(f'BAM files are mixed with both paired-end and single-end data, cannot proceed')
            sys.exit(1)
        
        meta = {k: v[1] for k, v in meta.items()}
        exe = f'featureCounts -p' if data == 'paired' else 'featureCounts'
        files, strands = [], []
        for k, v in meta.items():
            files.append(f'{k}')
            strands.append(v)
    
        files, strand = ' \\\n  '.join(files), ','.join(strands)
        cmd = f'{exe} \\\n  -a {args.gtf} \\\n  -s {strand} \\\n  -o {out} \\\n  {files} \\\n  &> {log}'
        _ = utility.logger.info(cmd) if args.dry else cmder.run(cmd, fmt_cmd=False)
    
        if out.exists():
            df = pd.read_csv(out, sep='\t', skiprows=1)
            df = df.drop(columns=['Chr', 'Start', 'End', 'Strand', 'Length'])
            df.columns = [Path(c).name.removesuffix('.bam') for c in df.columns]
            df.to_csv(output, index=False, sep='\t')
            utility.logger.info(f'Feature count output file created: {output}')

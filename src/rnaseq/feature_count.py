#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Count features using BAM files and parse the results.
"""

import os
import sys
import argparse

import pandas as pd
from seqflow import logger, task, Flow
import cmder

import utility

parser = argparse.ArgumentParser(description=__doc__, prog='feature_count')
parser.add_argument('bam', nargs='+', help=f'Path to BAM file(s)')
parser.add_argument('-i', '--ids', help=f'Unique identifier for each BAM file.', nargs='*')
parser.add_argument('-g', '--gtf', help=f'Path to annotation GTF file.', required=True)
parser.add_argument('-n', '--cpus', help=f'Number of CPUs can be used, default: %(default)s', default=16, type=int)
parser.add_argument('-c', '--counter', choices=('featureCounts', 'htseq-count'),
                    help=f'Name of the feature count program, default: %(default)s', default='featureCounts')
parser.add_argument('-o', '--outdir', help="Path to the output directory, default: current work directory.")
parser.add_argument('--dryrun', action='store_true',
                    help='Print out steps and files involved in each step without actually processing data.')
parser.add_argument('--debug', action='store_true',
                    help='Invoke debug mode and print out more logs.')

args = parser.parse_args()
if os.path.isfile(args.gtf):
    GTF = args.gtf
else:
    logger.error(f'GTF {args.gtf} may not be a file does not exist')
    sys.exit()

setattr(args, 'outdir', args.outdir or os.getcwd())
try:
    os.makedirs(args.outdir, exist_ok=True)
except OSError as e:
    logger.error(f'Create outdir failed: {e}.')
    sys.exit(1)
# os.chdir(args.outdir)

BAM2ID, BAM2COUNT, COUNT2ID = {}, {}, {}
ids = args.ids or [os.path.basename(bam)[:-4] for bam in args.bam]
if len(ids) != len(args.bam):
    logger.error('The number of BAM files and unique identifies does equal')
    sys.exit(1)
else:
    BAM2ID = {bam: uid for bam, uid in zip(args.bam, ids)}
    name = "feature" if args.counter == "featureCounts" else "htseq"
    for bam, uid in zip(args.bam, ids):
        count = os.path.join(args.outdir, f'{uid.split("___")[0]}.{name}.count')
        BAM2COUNT[bam] = count
        COUNT2ID[count] = uid


def feature_count(counter, bam, gtf, out):
    read_type, strand = utility.guess_strand(bam, gtf.replace('.gtf', '.bed'))
    if counter == 'featureCounts':
        # exe = 'featureCounts -pBP' if read_type == 'paired' else 'featureCounts'
        exe = 'featureCounts -p' if read_type == 'paired' else 'featureCounts'
        cmd = f'{exe} -t exon -g gene_id -a {gtf} -s {strand} -o {out} {bam} 2> {out}.log'
    elif counter == 'htseq-count':
        strand = {'0': 'no', '1': 'yes', '2': 'reverse'}[strand]
        exe = 'htseq-count -r pos' if read_type == 'paired' else 'htseq-count'
        cmd = f'{exe} -f bam -s {strand} -q {bam} {gtf} > {out}'
    else:
        logger.error(f'Feature count using {counter} not support yet')
        sys.exit(1)
    cmder.run(cmd, msg='Counting features using featureCounts ...', fmt_cmd=False)


def parse_count(out):
    names = ['feature_id', COUNT2ID[out]]
    if out.endswith('htseq.count'):
        df = pd.read_csv(out, sep='\t', header=None, dtype=str, names=names)
    else:
        df = pd.read_csv(out, sep='\t', skiprows=1)
        df = df.drop(columns=['Chr', 'Start', 'End', 'Strand', 'Length'])
        df.columns = names
    return df


def annotate_count(df, out):
    dd = pd.read_csv(out, sep='\t', dtype=str).rename(columns={'Unnamed: 0': 'feature_id'})
    if out.endswith('.tsv'):
        out = f'{out[:-4]}.annotate.tsv'
    else:
        dd.to_csv(f'{out}.tsv', sep='\t', index=False)
        cmder.run(f'rm {out}')
        out = f'{out}.annotate.tsv'
    dd = pd.merge(dd, df, on='feature_id', how='left')
    dd.to_csv(out, sep='\t', index=False)


@task(inputs=args.bam, cpus=args.cpus, outputs=lambda i: BAM2COUNT[i])
def count(bam, txt):
    feature_count(args.counter, bam, GTF, txt)


@task(inputs=[], parent=count, outputs=[os.path.join(args.outdir, 'feature.count.tsv')])
def merge_count(inputs, output):
    outputs = [BAM2COUNT[bam] for bam in args.bam]
    df, ss = pd.DataFrame(), pd.DataFrame()
    for out in outputs:
        dd = parse_count(out)
        df = dd if df.empty else pd.merge(df, dd, on='feature_id')
        summary = f'{out}.summary'
        if os.path.exists(summary):
            ds = pd.read_csv(summary, sep='\t')
            ds.columns = ['Status', os.path.basename(out).replace('.feature.count', '')]
            ss = ds if ss.empty else pd.merge(ss, ds, on='Status')
    out = os.path.join(args.outdir, 'feature.count.tsv')
    df.to_csv(out, sep='\t', index=False)
    if not ss.empty:
        ss.to_csv(os.path.join(args.outdir, 'feature.count.summary'), sep='\t', index=False)

    cpm, tpm = os.path.join(args.outdir, 'feature.cpm'), os.path.join(args.outdir, 'feature.tpm')
    cmd = f'rnanorm {out} --cpm-output={cpm} --tpm-output={tpm} --annotation={GTF}'
    cmder.run(cmd, msg="Normalizing feature counts using rnanorm ...")
    
    anno = utility.parse_gtf(GTF)
    annotate_count(anno, cpm)
    annotate_count(anno, tpm)
    annotate_count(anno, output)


@logger.catch()
def main():
    flow = Flow('feature-count', description=__doc__.strip())
    flow.run(dry_run=args.dryrun, cpus=args.cpus)


if __name__ == '__main__':
    main()


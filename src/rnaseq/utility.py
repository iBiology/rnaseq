#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from pathlib import Path

import cmder
import pandas as pd
from loguru import logger

VERSION = 1.0
HEADER = fr"""
             _____   _   _
            |  __ \ | \ | |    /\
            | |__) ||  \| |   /  \    ___   ___   __ _
            |  _  / | . ` |  / /\ \  / __| / _ \ / _` |
            | | \ \ | |\  | / ____ \ \__ \|  __/| (_| |
            |_|  \_\|_| \_|/_/    \_\|___/ \___| \__, |
                                                    | |
                           VERSION {VERSION}              |_|
"""
FOOTER = """
+======================================================================================================================+
|                                                                                                                      |
|                                                 MISSION ACCOMPLISHED                                                 |
|{hh_mm_ss}|
|                                                                                                                      |
+======================================================================================================================+
"""


logger.remove()
logger.add(sys.stdout, format="[{time:HH:mm:ss}] <level>{message}</level>", colorize=True, level='INFO')


def mkdir(args):
    outdir = args.outdir.resolve()
    if not outdir.exists():
        if args.dry:
            logger.info(f'Will creat {outdir}')
        else:
            logger.info(f'Creating {outdir}')
            outdir.mkdir(exist_ok=True, parents=True)
            logger.info(f'Directory created: {outdir}')
    return outdir


def guess_strand(bam, bed):
    p = cmder.run(f'infer_experiment.py -r {bed} -i {bam}')
    data, strand, x, y = 'paired', '', 0, 0
    for line in p.stdout.readlines():
        cmder.logger.info(line.strip())
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
    return data, strand


def error_and_exit(message):
    logger.error(message)
    sys.exit(1)


class File:
    def __init__(self, path):
        try:
            if not Path(path).is_file():
                error_and_exit(f'File {path} is not a file or does not exist')
        except Exception as e:
            error_and_exit(f'Check file {path} failed due to {e}')
        self.path = Path(path).resolve()
    
    def __call__(self, *args, **kwargs):
        return self.path
    
    def __getattr__(self, item):
        return getattr(self.path, item)


class Files:
    def __init__(self, paths):
        self.paths = [File(path) for path in paths]
    
    def __call__(self, *args, **kwargs):
        return self.paths


class Dir:
    def __init__(self, path):
        try:
            if not Path(path).is_dir():
                error_and_exit(f'Directory {path} is not a directory or does not exist')
        except Exception as e:
            error_and_exit(f'Check directory {path} failed due to {e}')
        self.path = Path(path).resolve()
    
    def __call__(self, *args, **kwargs):
        return self.path
    
    def __getattr__(self, item):
        return getattr(self.path, item)


def load_manifest(args):
    df = pd.read_csv(args.manifest, comment='#')
    if 'SampleName' not in df.columns:
        logger.error(f'No SampleName column found in manifest file')
        sys.exit(1)
    if 'FASTQ1' not in df.columns:
        logger.error(f'No FASTQ1 column found in manifest file')
        sys.exit(1)
        
    df['FASTQ1'] = df['FASTQ1'].apply(lambda x: File(x, root=args.manifest.parent))
    if 'FASTQ2' in df.columns:
        df['FASTQ2'] = df['FASTQ2'].apply(lambda x: File(x, root=args.manifest.parent))
    else:
        df['FASTQ2'] = None
    return df

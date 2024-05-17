#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Run the pipeline
"""

from rnaseq import utility, fastp, star, count


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

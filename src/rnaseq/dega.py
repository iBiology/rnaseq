#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Different Expression Analysis with pyDEseq2
"""

from itertools import combinations

import cmder
import pandas as pd
from rnaseq import utility


def dega(args):
    df = pd.read_csv(args.table, sep="\t", index_col=0)
    df = df[df.sum(axis=1) >= 10]
    print(df)
    
    groups = set(c.rsplit("_", 1)[0] for c in df.columns)
    comparisons = {c: c.split('.vs.') for c in args.comparison}
    for comparison, (group1, group2) in comparisons.items():
        if group1 not in groups:
            utility.logger.error(f'Invalid comparison {comparison}, group {group1} not in {groups}')
            continue
        elif group2 not in groups:
            utility.logger.error(f'Invalid comparison {comparison}, group {group2} not in {groups}')
            continue
        utility.logger.info(f'Performing differential expression analysis for {comparison}')
    

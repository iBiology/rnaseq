#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Process FASTQ files using fastp program
"""

import cmder
from rnaseq import utility


def fastp(args):
    files, outdir, df = [], utility.mkdir(args), utility.load_manifest(args)
    for row in df.itertuples():
        name, f1, f2 = row.SampleName, row.FASTQ1, row.FASTQ2
        data, html = outdir / f'{name}.fastp.json', outdir / f'{name}.fastp.html'
        o1, o2 = outdir / f'{name}.r1.fastq.gz', outdir / f'{name}.r2.fastq.gz'
        cmd = f'fastp -5 -3 -W 4 -M 20 -l 15 -x -n 5 -z 9 -w 8 \\\n  -j {data} \\\n  -h {html} '
        if o1.exists():
            if f2.exists():
                if o2.exists():
                    cmd = ''
                    utility.logger.info(f'Output files for paired-end sample {name} already exists')
                else:
                    utility.logger.info(f'Processing paired-end FASTQ files for sample {name} (due to missing {o2})')
                    cmd = f'{cmd} \\\n  -i {f1} \\\n  -o {o1} \\\n  -I {f2} \\\n  -O {o2}'
            else:
                cmd = ''
                utility.logger.info(f'Output file for single-end sample {name} already exists')
        else:
            if f2.exists():
                utility.logger.info(f'Processing paired-end FASTQ files for sample {name}')
                cmd = f'{cmd} \\\n  -i {f1} \\\n  -o {o1} \\\n  -I {f2} \\\n  -O {o2}'
            else:
                utility.logger.info(f'Processing single-end FASTQ files for sample {name}')
                cmd = f'{cmd} \\\n  -i {f1} \\\n  -o {o1}'
        
        if cmd:
            _ = utility.logger.info(cmd) if args.dry else cmder.run(cmd, fmt_cmd=False)
            
        files.append(o1)
    return files

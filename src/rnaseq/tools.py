#!/usr/bin/env python
# -*- coding: utf-8 -*-

# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Useful functions and Classes as utilities for building virtual screening pipeline
"""

import os
import sys
from pathlib import Path

import cmder
import pandas as pd
from tap import Tap
from loguru import logger


class Pipeline(Tap):
    nodes: int = 0  # Number of nodes could be requested
    hold: bool = False  # Flag to hold the submission
    cpu: int = 8  # Number of CPUs could be used
    outdir: Path = None  # Path to a directory for saving output(s)
    
    job_name: str = 'RNA-Seq'  # Meaningful name for the job
    email: str = ''  # Email address for notify user the status of the job
    dryrun: bool = False  # Flag to enable dryrun mode to only messaging out processing steps
    quiet: bool = False  # Flag to enable quiet mode to suppress verbose messages
    verbose: bool = False  # Flag to enable verbose mode to emit verbose messages
    develop: bool = False  # Flag to enable develop mode (for development use only)
    
    # @property
    # def name(self):
    #     return self.__class__.__name__
    #
    # @property
    # def output(self):
    #     return []
    #
    # def check_output(self):
    #     for output in self.output:
    #         if output:
    #             if Path(output).exists():
    #                 self.logger.debug(f"Output file {output} already exists")
    #             else:
    #                 self.logger.error(f"Output file {output} does not exist")
    #
    # def check_fastq_manifest(self):
    #     manifest, fastq1, fastq2, sample_name = [getattr(self, x, None)
    #                                              for x in ('manifest', 'fastq1', 'fastq2', 'sample_name')]
    #     if manifest:
    #         pass
    #     elif fastq1:
    #         if sample_name:
    #             if fastq2:
    #                 if len(fastq1) == len(fastq2) == len(sample_name):
    #                     pass
    #                 else:
    #                     self.error(f'The number of FASTQ files(s) does not match the number of sample name(s)')
    #             else:
    #                 if len(fastq1) == len(sample_name):
    #                     fastq2 = [None] * len(fastq1)
    #                 else:
    #                     self.error(f'The number of forward FASTQ files(s) does not match the number of sample name(s)s')
    #         else:
    #             self.error('FASTQ1 file(s) was provided, but no sample name was provided')
    #     else:
    #         self.error('Neither manifest nor FASTQ file(s) was provided')
    
    def pre_process(self):
        pass
    
    def process_args(self) -> None:
        self.outdir = self.mkdir(Path(self.outdir or Path.cwd()))
        # self.check_output()
        self.pre_process()
    
    @property
    def logger(self):
        logger.remove()
        level = "DEBUG" if self.verbose or self.develop else ("ERROR" if self.quiet else "INFO")
        formatter = '<level>[{time:YYYY-MM-DD HH:mm:ss}] {message}</level>'
        logger.add(sys.stdout, colorize=True, format=formatter, level=level)
        return logger
    
    def debug(self, message, terminate=False):
        self.logger.debug(message)
        if terminate:
            sys.exit(0)
    
    def info(self, message, terminate=False):
        self.logger.info(message)
        if terminate:
            sys.exit(0)
    
    def error(self, message):
        self.logger.error(message)
        sys.exit(1)
    
    @property
    def qvd(self):
        s = ' '.join([f'--{k}' for k in ('quiet', 'verbose', 'develop') if getattr(self, k)])
        if self.task:
            s = f'{s} --task {self.task} --status {self.status}'
        return f' {s}' if s else ''
    
    @property
    def command(self):
        raise NotImplementedError
    
    @property
    def cwd(self):
        return self.outdir
    
    def process(self):
        raise NotImplementedError
    
    def mkdir(self, p):
        try:
            p = Path(p).resolve()
            p.mkdir(parents=True, exist_ok=True)
        except IOError as e:
            self.error(f'Failed to make directory {p} due to {e}')
        return p
    
    def run(self):
        try:
            self.info(f'{self.name} ...')
            self.process()
            self.info(f'{self.name} complete.')
        except SystemExit as e:
            if e.code:
                self.error(f'Exit early due to {e}')
    
    def submit(self, **kwargs):
        nodes = 1
        if self.nodes > 1:
            self.debug(f'The current version only support 1 node, set nodes=1')
            
        script = self.outdir / 'submit.sh'
        with script.open('w') as o:
            o.write('#!/usr/bin/env bash\n')
            o.write('\n')
            o.write('#BSUB -W 48:00\n')
            o.write('#BSUB -q long\n')
            o.write(f'#BSUB -cwd {self.outdir}\n')
            o.write(f'#BSUB -nnodes {nodes}\n')
            o.write('#BSUB -n 28\n')
            o.write('#BSUB -M 64GB\n')
            o.write('#BSUB -R rusage[mem=2GB]\n')
            o.write(f'#BSUB -J {self.job_name}\n')
            if self.email:
                o.write(f'#BSUB -u {self.email}\n')
            o.write('\n')
            o.write(f'source {Path(sys.executable).parent}/activate\n')
            o.write(f'wd={self.outdir}\n')
            o.write('cd "${wd} || { echo "Failed cd into ${wd}"; exit 1; }\n\n')
            o.write(f'{self.command}\n')
            o.write('\n')
        logger.info(f'Submission script was saved to {script}')
        
        if self.hold or self.dryrun:
            self.info(f'Job was not submitted yet due to hold or dryrun flag was set to True ...')
        else:
            cmder.run(f'bsub < {script}')
    
    def fire(self):
        _ = self.submit() if self.nodes or self.dryrun else self.run()
        
        
def fastq_basename(fn: Path) -> str:
    known_suffixes, fn = ['.fastq.gz', '.fq.gz', '.fastq', '.fq'], fn.name
    for suffix in known_suffixes:
        fn = fn.removesuffix(suffix)
    return fn


def load_manifest(manifest):
    df = pd.read_csv(manifest, comment='#')
    if 'SampleName' not in df.columns:
        logger.error(f'No SampleName column found in manifest file')
        sys.exit(1)
    if 'FASTQ1' not in df.columns:
        logger.error(f'No FASTQ1 column found in manifest file')
        sys.exit(1)
    
    df['FASTQ1'] = df['FASTQ1'].apply(lambda x: cmder.File(x).path)
    if 'FASTQ2' in df.columns:
        df['FASTQ2'] = df['FASTQ2'].apply(lambda x: cmder.File(x).path)
    else:
        df['FASTQ2'] = None
    return df


if __name__ == '__main__':
    pass

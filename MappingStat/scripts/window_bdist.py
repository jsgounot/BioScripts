# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2025-03-11 15:25:49
# @Last modified by:   jsgounot
# @Last Modified time: 2025-03-11 15:28:20

import sys, subprocess, io
import pandas as pd

def process_chunk(sdf):
    sdf = sdf.groupby('contig')['depth'].value_counts().rename('count')
    return sdf.reset_index()

def make_depth_count(bamfile):
    cmdline = f'samtools depth -aa {bamfile}'
    sp = subprocess.Popen(cmdline.split(), stdout=subprocess.PIPE)
    sio = io.TextIOWrapper(sp.stdout)
    names = ['contig', 'position', 'depth']
    iterator = pd.read_csv(sio, sep='\t', names=names, iterator=True, chunksize=100000)
    df = pd.concat(process_chunk(sdf) for sdf in iterator)
    df = df.groupby(['contig', 'depth'])['count'].sum()
    return df.reset_index()

bamfile = snakemake.input[0]
df = make_depth_count(bamfile)
outfile = snakemake.output[0]
df.to_csv(outfile, sep='\t', index=False)
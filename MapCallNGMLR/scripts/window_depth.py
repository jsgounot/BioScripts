# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2020-12-09 09:56:48
# @Last Modified by:   jsgounot
# @Last Modified time: 2020-12-09 10:13:45


import pandas as pd

WINDOWSIZE = 10000
QUANTILE = .95 # < 1
HEADER = ('contig', 'position', 'coverage')

depthfile = snakemake.input[0]
outfile   = snakemake.output[0]

if depthfile.endswith(".gz") :
	df = pd.read_csv(depthfile, compression='gzip', sep="\t", names=HEADER)
else :
	df = pd.read_csv(depthfile, sep="\t", names=HEADER)

df["window"] = (df["position"] // WINDOWSIZE) * WINDOWSIZE
df = df.groupby(["contig", "window"])["coverage"].mean()
df = df.rename("mean coverage").reset_index()

if QUANTILE :
	quantile = df["mean coverage"].quantile(QUANTILE)
	df.loc[df["mean coverage"] > quantile, "mean coverage"] = quantile

df.to_csv(outfile, sep="\t")
# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-05-31 15:42:07
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-05-31 15:53:02

import os
import pandas as pd

fname = snakemake.input[0]
check = snakemake.output[0]

df = pd.read_csv(fname, sep="\t", index_col=0)
df = df[df["eigen_index"] == 0]

outdir = os.path.dirname(check)
for path, cid in zip(df["sequence"], df["clusterID"]) :
	outfname = os.path.join(outdir, "cluster" + str(cid) + ".fasta")
	os.symlink(path, outfname)

with open(check, "w") : 
	pass
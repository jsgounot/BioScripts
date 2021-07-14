# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-07-12 15:17:46
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-07-12 15:38:38

import pandas as pd

fname = snakemake.input[0]

names = ["query-ID", "reference-ID", "distance", "p-value", "shared-hashes"]
df = pd.read_csv(fname, sep="\t", names=names, compression="gzip")

df = df[df.groupby("query-ID")["distance"].transform(min) == df["distance"]]

outfile = snakemake.output[0]
df.to_csv(outfile, sep="\t")
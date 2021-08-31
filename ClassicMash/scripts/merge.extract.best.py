# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-07-12 15:17:52
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-07-12 15:34:55

import pandas as pd

fnames = sorted(snakemake.input)

df = pd.concat(
	(pd.read_csv(fname, sep="\t")
		for fname in fnames))

df = df[df.groupby("query-ID")["distance"].transform(min) == df["distance"]]

outfile = snakemake.output[0]
df.to_csv(outfile, sep="\t", index=False)
# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-08-30 17:40:50
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-08-31 09:13:12

import pandas as pd

depth_in  = snakemake.input[0]
idxstats  = snakemake.input[1]
depth_out = snakemake.output[0]

df = pd.read_csv(depth_in, sep="\t", compression="gzip")
df["start"] = (df["window"] // 1000) * 1000
df = df.groupby(["contig", "start"])["coverage"].sum().rename("tcov").reset_index()
df["stop"] = df["start"] + 1000

names = ["contig", "length", "mapped", "unmapped"]
idxstatdf = pd.read_csv(idxstats, sep="\t", names=names, skipfooter=1, usecols=[0, 1])
csize = idxstatdf.set_index("contig")["length"].to_dict()

df["csize"] = df["contig"].map(csize)
df["stop"] = df[["stop", "csize"]].min(axis=1)
df = df.drop("csize", axis=1)
df["mcov"] = df["tcov"] / (df["stop"] - df["start"])

df.to_csv(depth_out, sep="\t", compression="gzip", index=False)
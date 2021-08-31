# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-08-30 10:12:19
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-08-30 15:49:34

import pandas as pd
import numpy as np

depthfile = snakemake.input[0]
idxstat   = snakemake.input[1]
depthw    = snakemake.output[0]

def iterate_depth(fname, window=50, chunksize=100000):
    kwargs = {
        "sep": "\t",
        "compression": "gzip",
        "names": ["contig", "position", "coverage"],
        "chunksize": chunksize
    }
    
    with pd.read_csv(fname, ** kwargs) as reader:
        for sdf in reader:
            sdf["window"] = (sdf["position"] // window) * window
            sdf = sdf.groupby(["contig", "window"])["coverage"].sum()
            sdf = sdf.reset_index()
            yield sdf

def yield_window_idxstat(idxstat, window=50):
    names = ["contig", "length", "mapped", "unmapped"]
    idxstatdf = pd.read_csv(idxstat, sep="\t", names=names, skipfooter=1, usecols=[0, 1])

    for contig, length in zip(idxstatdf["contig"], idxstatdf["length"]):
        windows = np.arange(0, length + 1, window)
        windows = pd.Series(windows).rename("window").to_frame()
        windows["contig"] = contig
        windows = windows[["contig", "window"]]
        yield windows

df = pd.concat(iterate_depth(depthfile))
df = df.groupby(["contig", "window"])["coverage"].sum()
df = df.reset_index()
  
windows = pd.concat(yield_window_idxstat(idxstat, 50))
df = windows.merge(df, on=["contig", "window"], how="left").fillna(0)
df.to_csv(depthw, index=False, compression="gzip", sep="\t")
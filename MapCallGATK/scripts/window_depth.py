# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2020-12-09 09:56:48
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-03-10 14:52:05

import gzip
from collections import defaultdict
import pandas as pd

WINDOWSIZE = 10000
QUANTILE = .95 # < 1

depthfile = snakemake.input[0]
outfile   = snakemake.output[0]

class Counter() :
    def __init__(self) :
        self.value = 0
        self.count = 0

    def add(self, value) :
        self.value += value
        self.count += 1

counts =  defaultdict(lambda : defaultdict(Counter))

with gzip.open(depthfile) as f :
    for line in f :
        line = line.decode("utf-8")
        line = line.strip().split()

        contig, position, coverage = line
        position, coverage = int(position), int(coverage)

        window = position // WINDOWSIZE
        counts[contig][window].add(coverage)

data = [{"contig" : contig, "window" : window, "tcov" : counter.value, "tpos" : counter.count}
        for contig, wvalues in counts.items() for window, counter in wvalues.items()]

df = pd.DataFrame(data)
df["mean coverage"] = df["tcov"] / df["tpos"]

if QUANTILE :
	quantile = df["mean coverage"].quantile(QUANTILE)
	df.loc[df["mean coverage"] > quantile, "mean coverage"] = quantile

df.to_csv(outfile, sep="\t")

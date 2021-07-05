# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-07-02 15:36:36
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-07-02 16:49:15

import os, glob
import pandas as pd
from collections import defaultdict

bname = os.path.basename
dname = os.path.dirname

checkfile = snakemake.input[0]

fnames = os.path.join(os.path.dirname(checkfile), "*.barrnap.txt")
fnames = glob.glob(fnames)

barrdata = {}

for fname in sorted(fnames) :
    print ("Read fname : %s" %(fname))
    mID = bname(fname)[:-12]
    barrdata[mID] = defaultdict(int)

    with open(fname) as f :
        for line in f :
            if line.startswith("#") : continue
            values = line.strip().split("\t")[-1]
            values = [value.split("=") for value in values.split(";")]
            values = dict(values)
            if "partial" in values["product"] : continue
            barrdata[mID][values["Name"]] += 1

found = sorted(set.union(* [set(keys) for keys in barrdata.values()]) | {"5S_rRNA", "16S_rRNA", "23S_rRNA"})
df = pd.DataFrame([{"ID" : mID, ** {key : values[key] for key in found}}
    for mID, values in barrdata.items()])

outfile = snakemake.output[0]
df.to_csv(outfile, sep="\t")
# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-07-02 15:36:26
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-07-05 16:59:48

import os, glob
import pandas as pd

bname = os.path.basename
dname = os.path.dirname

checkfile = snakemake.input[0]

fnames = os.path.join(os.path.dirname(checkfile), "*.trnascanse.*.txt")
fnames = glob.glob(fnames)

tRNAStrict = set([
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        ])

"""
tRNATCount : Total number of tRNA including undefined results
tRNACount : Number of tRNA without undefined results but with redundancy
tRNAUnique : Unique number of tRNA including non-strict tRNA
tRNAUniqueStrict : Unique number of tRNA considering only strict ones (the one which should be considered for MIMAG analysis)
"""

trnadata = []

for fname in sorted(fnames) :
    mID  = bname(fname)[:-15]
    kind = bname(fname).split(".")[-2]

    df = pd.read_csv(fname, sep="\t", usecols=[0,4], names=["name", "RNAType"], skiprows=3)
    totalu = len(df)
    
    df = df[df["RNAType"] != "Undet"] # pseudo rna
    
    if df.empty :
        total = nunique = nunique_strict = 0

    else :
        total, nunique = len(df), df["RNAType"].nunique()
        nunique_strict = len(set(df["RNAType"].apply(lambda name : name.upper())) & tRNAStrict)
    
    trnadata.append(
        {"ID" : mID,
        "tRNAKind" : kind,
        "tRNACount" : total, 
        "tRNAUnique" : nunique, 
        "tRNAUniqueStrict" : nunique_strict, 
        "tRNATCount" : totalu}
        )
    
df = pd.DataFrame(trnadata)
outfile = snakemake.output[0]
df.to_csv(outfile, sep="\t")


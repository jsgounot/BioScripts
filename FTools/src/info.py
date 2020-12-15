#!/usr/bin/python
import os, sys
import argparse

import pandas as pd
pd.options.display.max_rows = 999
pd.set_option('display.width', 1000)

import utils

def seqinfo(seq) :
    seq = seq.upper()
    gc, n, length = seq.count("G") + seq.count("C"), seq.count("N"), len(seq)
    return {"Length" : length, "#GC" : gc, "#N" : n}

def run(files, sort="ID", head=0, ** fkwargs) :
    for fname, fdata in utils.iter_fdata(files, ** fkwargs) :
        fdata = ({"ID" : record.id, ** seqinfo(record.seq)} for record in fdata)

        df = pd.DataFrame(fdata)

        if not df.empty :
            df = df.sort_values(sort)
            
            s = df.sum(axis=0)
            s["ID"] = "Total"
            
            if head : df = df.head(head)

            df = df.append(s, ignore_index=True)

            df["%N"] = df["#N"] * 100 / df["Length"]
            df["%GC"] = df["#GC"] * 100 / df["Length"]

            df = df[["ID", "Length", "%GC", "%N"]]
            df["Length"] = df["Length"].apply(lambda x : '{:,}'.format(x))

            print (fname)
            print (df)
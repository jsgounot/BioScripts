#!/usr/bin/python
import os, sys
import argparse
from collections import Counter
import pandas as pd
import utils

def run_basecount(files, sort="ID", head=0, nsep=False, ** fkwargs) :
    for fname, fdata in utils.iter_fdata(files, ** fkwargs) :
        fdata = ({"ID" : record.id, ** Counter(record.seq)} for record in fdata)
        df = pd.DataFrame(fdata).fillna(0)

        for column in df.columns :
            if column == "ID" : continue
            df[column] = df[column].astype(int)
        
        columns = sorted([column for column in df.columns if column != "ID"])
        columns.insert(0, "ID")
        df = df[columns]

        df = df.sort_values(sort)
        s = df.sum(axis=0)
        s["ID"] = "Total"
            
        if head : df = df.head(head)
        df = df.append(s, ignore_index=True)

        print (fname)
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  
            print (df)

def run_default(files, sort="ID", head=0, nsep=False, ** fkwargs) :
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
            if nsep : 
                df["Length"] = df["Length"].apply(lambda x : '{:,}'.format(x))

            print (fname)
            with pd.option_context('display.max_rows', None, 'display.max_columns', None):  
                print (df)


def seqinfo(seq) :
    seq = seq.upper()
    gc, n, length = seq.count("G") + seq.count("C"), seq.count("N"), len(seq)
    return {"Length" : length, "#GC" : gc, "#N" : n}
 
def run(files, sort="ID", head=0, basecount=False, nsep=False, ** fkwargs) :
    fun = run_basecount if basecount else run_default
    fun(files, sort, head, nsep, ** fkwargs)
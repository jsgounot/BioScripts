#!/usr/bin/python
import os, sys
import argparse
from collections import Counter
import pandas as pd
import utils

def run_basecount(files, sort="ID", head=0, nsep=False, ** fkwargs) :
    for fname, fdata in utils.iter_fdata(files, ** fkwargs) :
        fdata = ({"ID" : record.id, ** Counter(record.seq.upper())} for record in fdata)
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

        yield fname, df

def run_default(files, sort="ID", head=0, nsep=False, ** fkwargs) :
    mdf = None
    for fname, fdata in utils.iter_fdata(files, ** fkwargs) :
        fdata = ({"ID" : record.id, ** seqinfo(record.seq)} for record in fdata)
        df = pd.DataFrame(fdata)

        if not df.empty :
            df = df.sort_values(sort)
            
            s = df.sum(axis=0)
            s["ID"] = "Total"
            
            if head : 
                df = df.head(head)


            df = df.append(s, ignore_index=True)

            df["%N"] = df["#N"] * 100 / df["Length"]
            df["%GC"] = df["#GC"] * 100 / df["Length"]

            df = df[["ID", "Length", "%GC", "%N"]]
            if nsep : 
                df["Length"] = df["Length"].apply(lambda x : '{:,}'.format(x))

            yield fname, df

def seqinfo(seq) :
    seq = seq.upper()
    gc, n, length = seq.count("G") + seq.count("C"), seq.count("N"), len(seq)
    return {"Length" : length, "#GC" : gc, "#N" : n}
 
def run(files, sort="ID", head=0, basecount=False, nsep=False, asone=False, fullname=False, outfile=None, ** fkwargs) :
    fun = run_basecount if basecount else run_default
    asone = True if outfile else asone
    mdf = None
    for fname, df in fun(files, sort, head, nsep, ** fkwargs):
        if asone:
            df['fname'] = fname if fullname else os.path.basename(fname)
            df = df[['fname'] + list(df.columns)[:-1]]
            mdf = df if mdf is None else pd.concat((mdf, df)).fillna(0)

        else :           
            print (fname)
            with pd.option_context('display.max_rows', None, 'display.max_columns', None):  
                print (df)

    if asone and mdf is not None and basecount:
        bases = sorted(set(mdf.columns) - {'ID', 'fname'})
        for base in bases: mdf[base] = mdf[base].astype(int)

    if asone and mdf is not None:
        mdf = mdf[mdf["ID"] != "Total"]
        mdf = mdf.reset_index(drop=True)

    if asone and mdf is None:
        print ("Empty table, nothing to be print ...")

    elif asone and outfile:
        mdf.to_csv(outfile, sep="\t", index=False)

    elif asone:
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  
            print (mdf)
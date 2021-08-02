# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-07-02 15:50:09
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-08-02 12:53:00

import os, glob, re
import pandas as pd

bname = os.path.basename
dname = os.path.dirname

def touch(fname):
    with open(fname, "w") as f:
        pass

def load_eval(evalFile) :
    data = []
    
    # Wanted column index
    columns = {
        "ID" : 1,
        "Completness" : 6,
        "Contamination" : 7,
        "Size" : 9,
        "#Contigs" : 11,
        "Contig_N50" : 14,
        "Longest_contig" : 18
        }
    
    with open(evalFile) as f :
        for line in f :
            line = re.split(" \s+", line)
            if len(line) == 1 or line[1] == "Bin Id" : continue
            data.append({column : line[idx] for column, idx in columns.items()})
                                       
    assert len(evalFile)
                                                
    df = pd.DataFrame(data)
    
    def set_completness(row) :
        completness = row["Completness"]
        contamination = row["Contamination"]
        if completness > 90 and contamination < 5 :
            return "HIGH"
        elif completness >= 50 and contamination <= 10 :
            return "MEDIUM"
        else :
            return "LOW"
                                                                                                
    df["Contamination"] = df["Contamination"].astype(float)
    df["Completness"] = df["Completness"].astype(float)
                                                                                                
    fun = lambda row : set_completness(row)
    df["CheckMStatus"] = df.apply(lambda row : set_completness(row), axis=1)

    return df

fname = snakemake.input[0]
outfile = snakemake.output[0]

if os.path.basename(fname) == "check.empty":
    touch(outfile)

else:
    df = load_eval(fname)
    df.to_csv(outfile, sep="\t")
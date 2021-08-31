# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-07-09 12:39:20
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-07-09 12:49:59

import os
import pandas as pd

fnames = snakemake.input
result = []

for fname in fnames :
    name = os.path.splitext(os.path.basename(fname))[0]
    initial_file = os.path.join("./input", name + ".fastq.gz")

    spring_size = os.path.getsize(fname)
    intial_size = os.path.getsize(initial_file)

    result.append({
        "name" : name,
        "initial_size" : intial_size,
        "spring_size" : spring_size
        })  

df = pd.DataFrame(result)

def sizeof_fmt(num, suffix='B'):
    # https://stackoverflow.com/questions/1094841/get-human-readable-version-of-file-size
    for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)

df["gain"] = df["initial_size"] - df["spring_size"]
df["ratio"] = df["initial_size"] / df["spring_size"]

df["initial_size"] = df["initial_size"].apply(sizeof_fmt)
df["spring_size"] = df["spring_size"].apply(sizeof_fmt)
df["gain"] = df["gain"].apply(sizeof_fmt)

outfile = snakemake.output[0]
df.to_csv(outfile, sep="\t", index=False)
# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-07-09 12:39:20
# @Last Modified by:   jsgounot
# @Last Modified time: 2024-09-25 20:27:44

import os
import pandas as pd

fnames = snakemake.input
result = []

for fname in fnames :
    
    name = os.path.basename(fname)[:-11]

    with open(fname) as f:
        lines = [line.strip().split() for line in f]
        assert len(lines) == 3 

        spring_size = intial_size = 0
        spring = False

        for line in lines:
            if line[1].endswith('.spring'):
                assert spring == False
                spring = True
                spring_size = int(line[0])
            else:
                intial_size += int(line[0])

        assert spring_size != 0
        assert intial_size != 0

    result.append({
        "name" : name,
        "initial_size" : intial_size,
        "spring_size" : spring_size
        })  

df = pd.DataFrame(result)

def sizeof_fmt(num, suffix='B'):
    # https://stackoverflow.com/questions/1094841/get-human-readable-version-of-file-size
    for unit in ['Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
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
# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2024-09-29 22:26:36
# @Last modified by:   jsgounot
# @Last Modified time: 2024-09-30 14:23:38

import os, glob, gzip, sys, json
import pandas as pd

bname = os.path.basename
dname = os.path.dirname

fnames = sys.argv[1:]

jdata = []
for fname in fnames:
    name = bname(fname)
    * name, mapper, ext = name.split('.')
    name = '.'.join(name)
    
    with open(fname) as f:
        fdata = json.load(f)

    fdata['name'], fdata['mapper'] = name, mapper
    jdata.append(fdata)

df = pd.DataFrame(jdata)

df['%fastp_base'] = df['fastp_bases'] * 100 / df['raw_bases']
df['%decont_bases'] = df['decont_bases'] * 100 / df['raw_bases']
df['decont_diff'] =  df['%fastp_base'] - df['%decont_bases']

print (df)
# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2024-09-26 09:20:13
# @Last modified by:   jsgounot
# @Last Modified time: 2024-09-26 15:51:48

import os, glob, gzip
import pandas as pd

bname = os.path.basename
dname = os.path.dirname

fnames = snakemake.input
outfile = snakemake.output[0]

res = []
for fname in fnames:
	with open(fname) as f:
		if any(line.strip() for line in f):
			res.append([fname, False])
		else:
			res.append([fname, True])

with open(outfile, 'w') as f:
	f.write(f'sample\tridx\tcheck\n')
	for fname, out in res:
		* name, ridx, _ = bname(fname).split('.')
		name = '.'.join(name)
		f.write(f'{name}\t{ridx}\t{out}\n')
		if not out:
			print (f'Warning: It seems that compression leads to a different output for: {name} - {ridx} ({fname})')

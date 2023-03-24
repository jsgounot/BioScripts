# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-03-24 14:13:34
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-03-24 14:21:23

import os
import pandas as pd

fname = snakemake.input[0]
df = pd.read_csv(fname, sep='\t', skiprows=1)
df = df[df['refseq_category'] == 'representative genome']

outfile = snakemake.output[0]
with open(outfile, 'w') as f:
	for line in df['ftp_path']:
		bname = os.path.basename(line)
		line = os.path.join(line, bname + '_genomic.fna.gz')
		f.write(line + '\n')
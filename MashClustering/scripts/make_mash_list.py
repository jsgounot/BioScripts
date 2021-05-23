# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-05-23 16:22:58
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-05-23 16:43:50

fnames  = snakemake.params.fnames
outfile = snakemake.output[0]

with open(outfile, "w") as f :
	for fname in fnames :
		f.write(fname + "\n")
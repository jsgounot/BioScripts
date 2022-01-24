# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-09-22 11:02:36
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-11-19 18:41:47

import json

outfile = snakemake.output[0]
fnames = snakemake.params[0]

fnames = {
	idx: fname
	for idx, fname in enumerate(fnames)
}

with open(outfile, 'w') as f :
	json.dump(fnames, f, indent=4)
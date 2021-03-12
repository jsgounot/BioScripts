# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-03-12 14:13:53
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-03-12 14:58:21

import os
import pandas as pd

fname   = snakemake.input[0]
outfile = snakemake.output[0]

class RName :
	def __init__(self) :
		self.d = {}

	def get(self, name) :
		value = self.d.get(name, None)
		if value is None : 
			value = os.path.splitext(os.path.basename(name))[0]
			self.d[name] = value
		return value

rname = RName()
gname = lambda name : rname.get(name)

names = ["rID", "qID", "distance", "p-value", "shared-hashes"]
df = pd.read_csv(fname, names=names, compression="gzip", sep="\t")

df["rID"] = df["rID"].map(rname.get)
df["qID"] = df["qID"].map(rname.get)

df = pd.pivot_table(df, index="rID", columns="qID", values="distance")
df.to_csv(outfile, compression="gzip", index_label="")
# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-07-02 16:51:41
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-08-02 13:47:45

import os
import pandas as pd

rmext = lambda fname : os.path.splitext(fname)[0]

checkm   = snakemake.input.checkm
barrnap  = snakemake.input.barrnap
trnascan = snakemake.input.trnascan

def isfileempty(fname):
	return os.stat(fname).st_size == 0

def load_df(fname) :
	df = pd.read_csv(fname, sep="\t", index_col=0)
	df["sample"] = os.path.basename(os.path.dirname(fname))
	return df

checkm   = pd.concat((load_df(fname) for fname in checkm if not isfileempty(fname)))
barrnap  = pd.concat((load_df(fname) for fname in barrnap if not isfileempty(fname)))
trnascan = pd.concat((load_df(fname) for fname in trnascan if not isfileempty(fname)))

barrnap["ID"] = barrnap["ID"].apply(rmext)

# Need to be applied twice to remove bother extension
# and trnakind name
trnascan["ID"] = trnascan["ID"].apply(rmext)
trnascan["ID"] = trnascan["ID"].apply(rmext)

trnascan = trnascan[["ID", "tRNAKind", "sample", "tRNAUniqueStrict"]]

# We select the tRNA scan results which display the highest number of tRNA
# this is not the best but the alternative would be to identify
# whether each bin is archea or bacteria, which would need phylogenetic placement
trnascan = trnascan[trnascan.groupby(["sample", "ID"])["tRNAUniqueStrict"].transform(max) == trnascan["tRNAUniqueStrict"]]

# In case of duplicates (# archea == # bacterial), we select the bacterial one
# Change nothing, just make more sense to assign bacteria
trnascan = trnascan.sort_values(["sample", "ID", "tRNAKind"])
trnascan = trnascan.drop_duplicates(subset=["sample", "ID"], keep="last")

df = checkm.merge(barrnap, on=["sample", "ID"], how="left")
df = df.merge(trnascan, on=["sample", "ID"], how="left")

df["RNAPASS"] = (df["5S_rRNA"] > 0) & (df["16S_rRNA"] > 0) & (df["23S_rRNA"] > 0) & (df["tRNAUniqueStrict"] >= 18)

def get_mimag(row) :
	cstat = row["CheckMStatus"]
	tstat = row["RNAPASS"]

	if cstat == "HIGH" :
		if not tstat : cstat = "MEDIUM"

	return cstat

df["MIMAG"] = df.apply(get_mimag, axis=1)

outfile = snakemake.output.bins
df.to_csv(outfile, sep="\t", index=False)

# MAGs softlink (MEDIUM AND HIGH)
outdir = "data/mags/"
os.makedirs(outdir, exist_ok=True)
cwd = os.getcwd()

df = df[df["MIMAG"].isin(("MEDIUM", "HIGH"))]

for sample, binid in zip(df["sample"], df["ID"]) :
	fname = os.path.join(cwd, "data/bins/", sample, binid + ".fa")
	outfname = os.path.join(outdir, binid + ".fa")
	if not os.path.isfile(outfname) : os.symlink(fname, outfname)

outfile = snakemake.output.mags
df.to_csv(outfile, sep="\t", index=False)
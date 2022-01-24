# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-12-16 11:58:52
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-12-16 12:03:54

from Bio import SeqIO

fname = snakemake.input[0]
fdata = list(SeqIO.parse(fname, 'fasta'))

sample = snakemake.wildcards['mid']

print (sample)

for record in fdata:
	nname = record.id.strip().split()[0]
	nname = sample + '_' + nname
	record.id = record.name = nname

fname = snakemake.output[0]
SeqIO.write(fdata, fname, 'fasta')
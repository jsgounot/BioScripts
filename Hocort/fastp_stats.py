# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2024-09-27 18:29:06
# @Last modified by:   jsgounot
# @Last Modified time: 2024-10-09 09:40:27

import json

first = snakemake.input['first']
second = snakemake.input['second']

res = {}

with open(first) as f:
	jdata = json.load(f)
	res['raw_bases'] = jdata['summary']['before_filtering']['total_bases']
	res['fastp_bases'] = jdata['summary']['after_filtering']['total_bases']
	res['raw_r1_n'] = jdata['read1_before_filtering']['total_reads']
	res['raw_r2_n'] = jdata['read2_before_filtering']['total_reads']
	res['fastp_r1_n'] = jdata['read1_after_filtering']['total_reads']
	res['fastp_r2_n'] = jdata['read2_after_filtering']['total_reads']

with open(second) as f:
	jdata = json.load(f)

	left = jdata['summary']['before_filtering']['total_bases']
	right = jdata['summary']['after_filtering']['total_bases']
	
	if left != right:
		prc = (left - right) * 100 / left
		print (f'WARNING: The number of bases before and after the second fastp round is different')
		print (f'WARNING: Before: {left:,} - After {right:,} - % lost: {prc:.2f}')
		print (f'WARNING: In theory, this should not happen: Fastp should not find any new reason to clean after decontamination')
		print (f'WARNING: Nevertheless, this can happen and can be ignored if the prc is low')
		print (f'WARNING: Will keep the number of bases AFTER second filtering: {right:,}')
		print (f'WARNING: This pipeline will crash if % lost > 1%')
		if prc > 1: raise Exception('Prc lost > 1% - see reason in log')

	res['decont_bases'] = jdata['summary']['after_filtering']['total_bases']
	res['decont_r1_n'] = jdata['read1_after_filtering']['total_reads']
	res['decont_r2_n'] = jdata['read2_after_filtering']['total_reads']

outfile = snakemake.output[0]
with open(outfile, 'w') as f:
	json.dump(res, f, indent=4, sort_keys=True)
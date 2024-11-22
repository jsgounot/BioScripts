# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2024-11-20 18:22:51
# @Last modified by:   jsgounot
# @Last Modified time: 2024-11-21 14:37:39

import os, glob, gzip, sys
import argparse, random

bname = os.path.basename
dname = os.path.dirname

parser = argparse.ArgumentParser(
    prog='extract_reads.py',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('output', help='Kraken fastq file')
parser.add_argument('taxid', help='Target taxonomic ID. Unclassified is 0', type=int)
parser.add_argument('fastq', nargs='*', help='Fastq file(s)')
parser.add_argument('--paired', action='store_true', help='If paired-end fastq, return reads only if both pairs are assigned to taxid')
parser.add_argument('--strict', action='store_true', help='Raise an error if one read is not found')
parser.add_argument('--count', nargs='?', help='Limit the number of output reads to n', type=int, default=0)
parser.add_argument('--seed', nargs='?', help='Random seed number if case of count', type=int, default=1)
parser.add_argument('--quiet', action='store_true', help='Might reduce RAM usage')
args = parser.parse_args()

# python /home/users/astar/gis/ericejs/scratch/softs/BioScripts/Kraken2/run_pipeline/extract_reads.py kraken2/new_hocort/reports/WHS629.hisat2/HRGMFungiHS.kraken.output.tsv.gz 7062 temp/WHS629.hisat2.spring.r1.fq.gz

smart_opener = lambda fname: gzip.open(fname, 'rt') if fname.endswith('.gz') else open(fname)
random.seed(args.seed)

if not os.path.isfile(args.output):
    raise Exception(f'File missing: {args.output}')

for fastq in args.fastq:
    if not os.path.isfile(fastq):
        raise Exception(f'File missing: {fastq}')

# --------------------------------------------------------------------------------

def extract_rids(fname, taxid):
    with smart_opener(fname) as f:
        for line in f:
            line = line.strip().split('\t')
            rid, l_taxid = line[1], int(line[2])
            if l_taxid == taxid:
                yield rid

def parse_rids_paired(fname, taxid, quiet):
    iterator = extract_rids(fname, taxid)
    left, right = set(), set()
    for rid in iterator:
        if rid.endswith('/1'): container = left
        elif rid.endswith('/2'): container = right
        else: raise Exception(f'rid does not end with /1 or /2: {rid}')
        container.add(rid[:-2])

    merged = left & right
    if not quiet:
        print (f'Pair 1 reads count: {len(left)}', file=sys.stderr)
        print (f'Pair 2 reads count: {len(right)}', file=sys.stderr)
        print (f'Merged reads count: {len(merged)}', file=sys.stderr)

    for rid in merged:
        yield rid + '/1'
        yield rid + '/2'

def parse_rids(fname, taxid, paired, count, quiet):
    iterator = parse_rids_paired(fname, taxid, quiet) if paired else extract_rids(fname, taxid)
    rids = set(iterator)    

    if not args.quiet:
        msg = f'Found {len(rids):,} reads'
        print(msg, file=sys.stderr)

    return set(random.sample(sorted(rids), k=count)) if count else rids

def iter_fastq(fastq, rids):
    record, found = False, set()
    with smart_opener(fastq) as f:
        for idx, line in enumerate(f):
            if idx % 4 == 0:
                rid = line.strip()[1:]
                record = rid in rids
                if record: found.add(rid)            

            if record:
                print(line.strip(), file=sys.stdout)

    return found

# --------------------------------------------------------------------------------

if not args.quiet:
    print (f'Search reads ID for taxid {args.taxid}', file=sys.stderr)

rids = parse_rids(args.output, args.taxid, args.paired, args.count, args.quiet)

found = set()
for fastq in args.fastq:
    found |= iter_fastq(fastq, rids)

missing = rids - found
if args.strict and missing:
    raise Exception(f'{len(missing):,} reads were not found in the fastq files, e.g: {next(iter(missing))}')





# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-06-29 14:28:14
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-07-05 09:24:15

import argparse
import pysam
import sys
import json
import gzip

'''
python compare_bam.py bamcompare --bam1 test/ERR7671878.ahadrus.spmp.bam --bam2 test/ERR7671878.ahadrus.uhgg.bam | python compare_bam.py readextract --rids - --reads /home/users/astar/gis/ericejs/scratch/projects/spmp/data/reads/spring/temp/ERR7671878_
2.fastq.gz --outreads1 ./test.reads.fastq.gz

python compare_bam.py bamcompare --bam1 test/ERR7671878.ahadrus.spmp.bam --bam2 test/ERR7671878.ahadrus.uhgg.bam --head | python compare_bam.py readstat --rids -
'''

def extract_readids(bamfile, min_coverage, min_identity, head=False):
    sys.stderr.write(f'Parse bamfile: {bamfile} ...\n')

    idthreshold = 1 - min_identity
    bamfile = pysam.AlignmentFile(bamfile, "rb")
    rids  = set()
    total = set()

    for read in bamfile.fetch(until_eof=True):

        total.add(read.query_name)

        if read.infer_query_length() != None and read.is_proper_pair:
            len_ops, num_ops = read.get_cigar_stats()

            # (M + I + D / read_length) and NM / (M + I + D)
            mid = sum(len_ops[0:3])

            if mid / read.infer_read_length() >= min_coverage and len_ops[10] / mid <= idthreshold:
                rids.add(read.query_name)

        if head and len(total) > 100000:
            break

    bamfile.close()     
    sys.stderr.write(f'Found {len(rids):,} reads passing filters - total {len(total):,} - prc {((len(rids) * 100) / len(total)):.2f}\n')
    return rids

def parse_reads(fnames, rids):
    for fname in fnames:
        sys.stderr.write(f'Parse fasta file: {fname} ...\n')
        with gzip.open(fname, 'rt') as f:
            for idx, line in enumerate(f):
                if idx % 4 == 0:
                    identifier = line.strip().split()[0][1:]
                    found = identifier in rids
                    strand = identifier[-1]
                if found:
                    yield strand, line


def write_single(reads, outreads1, rids):
    found = set()

    with gzip.open(outreads1, 'wt') as of:
        for strand, line in parse_reads(reads, rids):
            of.write(line)

def write_paired(reads, outreads1, outreads2, rids):
    found = set()
    print (rids)

    with open(outfile, 'w') as of:
        for read in reads:
            with gzip.open(read, 'rt') as f:
                for line in f:
                    identifier = line.strip().split()[0][1:]
                    record = identifier in found
                    if record: found.add(identifier)

                    print (line)
                    exit()

                    for i in range(4):
                        line = next(f)

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help')
    
    sub_bamcompare = subparsers.add_parser('bamcompare', help='Compare bam and extract reads IDS')
    sub_bamcompare.add_argument("--bam1", dest="bam1", type=str, help="input bam file 1", required=True)
    sub_bamcompare.add_argument("--bam2", dest="bam2", type=str, help="input bam file 2", required=True)
    sub_bamcompare.add_argument("--min_coverage", dest="min_coverage", default=0.90, type=float, help="a minimum read coverage [0,1] (default: 0.9)", required=False)
    sub_bamcompare.add_argument("--min_identity", dest="min_identity", default=0.99, type=float, help="a minimum alignement identity [0,1] (default: 0.99)", required=False)
    sub_bamcompare.add_argument("--head", dest="head", action='store_true', help="just use the first 100k bam alignements", required=False)

    sub_stat = subparsers.add_parser('readstat', help='Extract statistic based on bamparser IDs')
    sub_stat.add_argument("--rids", dest="rids", default='-', type=str, help="reads IDs file (default stdin)", required=True)

    sub_readextract = subparsers.add_parser('readextract', help='Extract reads based on bamparser IDs')
    sub_readextract.add_argument("--rids", dest="rids", default='-', type=str, help="reads IDs file", required=True)
    #sub_readextract.add_argument("--reads", dest="reads", nargs='+', help="reads file", required=True)
    #sub_readextract.add_argument("--outreads1", dest="outreads1", default='', type=str, help="output for pair1", required=True)
    #sub_readextract.add_argument("--outreads2", dest="outreads2", default='', type=str, help="output for pair2", required=False)
    sub_readextract.add_argument("--kind", dest="kind", default='intersection', type=str, choices=['bam1', 'bam2', 'intersection', 'union'] , help="reads group", required=False)

    args = parser.parse_args()

    if 'bam1' in args:
        bamcompare(args)

    elif 'rids' in args and 'kind' in args:
        readextract(args)

    elif 'rids' in args:
        readstat(args)

def bamcompare(args):
    reads_id_1 = ri1 = extract_readids(args.bam1, args.min_coverage, args.min_identity, args.head)
    reads_id_2 = ri2 = extract_readids(args.bam2, args.min_coverage, args.min_identity, args.head)

    jdata = {
        'bam1': list(ri1),
        'bam2': list(ri2)
    }

    json.dump(jdata, sys.stdout)

def readstat(args):
    if args.rids == '-':
        jdata = str(sys.stdin.read())
        jdata = json.loads(jdata)

    else:
        with open(args.rids, 'r') as f:
            jdata = json.load(f)

    ri1 = set(jdata['bam1'])
    ri2 = set(jdata['bam2'])

    print ('Statistics ...')
    print (f'- Intersection: {len(ri1 & ri2):,}')
    print (f'- Union: {len(ri1 | ri2):,}')
    print (f'- Unique R1: {len(ri1 - ri2):,}')
    print (f'- Unique R2: {len(ri2 - ri1):,}')

def readextract(args):
    if args.rids == '-':
        jdata = str(sys.stdin.read())
        jdata = json.loads(jdata)

    else:
        with open(args.rids, 'r') as f:
            jdata = json.load(f)

    if args.kind == 'bam1':
        rids = set(jdata['bam1'])
    elif args.kind == 'bam2':
        rids = set(jdata['bam2'])
    elif args.kind == 'intersection':
        rids = set(jdata['bam1']) & set(jdata['bam2'])
    elif args.kind == 'union':
        rids = set(jdata['bam1']) | set(jdata['bam2'])
    else:
        raise Exception()

    sys.stdout.write('\n'.join(('^@' + rid for rid in rids)))

if __name__ == '__main__':
    main()
"""
Long Read Connectivity SNP
https://github.com/jsgounot/BioScripts/tree/main/LRCS
"""

import os, gzip
from collections import defaultdict
from itertools import combinations

import click
import pysam
import pandas as pd

try : import tqdm; USETQDM = True
except ImportError : USETQDM = False

def extract_snps(vcffile) :
    """
    Extract SNPs positions based on a vcf file
    Arguments:
        vcffile {str} -- VCF file path
    """
    print (f"Read vcf file : {vcffile}")

    data = defaultdict(lambda : defaultdict(set))
    vcf = pysam.VariantFile(vcffile, mode="rb")
    for idx, record in enumerate(vcf.fetch()) :
        ref = record.ref[0]
        contig, pos = record.contig, record.pos
        alts = record.alts

        for alt in alts :
            # deletion of insertion (then ref[0] can be equals to alt)
            if len(alt) > 1 or alt == ref : continue
            data[contig][pos].add(alt)

    for contig, sites in data.items() :
        nsite = len(sites)
        nsnp = sum(len(snps) for snps in sites.values())
        print (f"Contig {contig} : {nsite} sites and {nsnp} snps found")

    return data

def get_windows(snps, bamfile) :
    """
    Select SNPs based on the rule descibe above
    Arguments:
        snps {dict} -- SNPs position and allele provided by extract_snps
        bamfile {str} -- Bam file path
    """
    print (f"Read bam file file : {bamfile}")

    results = defaultdict(list)
    bamdata = pysam.AlignmentFile(bamfile, "rb")
    for contig, sites in snps.items() :
        snps_positions = set(sites)
        found = set()

        if USETQDM : iterator = tqdm.tqdm(bamdata.fetch(contig), total=bamdata.count(contig))
        else : iterator = bamdata.fetch(contig)

        for read in iterator :
            
            # we search for SNPs position
            read_positions = set(read.get_reference_positions())
            overlapp = read_positions & snps_positions
            if len(overlapp) < 2 : continue

            try :
                if results[contig][-1] & overlapp :
                    results[contig][-1] |= overlapp 
                else :
                    results[contig].append(overlapp)
            except IndexError :
                results[contig].append(overlapp)

    return results

def merge(contig, sets):
    # https://stackoverflow.com/questions/9110837/python-simple-list-merging-based-on-intersections
    print (f"Merge contig for : {contig}")

    merged = True
    while merged:
        merged = False
        results = []
        while sets:
            common, rest = sets[0], sets[1:]
            sets = []
            for x in rest:
                if x.isdisjoint(common):
                    sets.append(x)
                else:
                    merged = True
                    common |= x
            results.append(common)
        sets = results
    return sets

def merge_windows(cwindows) :
    # Merge overlapping windows in case of unsorted bam
    return {contig : merge(contig, windows) for contig, windows in cwindows.items()}

def w2df(cwindows) :
    data = ({"Contig" : contig, "Start" : min(window), "Stop" : max(window), "#SNPs" : len(window)}
            for contig, windows in cwindows.items() for window in windows)
    
    df = pd.DataFrame(data)
    df["Size"] = df["Stop"] - df["Start"]

    return df

@click.command()
@click.option("--vcf", required=True, type=str)
@click.option("--bam", required=True, type=str)
@click.option("--sortedbam", required=True, is_flag=True, default=False)
def run(vcf, bam, sortedbam) :
    print (f"Work with : {vcf}")
    snps = extract_snps(vcf)
    cwindows = get_windows(snps, bam)
    
    if not sortedbam :
        cwindows = merge_windows(cwindows)

    df = w2df(cwindows)

    with pd.option_context('display.max_rows', None, 'display.max_columns', None) :
        print(df)

if __name__ == "__main__" :
    run()
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

class ReadsGroup(object) :
    def __init__(self, start=None, stop=None, snps=None) :
        self.start = start
        self.stop = stop
        self.snps = snps or set()

    @property   
    def nsnps(self) :
        return len(self.snps)

    def isdisjoint(self, element) :
        return self.snps.isdisjoint(element.snps)

    @staticmethod
    def getminmax(v1, v2, fun) :
        value = None
        if v1 is not None :
            value = v1
        if v2 is not None :
            if value is None :
                value = v2
            else :
                value = fun((v1, v2))
        return value

    def __ior__(self, other) :
        start = ReadsGroup.getminmax(self.start, other.start, min)
        stop  = ReadsGroup.getminmax(self.stop, other.stop, max)
        snps  = self.snps | other.snps
        return ReadsGroup(start, stop, snps)

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

def iter_reads(snps, bamfile) :
    print (f"Read bam file file : {bamfile}")
    bamdata = pysam.AlignmentFile(bamfile, "rb")

    for contig, sites in snps.items() :
        snps_positions = set(sites)

        if USETQDM : iterator = tqdm.tqdm(bamdata.fetch(contig), total=bamdata.count(contig))
        else : iterator = bamdata.fetch(contig)

        for read in iterator :
            # we search for SNPs position
            read_positions = set(read.get_reference_positions())
            overlapp = read_positions & snps_positions
            if len(overlapp) < 2 : continue

            name = read.query_name

            start = min(read_positions)
            stop  = max(read_positions)
            rg = ReadsGroup(start, stop, overlapp)

            yield contig, name, rg


def get_windows(snps, bamfile) :
    """
    Select SNPs based on the rule descibe above
    Arguments:
        snps {dict} -- SNPs position and allele provided by extract_snps
        bamfile {str} -- Bam file path
    """
    results = defaultdict(list)
    
    for contig, name, rg in iter_reads(snps, bamfile) :   
        try :
            if not results[contig][-1].isdisjoint(rg) :
                results[contig][-1] |= rg 
            else :
                results[contig].append(rg)
        except IndexError :
            results[contig].append(rg)

    return results

def get_windows_chimerics(snps, bamfile) :
    """
    Return all position for all read name
    Does not perform any merging operation between reads
    Need to use merge function for this
    """
    data = defaultdict(lambda : defaultdict(ReadsGroup))
    for contig, name, rg in iter_reads(snps, bamfile) :
        data[contig][name] |= rg

    data = {contig : list(value.values()) for contig, value in data.items()}

    return data

def merge(contig, rgs):
    # https://stackoverflow.com/questions/9110837/python-simple-list-merging-based-on-intersections
    print (f"Merge contig for : {contig}")

    merged = True
    while merged:
        merged = False
        results = []
        while rgs:
            common, rest = rgs[0], rgs[1:]
            rgs = []
            for x in rest:
                if x.isdisjoint(common):
                    rgs.append(x)
                else:
                    merged = True
                    common |= x
            results.append(common)
        rgs = results
    return rgs

def merge_windows(crgs) :
    # Merge overlapping windows in case of unsorted bam
    return {contig : merge(contig, rgs) for contig, rgs in crgs.items()}

def w2df(crgs) :
    data = ({"Contig" : contig, "Start" : rg.start, "Stop" : rg.stop, "#SNPs" : rg.nsnps,
             "MinSNP" : min(rg.snps), "MaxSNP" : max(rg.snps)}
            for contig, rgs in crgs.items() for rg in rgs)
    
    df = pd.DataFrame(data)
    df["Size"] = df["Stop"] - df["Start"]

    return df

@click.command()
@click.option("--vcf", required=True, type=str)
@click.option("--bam", required=True, type=str)
@click.option("--sortedbam", required=True, is_flag=True, default=False)
@click.option("--chimeric", required=True, is_flag=True, default=False)
def run(vcf, bam, sortedbam, chimeric) :
    print (f"Work with : {vcf}")

    snps = extract_snps(vcf)
    
    if chimeric :
        crgs = get_windows_chimerics(snps, bam)
    else :
        crgs = get_windows(snps, bam)

    if sortedbam == False or chimeric == True :
        crgs = merge_windows(crgs)

    df = w2df(crgs)

    with pd.option_context('display.max_rows', None, 'display.max_columns', None) :
        print(df)

if __name__ == "__main__" :
    run()
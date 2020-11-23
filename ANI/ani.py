import glob, os
from itertools import permutations

from Bio import SeqIO
import pandas as pd
import fire

from pyblast import BCLine6, utils

"""
A local python script to quickly calculate ANI (Average Nucleotidic Diversity)
values between samples
https://pubmed.ncbi.nlm.nih.gov/26576653/
https://pubmed.ncbi.nlm.nih.gov/17220447/
http://jspecies.ribohost.com/jspeciesws/#analyse
https://github.com/chjp/ANI

Internaly use my pyblast librairy for fast processing
https://github.com/jsgounot/PyBlast

Each fasta file should contain data for only one species
"""

def chop(fname, size=1020) :
    """
    Chop a fasta file sequences to segments of size `size` 
    into a pyblast temporary fasta file.
    """

    fdata = SeqIO.parse(fname, "fasta")

    def iterator(fdata, size) :
        for element in fdata :
            for i in range(0, len(element), size) :
                subelement = element[i:i+size]
                subelement.id = subelement.name = subelement.name + "." + str(i)
                subelement.description = ""
                if len(subelement) < 100 : continue
                yield subelement

    fdata = iterator(fdata, size)
    return utils.TMPFasta(fdata, delete=True)

def blastduo(blastn, duo, f1, f2, qlen) :
    """
    Blast two sequences and calculate the ANI value and coverage

    https://www.ncbi.nlm.nih.gov/Class/BLAST/blastallopts.txt
    https://github.com/widdowquinn/pyani/blob/8bf51a64b92d3fa787a3b2b106fda93505cba44c/pyani/anib.py
    """

    cmdline = {"xdrop_gap_final" :  150, "evalue" : "1e-15", "dust" : "no", "max_target_seqs" : 1}

    bcl = BCLine6("blastn", cmd=blastn, query=f1, subject=f2, **cmdline)
    df = bcl.run()

    df = df[df.groupby(['qseqid'])['qide'].transform(max) == df['qide']]
    df = df[df["qide"] > 30]
    df = df[df["qcov"] > 70]

    return {
        "F1" : duo[0], 
        "F2" : duo[1],
        "ANI" : df["qide"].mean(),
        "Cov" : df["qlen"].sum() * 100 / qlen
        }

def pairwise_blastn(fnames, blastn, outfile, ncore=1, mdb=False, bname=False) :
    """
    Run the pairwise blastn process
    """

    if mdb :
        for fname in fnames :
            BCLine6.mbdb(fname)
    
    qlen = lambda fname : sum(len(record) for record in SeqIO.parse(fname, "fasta"))
    chopDic = {fname : chop(fname) for fname in fnames}
    qlens = {fname : qlen(fname) for fname in fnames}

    args = [utils.FunArgs(blastduo, blastn, duo, str(chopDic[duo[0]]), duo[1], qlens[duo[0]]) 
            for duo in permutations(fnames, 2)]
    
    results = utils.mproc(args, ncore)
    df = pd.DataFrame(results)

    if bname :
        df["F1"] = df["F1"].apply(os.path.basename)
        df["F2"] = df["F2"].apply(os.path.basename)

    if outfile :
        df.to_csv(outfile, sep="\t")
    else :
        print (df)

def run(fastas, blastn="blastn", outfile=None, ncore=1, mdb=False, bname=False) :
    fnames = glob.glob(fastas)
    if len(fnames) < 2 :
        raise IOError("Less than two fasta has been found")

    results = pairwise_blastn(fnames, blastn, outfile, ncore, mdb, bname)
    

if __name__ == "__main__" :
    fire.Fire(run)
import os, glob
from itertools import combinations
import pandas as pd
from utils import flat_fdata, iter_fdata

def run(fasta1, fasta2, ** fkwargs) :
    fd1 = iter_fdata([fasta1], ** fkwargs)
    fd2 = iter_fdata([fasta2], ** fkwargs)

    fd1 = list(flat_fdata(fd1))
    fd2 = list(flat_fdata(fd2))

    r1found = []
    r2found = []

    for r1 in fd1 :
        for r2 in fd2 :
            if r1.seq == r2.seq :
                r1found.append(r1.id)
                r2found.append(r2.id)
                print (f"Found match between {r1.id} and {r2.id}")
                break

    print (fasta1 + ":")
    print (f"{len(fd1)} records - {len(r1found)} found",
           f"(prc : {len(r1found) * 100 / len(fd1)})")

    print (fasta2 + ":")
    print (f"{len(fd2)} records - {len(r2found)} found",
           f"(prc : {len(r2found) * 100 / len(fd2)})")

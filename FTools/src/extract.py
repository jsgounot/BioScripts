# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2020-12-08 11:28:53
# @Last modified by:   jsgounot
# @Last Modified time: 2020-12-08 12:51:57

import argparse
from itertools import islice, chain
from random import sample
from Bio.Seq import Seq
from utils import iter_fdata

def show_fdata(fdata, revcomp=False, ncount=0, rand=False) :
    if ncount and not rand : fdata = islice(fdata, 0, ncount)
    if ncount and rand : fdata = sample(list(fdata), ncount)

    for record in fdata :
        rseq = record.seq.reverse_complement() if revcomp else record.seq
        rname = record.id

        print (">%s" %(rname))
        for i in range(0, len(rseq), 60) :
            print (rseq[i:i+60])

def subseq(record, start, end) :
     start = start if start != -1 else 0
     end = end if end != -1 else len(record.seq)
     record.seq = record.seq[start:end]
     return record

def qcor(record) :
     record.seq = Seq(str(record.seq).replace("?", "N"))
     return record

def run(files, minsize=0, maxsize=0, start=-1, end=-1, names=[], revcomp=False, ncount=0, rand=False, 
    bcorrect=False, sort=None, reverse=False) :

    fdata = chain(* (fdata for fname, fdata in iter_fdata(files)))

    if minsize :
        fdata = (record for record in fdata if len(record.seq) >= minsize)

    if maxsize :
        fdata = (record for record in fdata if len(record.seq) <= maxsize)            

    if names :
        fdata = (record for record in fdata if record.id in names)

    if bcorrect :
        fdata = (qcor(record) for record in fdata)

    if start > 0 or end > 0 :
        fdata = (subseq(record, start, end) for record in fdata)

    if sort == "size" :
        fdata = sorted(fdata, key = lambda record : len(record.seq), reverse=reverse)

    if sort == "name" :
        fdata = sorted(fdata, key = lambda record : record.id, reverse=reverse)            

    show_fdata(fdata, revcomp, ncount, rand)


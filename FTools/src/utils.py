# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2020-12-08 12:30:50
# @Last modified by:   jsgounot
# @Last Modified time: 2020-12-08 13:20:34

import glob, os
from itertools import chain
from Bio import SeqIO

class FastaToolsError(Exception) :
    pass

def iter_fdata(fnames, ** filter_kwargs) :
    
    if not fnames :
        raise FastaToolsError("Please provide at least one fasta file")
    
    for fname in fnames :
        
        if not os.path.isfile(fname) :
            raise FastaToolsError("File not found : " + fname)
        
        if not (fname.endswith(".fa") or fname.endswith(".fasta")) :
            if not filter_kwargs.get("force") :
                print ("WARNING : No proper fasta extension for : " + fname)
                print ("File ignored (use --force to still use it)")
                continue

        fdata = SeqIO.parse(fname, "fasta")
        yield fname, filter(fdata, ** filter_kwargs)

def flat_fdata(fdata) :
    return (record for fname, records in fdata
        for record in records)

def filter(fdata, ** filter_kwargs) :
    minsize = filter_kwargs.get("minsize", None)
    maxsize = filter_kwargs.get("maxsize", None)
    names   = filter_kwargs.get("names", None)

    for record in fdata :
        if minsize and len(record.seq) < minsize :
            continue
        
        if maxsize and len(record.seq) > maxsize :
            continue

        if names and record.id not in names :
            continue

        yield record


# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2020-12-08 12:30:50
# @Last modified by:   jsgounot
# @Last Modified time: 2021-02-17 14:31:08

import gzip
import glob, os
from itertools import chain

from Bio import SeqIO

from processing import FunArgs, mproc

class FastaToolsError(Exception) :
    pass

class FastaToolsException(Exception) :
    pass

def iter_fdata(fnames, ncore=1, post_fun=None, pkwargs={}, ** iter_kwargs) :
    allowed_exts = [".fa", ".fasta", "fna"]
    allowed_exts += [ext + ".gz" for ext in allowed_exts]

    if not fnames :
        raise FastaToolsError("Please provide at least one fasta file")
    
    fargs = [FunArgs(process_apply, fname, allowed_exts, 
            post_fun=post_fun, pkwargs=pkwargs, ** iter_kwargs)
            for fname in fnames]

    return (result for result in mproc(fargs, ncore) if result is not None)

def process_apply(* args, post_fun=None, pkwargs={}, ** kwargs) :
    fname, fdata = process_fasta(* args, ** kwargs)

    if fname is not None and post_fun is not None : 
        return post_fun(fname, fdata, ** pkwargs)
    
    return fname, fdata

def process_fasta(fname, allowed_exts, gzipped=False, ** filter_kwargs) :
    if not os.path.isfile(fname) :
        raise FastaToolsError("File not found : " + fname)
    
    if not any(fname.endswith(ext) for ext in allowed_exts) :
        if not filter_kwargs.get("force") :
            print ("WARNING : No proper fasta extension for : " + fname)
            print ("File ignored (use --force to still use it)")
            return None

    if fname.endswith(".gz") or gzipped :
        with gzip.open(fname, "rt") as handle :
            fdata = SeqIO.parse(handle, "fasta")
            fdata = list(filter(fdata, ** filter_kwargs))
            return fname, fdata
            
    else :
        fdata = SeqIO.parse(fname, "fasta")
        return fname, filter(fdata, ** filter_kwargs)

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


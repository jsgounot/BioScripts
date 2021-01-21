# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2020-12-08 11:27:43
# @Last modified by:   jsgounot
# @Last Modified time: 2021-01-21 17:31:28

import os, glob
import pandas as pd
pd.set_option('display.width', 1000)

bname = os.path.basename
dname = os.path.dirname

from utils import iter_fdata, FastaToolsError

def ace_info(fsizes, nvalue, refsize=0) :

    try : nvalue = int(nvalue)
    except ValueError :
        FastaToolsError("Nvalue incorrect : %s" %(str(nvalue)))

    data = {}

    totalN_length = sum(fsizes.values())
    prcN = totalN_length * (nvalue / 100.)
    prcNG = refsize * (nvalue / 100.)

    nfound = ngfound = False
    sum_size = 0
    for idx, name in enumerate(sorted(fsizes, key = lambda x : fsizes.get(x), reverse = True)) :
        sum_size += fsizes[name]
        
        if sum_size >= prcN and nfound == False :
            data["L%i" %(nvalue)] = idx + 1
            data["N%i" %(nvalue)] = fsizes[name]
            nfound = True

        if refsize and sum_size >= prcNG and ngfound == False :
            data["LG%i" %(nvalue)] = idx + 1
            data["NG%i" %(nvalue)] = fsizes[name]
            ngfound = True

    return data


def run(files, ref, refsize, nvalue=50, ** fkwargs) :
    data = []

    if refsize == 0 and ref :
        refsize = sum(len(record.seq) 
            for _, fdata in iter_fdata((ref,))
            for record in fdata)
    
    for fname, fdata in iter_fdata(files, ** fkwargs) :
        fsizes = {record.id : len(record.seq) for record in fdata}

        try : name = bname(fname.name)
        except AttributeError : name = bname(str(fname))

        finfo = {
            "basename" : name, 
            "averageSize" : int(sum(fsizes.values()) / len(fsizes)),
            "seqNumber": len(fsizes), 
            "maxSize"  : max(fsizes.values()), 
            "minSize"  : min(fsizes.values())
        }


        finfo = dict(** finfo, ** ace_info(fsizes, nvalue, refsize))
        data.append(finfo)

    df = pd.DataFrame(data)

    cn = lambda x, i : x + str(i)
    nlcols = [cn("N", nvalue), cn("L", nvalue), cn("NG", nvalue), cn("LG", nvalue)]
    nlcols = [col for col in nlcols if col in df.columns]

    if not df.empty :
        for column in ["minSize", "maxSize", "averageSize", "seqNumber"] + nlcols :
            df[column] = df[column].apply(lambda x : '{:,}'.format(x))
    
        df = df[["basename", "seqNumber", "minSize", "maxSize", "averageSize"] 
            + nlcols]
        
        with pd.option_context('display.max_rows', None, 'display.max_columns', None) :
            print(df)
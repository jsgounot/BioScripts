# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-01-24 13:37:16
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-01-24 13:49:12

import numpy as np

def tss_transform(df):
    sij = df.sum()
    return (df / sij) * 1000

def cumNormStat(matrix, rel=.1, default=.5):
    sa = np.sort(a, axis=0)
    ref = sa.mean(axis=1)
    refS = np.sort(ref)

    k = np.argmax(refS > 0)
    srefS = refS[k:]

    lo = (len(refS) - k)
    s = np.linspace(0, 1, lo)
    fun = lambda column: np.quantile(column[column != 0], q=s)

    diffr = np.apply_along_axis(fun, 0, a)
    diffr = np.reshape(srefS, (1, -1)).T - diffr

    diffr2 = np.median(np.absolute(diffr), axis=1)
    x = (np.argmax((np.absolute(np.diff(diffr2)) / diffr2[1:]) > rel) + 1) / len(diffr2)
    
    if x < default:
        x = default
        print ('Default value being used.')
    
    return x

def css_transform(df):
    a = df.values
    x = cumNormStat(a)

    q = df[df != 0].quantile(q=x)
    sij = df[df <= q].sum()
    return (df / sij) * 1000
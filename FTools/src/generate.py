# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2020-12-08 14:25:26
# @Last modified by:   jsgounot
# @Last Modified time: 2020-12-08 15:20:34

import random

from utils import FastaToolsError

def run(basename, start, number, size, uppersize, bases, weights, seed) :
    if seed : random.seed(seed)

    weights = weights.split("/")
    weights = [float(weight) for weight in weights]

    if len(weights) != len(bases) :
        raise FastaToolsError("""The number of weights must be equals
            to the number of bases""")

    for i in range(number) :
        i += start
        name = basename + str(i)

        try : random.choices
        except AttributeError : 
            raise FastaToolsError("""Unable to load random choice, please upgrade your
                python version to be >= 3.6""")

        seqSize = random.randint(size, uppersize) if uppersize else size 
        seq = random.choices(bases, weights=weights, k=seqSize)
        seq = "".join(seq)

        print (">" + name)
        for i in range(0, len(seq), 60) :
            print (seq[i:i+60])
# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-03-12 14:47:50
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-03-12 15:00:40

import gzip
from dendropy import PhylogeneticDistanceMatrix as PDM

fname   = snakemake.input[0]
outfile = snakemake.output[0]

# Load
with gzip.open(fname, "rt") as f :
    pdm = PDM.from_csv(f)

# Create tree
tree = pdm.nj_tree()

try :
    tree.reroot_at_midpoint(update_bipartitions=False, suppress_unifurcations=False)
except AssertionError:
    print ("Unable to root the tree, ignore")   


# Save tree
with open(outfile, "w") as f :
    f.write(tree.as_string("newick"))
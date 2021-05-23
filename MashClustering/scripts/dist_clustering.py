# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-05-23 16:52:04
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-05-23 17:32:54

import glob, os

from itertools import combinations

from sklearn.cluster import AgglomerativeClustering as AC
import networkx as nx

import pandas as pd
import numpy as np

fname = snakemake.input[0]
out_matrix = snakemake.output.matrix
out_clusters = snakemake.output.clusters
clustering_threshold = snakemake.params.threshold

names = ["reference-ID", "query-ID", "distance", "p-value", "shared-hashes"]

print ("read ", fname)
df = pd.read_csv(fname, sep="\t", names=names, compression="gzip")
df["reference-ID"] = df["reference-ID"].apply(lambda path : os.path.basename(path))
df["query-ID"] = df["query-ID"].apply(lambda path : os.path.basename(path))

print ("make matrix ...")
matrix = pd.pivot_table(df, index="reference-ID", columns="query-ID", values="distance", fill_value=1)

print ("launch clustering ...")
clustering = AC(linkage="single", n_clusters=None, compute_full_tree=True, 
                distance_threshold=clustering_threshold, affinity="precomputed").fit(matrix)

identifiers = sorted(df["reference-ID"].unique())
data = []

for idx, identifier in enumerate(identifiers) :
    clusterID = clustering.labels_[idx]
    data.append({"sequence" : identifier, "clusterID" : clusterID})
                            
cdf = pd.DataFrame(data)

print ("%i clusters found" %(cdf['clusterID'].nunique()))
print ("Clusters size stats :")
print (cdf.groupby("clusterID").size().describe())

# Looking for central node
# based on mash distance

print ("search for center nodes ...")
dist = matrix.to_dict()

all_nodes_order = {}
all_eigen_values = {}

for cid, sdf in cdf.groupby("clusterID") :
    bins = set(sdf["sequence"])

    if len(bins) > 1 :
        graph = nx.Graph()   
        for s1, s2 in combinations(bins, 2) :
            mdist = 1 - dist[s1][s2]
            graph.add_edge(s1, s2, weight=mdist)

        eigen_values = nx.eigenvector_centrality(graph, weight="weight") # called cnode before   
        nodes_order = sorted(eigen_values, key=eigen_values.get)[::-1]
        nodes_order = {bin : idx for idx, bin in enumerate(nodes_order)}
        # cnode = max(cnode, key=cnode.get)

    else :
        bin = list(bins)[0]
        eigen_values = {bin : np.nan}
        nodes_order = {bin : 0}

    assert len(set(all_nodes_order) & set(nodes_order)) == 0
    all_nodes_order.update(nodes_order)
    all_eigen_values.update(eigen_values)

cdf["eigen_index"] = cdf["sequence"].apply(lambda name : all_nodes_order[name])
cdf["eigen_weight"] = cdf["sequence"].apply(lambda name : all_eigen_values[name])

matrix.to_csv(out_matrix, sep="\t")
cdf.to_csv(out_clusters, sep="\t")
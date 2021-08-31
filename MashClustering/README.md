# Mash Clustering

MAGs derepliction based on [Mash](https://github.com/marbl/Mash) distances.

1. Sketch sequences
2. Calculate mash distances 
3. Cluster MAGs using an [agglomerative clustering](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.AgglomerativeClustering.html) (sklearn)
4. Define central nodes based on [eigen centrality](https://networkx.org/documentation/networkx-1.10/reference/generated/networkx.algorithms.centrality.eigenvector_centrality.html) (networkx).

# Setup

```
mamba create -n mash_clustering -c bioconda -c anaconda -c conda-forge mash networkx scikit-learn pandas snakemake
```

You will need to replace the files paths directly inside the SnakeFile.

# See also

* [Galah](https://github.com/wwood/galah) : Another dereplication method using ANI distances and define central nodes based on checkM results. The Mash pipeline remains must faster.
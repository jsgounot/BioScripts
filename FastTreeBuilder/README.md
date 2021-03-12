# FastTreeBuilder

snakemake pipeline to produce UPGMA and NJ trees from group of fasta file.

1. Runs [MASH](https://github.com/marbl/Mash) distance between fasta files and produces matrix
2. Produce rooted trees with [dendropy](https://github.com/jeetsukumaran/DendroPy)

# Dependancies

Python

* pandas
* dendropy

Softwares

* snakemake
* mash

`conda create --name ftb -c bioconda mash dendropy pandas snakemake`

# Notes

A couple of mash options, including mash exectuable path, can be modified in the `config.json` file. Final files can be found in the `tree` directory. This pipeline has been designed to consider each subdirectory in `/fasta` as a sample, and each fasta files inside this subdirectory as a final node of the tree.

Inspired by [this github directory](https://github.com/afelten-Anses/QuickPhylo).
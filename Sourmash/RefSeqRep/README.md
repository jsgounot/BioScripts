A snakemake pipeline to build sourmash database for RefSeq representative genomes. 

##### Why not sourmash db?

You can download precomputed databases with all genbank genomes for each collection directly on the [sourmash website](https://sourmash.readthedocs.io/en/latest/databases.html#genbank-bacterial). Those contain all genomes available, which is much more complete but also make the search much slower. Additionally, the sourmash database need to be indexed, which need time and ressources.

##### How

I recommend to use a lot more CPUs than what you have with a greedy scheduler. You can try with:

```bash
snakemake -s sourmash_index.snk -c 100 --use-conda --scheduler greedy --until protozoa/protozoa.sbt.zip
```

You need `sourmash` (saved here in the conda environment named with the same name) and pandas in the snakemake environment.

Does not work with viral, there are no representative genomes for this collection. You can use the sourmash database for this, it's quite small.
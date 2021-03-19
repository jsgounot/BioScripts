# Kraken2 / Bracken pipeline

Pipeline to quickly run [kraken2](http://ccb.jhu.edu/software/kraken2/) and [bracken](https://ccb.jhu.edu/software/bracken/) for several samples and database. Some arguments (`paired`, `gzip-compressed` and `bzip2-compressed`) are infered based on your input data. You can custom all command lines with config file parameters (see configuration file section).

# Dependancies

Python

* snakemake

Softwares

* kraken2
* bracken

`conda create -n kraken2 kraken2 bracken snakemake`

# Configuration file

Pipeline use a configuration file (see example in this directory). By default, only reads and database are mandatory. You can customize kraken2 and bracken calls by adding options to a each experiments with keys `kraken_adding` and `bracken_adding` (see example config). `reads`, `kraken_adding` and `bracken_adding` are expected to be list. `bracken_ranks` defines bracken calls (`-l` option) : `{"species" : "S", "genus" : "G"}` means that 2 output files will be produced : `results/sample/name.bracken.species.report.tsv` (`-l S`) and `results/sample/name.bracken.genus.report.tsv` (`-l G`).

By default, the number of threads (4), the memory (1000MB) for each job and the used bracken ranks are defined in the header of the `Snakefile`. You can also defined the number of threads and memory for each experiment using json keys `threads` and `memory` (see example config file).
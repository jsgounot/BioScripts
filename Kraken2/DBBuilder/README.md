
# Kraken2 DBBUILDER

A simple script to help to quickly produce a [kraken](http://ccb.jhu.edu/software/kraken2/) database from scratch. This script is basicaly following recommandation from the [documentation](https://github.com/DerrickWood/kraken2/wiki/Manual#custom-databases) with basic options. Please refere to this document for additionnal parameters. Additionnaly, the script will also create in the same directory the [bracken](https://ccb.jhu.edu/software/bracken/) database.

# Dependancies

Python

* click
* biopython

Softwares

* kraken2
* bracken

Tested on :

* kraken2 v.2.1.1
* bracken v.2.6.0

`conda create -n kraken -c bioconda kraken2 bracken click biopython`

# How to

```
python make_db.py --help
Usage: make_db.py [OPTIONS] FASTAS OUTDIR

Options:
  --rnode INTEGER          mother taxid node
  --taxidsfile TEXT        fname to taxid map, see readme
  --ignore_missing_taxids  ignore missing taxids
  --taxonomy TEXT          taxonomy directory
  --krakenb TEXT           kraken2-build path
  --krakeni TEXT           kraken2-inspect path
  --brackenb TEXT          bracken-build path
  --startid INTEGER        The starting value for your taxonomic ids
  --threads INTEGER
  --help                   Show this message and exit.
```

Most of the options should not be used if you're using a conda install. You should try with a subset of you sequences first (< 10 fasta)

## taxidsfile

In case you have genomes from the same species or cluster you want to group together under the same node, you need to provide a simple csv file (without header) with two columns, first one being the fasta path (as provided in the software input), the second being the clusterID (an integer). Note that this integer will be added to the higher taxid from the current database to avoid duplicated taxid in the database. Node will be named after this taxid and the first MAG found under this taxid, example `taxid_2819974:your.bin.name.fa`, which will be reported in the kraken result. Option `ignore_missing_taxids` (default value = `True`) defines if the script should raise an Error if one genome is not found in the `taxidsfile`.

## Memory usage

Please be aware that kraken2 database building may need a lot of memory, especially if you use thousand of sequences.

## names and nodes files

You will need an existing taxonomic directory to produce your database and more precisely the `names.dmp` and `nodes.dmp` files inside it. You can download such directory with this command line `kraken2-build --download-taxonomy --db $DBNAME` (such as `bacteria`). However you can also unzip the names and nodes files I provide here and the script will automaticaly use them if no taxonomic directory is provided (downloaded 18 March 2021).

## Rnode parameter

The `rnode` parameter is the taxid node where your sequence is going to refere to. You can use a value of `2` if you're using the provided nodes and names files (corresponding to bacteria) but you might need to change this value if you're using another taxonomy database.

## What is going to reported in your report

The basename of your files, if multiple files are used for one single node, one of these files will be pick (should be the first one, but might be random).
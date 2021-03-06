# MapCallNGMLR

Basic snakemake mapping and calling process for long reads using [NGMLR](https://github.com/philres/ngmlr). This pipeline manages on its own reference index. You can specify one or multiple set of reads for one sample. Similar to [MapCallGATK](https://github.com/jsgounot/BioScripts/tree/main/MapCallGATK).

**Note** : Please be sure to add proper right to scripts/mapped.stat.sh before running this pipeline !

# Dependancies

Python

* snakemake
* pandas

Softwares

* ngmlr
* gatk
* picard
* samtools
* rtg (tools)

All softwares listed above must be executable with their name as indicated and can be downloaded as such using conda. Note that this pipeline is done for GATK 4.

Tested on :
* bwa  0.2.7
* gatk v4.1.9.0
* picard 2.23.3
* samtools 1.7 (using htslib 1.7)
* rtg tools 3.11

# Pipeline workflow

See [MapCallGATK](https://github.com/jsgounot/BioScripts/tree/main/MapCallGATK) for an overview of the workflow (only `bwa mem` is replaced).

# Configuration file

This pipeline uses a json configuration file. For each sample, you have to specify at least its reference and at least one group of reads. Note that if you only have one read or paired-read, you still have to put them between two brackets (see config file).

You can specify either one ploidy (gatk variant calling) by specifying an integer or multiple ploidies using a list.

# Additionnal note on stats

This pipeline produces some stats :
* Prc of mapped reads
* Windows coverage depth
* VCF stats

See [this page](https://github.com/jsgounot/BioScripts/tree/main/MappingStat) for more informations and see how to configure the first two scripts. 

For VCF statistic, the first file is produced using rtg tool `vcfstats` command. This file is further extended using a python script to show the number of variants per kb.
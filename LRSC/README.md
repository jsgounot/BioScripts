# LRCS

A python script to determine (Long Reads) connectivity with SNPs. While this script has been designed for long reads, it can be used for any couple of bam / vcf files.

## Context

Some phasing tools such as [gretel](https://github.com/SamStudio8/gretel) or [nPhase](https://github.com/OmarOakheart/nPhase) rely on shared polymorphic sites between reads to either determine a path along the reference or to fuse reads together. This tools provides a way to get all best scenario paths, meaning all windows of your genome which are covered by reads sharing at least one polymorphic site. While the results does not give you information on the actual results of these tools, it will let you know if your input file is suitable for such analyses.

## Dependancies

* click
* pysam
* pandas
* tqdm (optional)

## Input

* bamfile
* vcffile

Note that your vcf file might not originated from your bamfile, for instance if you want to use short read SNPs for long reads clustering.

## Basic command :

`python lrcs.py --help`

`python lrcs.py --bam bamfile --vcf vcffile`

It is recommanded to sort your bamfile before this script. By default, lrcs will consider it is not, leading to a supplementary merging step. If you bamfile is sorted, you can avoid this step by using the `--sortedbam` flag.

You can also consider chimeric reads as a single read (as used by `nphase` for example) with the `--chimeric` option. Note that enabling chimeric flag will necessarily lead to a merging step similar to `--sortedbam`.
# Mapping stat

Snakefile to get basic mapping stats from sam file
Can be easily merged to another pipeline.

# Process

* Sort reads (samtools)
* Define prc mapped and unmapped reads
* Extract mapped reads (samtools)
* Calculate depth (samtools)
* Define average coverage for window size N (default 10000)

# Dependancies

* snakemake
* python pandas

# Additionnal note

Window size can be modified directly into the scripts/window_depth.py file : `WINDOW`. Additionnaly, you can modify `QUANTILE` corresponding to the upper quantile of the mean coverage value, used as the maximum coverage value to remove extreme values.

# Mapping stat

Snakemake pipeline to get basic mapping stats from sam file
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
* samtools

# Additionnal note

Be sure to give proper right to shell script on scripts/mapped.stat.sh (chmod +x) !

Window size can be modified directly into the scripts/window_depth.py file : `WINDOW`. Additionnaly, you can modify `QUANTILE` corresponding to the upper quantile of the mean coverage value, used as the maximum coverage value to remove extreme values. Note that this pipeline does not remove intermediate files.
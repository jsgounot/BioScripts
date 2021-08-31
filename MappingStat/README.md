# Mapping stat

Snakemake pipeline to map reads and get window coverage.
Can be easily merged to another pipeline.

# Dependancies

* snakemake
* python pandas
* minimap2
* samtools

# Additionnal note

Be sure to give proper right to shell script on scripts/mapped.stat.sh (chmod +x) !
You will obtain 2 window files for each sample. The first one is a 50bp window and the second one is 1kb window.
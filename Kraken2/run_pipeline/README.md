# Kraken2 / Bracken pipeline

Snakemake pipeline for [kraken2](http://ccb.jhu.edu/software/kraken2/) and [bracken](https://ccb.jhu.edu/software/bracken/) for several samples and database. Some arguments (`paired`, `gzip-compressed` and `bzip2-compressed`) are infered based on your input data. You can custom all command lines with config file parameters (see configuration file section). Handle optional spring files input (Illumina paired-end only).

```bash
mamba create -n kraken2 kraken2 bracken snakemake
# Make your configuration file
snakemake -s kraken2.snk -c {cores} --configfile {config.path} (--use-conda) (--keep-going)
```

# Configuration file

Pipeline use a configuration file (see example in this directory). By default, only reads and database are mandatory. You can customize kraken2 and bracken calls by adding options to a each experiments with keys `kraken_adding` and `bracken_adding` (see example config). `reads`, `kraken_adding` and `bracken_adding` are expected to be list. `bracken_ranks` defines bracken calls (`-l` option) : `{"species" : "S", "genus" : "G"}` means that 2 output files will be produced : `results/sample/name.bracken.species.report.tsv` (`-l S`) and `results/sample/name.bracken.genus.report.tsv` (`-l G`).

By default, the number of threads (4), the memory (1000MB) for each job and the used bracken ranks are defined in the header of the `Snakefile`. You can also defined the number of threads and memory for each experiment using json keys `threads` and `memory` (see example config file).

**Basic json file**

```json
{ 
	"sample_name" : {
		"spring": "path/to/reads.spring",
		"db1" : {
			"database" : "path/to/database_1",
			"bracken_ranks" : {"species" : "S"},
		}
	}
}
```

**A bit more complexe**

```json
{ 
	"sample_name" : {
		"reads": [
			"path/to/read_1.fastq.gz",
			"path/to/read_2.fastq.gz"
		],
		"db2" : {
			"database" : "path/to/database_2",
			"memory" : 1200,
			"threads" : 8,
			"kraken_adding" : ["--minimum-hit-groups 2"],
			"bracken_adding" : ["-t 0"],
			"bracken_ranks" : {"species" : "S"},
			"paired" : "False"
		}
	}
}
```

# Utility script(s)

## `extract_reads.py`

Extract the reads associated to a taxonomic ID, based on Kraken2 output files.

```bash
python extract_reads.py --help
usage: extract_reads.py [-h] [--paired] [--strict] [--count [COUNT]] [--seed [SEED]] [--quiet] output taxid [fastq ...]

positional arguments:
  output           Kraken output file
  taxid            Target taxonomic ID. Unclassified is 0
  fastq            Fastq file(s) (default: None)

options:
  -h, --help       show this help message and exit
  --paired         If paired-end fastq, return reads only if both pairs are assigned to taxid (default: False)
  --strict         Raise an error if one read is not found (default: False)
  --count [COUNT]  Randomly pick and output COUNT reads. 0 = all reads (default: 0)
  --seed [SEED]    Seed number when using count (default: 1)
  --quiet          Might reduce RAM usage (default: False)
```

Example:

```bash
python extract_reads.py kraken2/groupname/reports/sample/database.kraken.output.tsv.gz 7062 sample.fastq.*.fq.gz --count 10 --strict | gzip > sample.7062.fastq.gz
```

Output is a fastq format, that you can easily convert to fasta on the fly

```bash
python extract_reads.py kraken2/groupname/reports/sample/database.kraken.output.tsv.gz 7062 sample.fastq.*.fq.gz --count 10 --strict | sed -n '1~4s/^@/>/p;2~4p' | gzip > sample.7062.fa.gz
```


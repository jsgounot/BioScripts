# Hocort pipeline

A [hocort](https://github.com/ignasrum/hocort) / [spring](https://github.com/shubhamchandak94/Spring) / [fastp](https://github.com/OpenGene/fastp) decontamination snakemake pipeline for paired-short reads (for now). The pipeline will: (1) uncompress your reads, (2) run fastp and hocort, (3) recompress reads and (4) provides statistic on the decontamination.

### Setup

A. Prepare a conda environment with spring, hocort, fastp and (eventually) snakemake.

```bash
mamba create -n hocort_pipe -c conda-forge -c bioconda spring hocort fastp snakemake
```

B. Index your database using hocort (check project's readme).

C. Edit the snakemake file header with your conda environment (if the conda environment is not the same than the snakemake one) and database path.

### Usage

The configuration file should be formatted like this:

```json
{
	"sample1": {"spring": "path/to/springfile"},
	"sample2": {"spring": "path/to/springfile"}
}
```

Or if you have fastq files

```json
{
	"sample1": {
		"r1": "path/to/fastq.r1.fq.gz",
		"r2": "path/to/fastq.r2.fq.gz"
	}
}
```

Run the pipeline:

```bash
snakemake -s hocort_spring.snk --configfile config.json -d hs_wdir --cores 4 -np
```

You can use `--use-conda` if you want to use a different conda environment.

Once over, you can have a quick view of the results using the `summary_stat.py` script:

```
python summary_stat.py hs_wdir/data/fastp_stats/*
```

with three important columns: `%fastp_base` = % reads after fastp, `%decont_bases` = % reads after fastp and decont (compared to raw again), and `decont_diff` = `%fastp_base` - `%decont_bases` = % reads lost from human decontamination only.

### Extension

Curent pipeline only uses bowtie2, but Hocort includes more classifers / mappers. You can add another one by duplicating the `hocort_bowties` rule.
# Reads compression using Spring

Snakemake pipeline to compress reads for storage purpose using [Spring](https://github.com/shubhamchandak94/Spring). This pipeline has been designed to compress **paired-end** illumina reads. To compress single end short-read, you will need to modify the pipeline. Note that [a specific version of Spring](https://github.com/qm2/NanoSpring) has been designed for nanopore reads.

Based on my results, the spring file is around 3 to 4 times smaller than the fastq.gz file. The major con is the time for compression and decompression, but both operation can be threaded.

### Usage

The pipeline will can either work on local `fastq` files or download the files for you with `wget`. In case of local `fastq`, files are **not** removed from their original location. You will need to remove the files yourself.

You need to feed the pipeline with a json configuration file:

```json
{
    "ERR7671874": {
        "r1": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR767/004/ERR7671874/ERR7671874_1.fastq.gz",
        "r2": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR767/004/ERR7671874/ERR7671874_2.fastq.gz"
    },
    "ERR7671875": {
        "r1": "/local/path/ERR7671875_1.fastq.gz",
        "r2": "/local/path/ERR7671875_2.fastq.gz"
    }
}
```

 The pipeline can be launched (dry-run) this way:

```bash
snakemake -s spring.snk --configfile spring_config.json -d spring --cores 4 -np
```

### Workflow

The pipeline will not only compress your fastq file but also check that the compressed file returns the initial file. Workflow for each sample can be summarized this way:

* Compression
* Decompression
* Check initial `fastq` files with decompressed `fastq` files (`zcmp`)
* If a diff is observed, remove the third line in all `fastq` files with `awk` (see note) and run the same comparison

All these operations are run with `gziped` `fastq` file. You can also skip the check phase running the pipeline with `--until compress`.

### About the third line

Spring does not conserve the optional part in the third line of the `fastq` file ([fastq format](https://en.wikipedia.org/wiki/FASTQ_format)). This should not be an issue for your analysis since this part is optional but to confirm that only the third line is an issue, an additional `zcmp` is done with `fastq` files without this third line.

### Outputs

You should have both the `spring` files, `cmp` files should be empty. You can double check that all files were checked (no differences observed between the initial fastq and decompressed fasta files) by checking `output/compression.check.txt`, all samples should have  a `True` value in their 3 column that you can check with this command line: `awk '{ print $3}' output/compression.check.txt | uniq -c`. Finally, a table `output/compression.stat.tsv` contains information about the compression rate for each sample.
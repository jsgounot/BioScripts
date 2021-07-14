# Reads compression using Spring

Snakemake pipeline to compress reads for storage purpose using [Spring](https://github.com/shubhamchandak94/Spring).
This pipeline has been designed for pre-compressed gzip **Illumina** fastq reads. See Spring option (such as `-l`) for long reads.
Based on my results, the spring file is around 3 time smaller than the fastq.gz file, both compression and decompression can be threaded.

# Safety check and usage

Since the pipeline will compressed all `*.fastq.gz` files found in `./input`. Initial fastq.gz files are conserved and will not be removed.

Once fastq files are compressed, the pipeline will check whether the same initial input can be retrieved using the compressed file. To do that, the spring file is uncompressed and both the new and the original `fastq.gz` files are compared using `zcmp` If both file are identical, the compression is considered successfull. 

If not, this can arise because of variation in the optional third line, which does not have an impact on further analyses. Therefore, an additionnal step is processed where the third line of each read is removed in both files using `awk` and reads are compared again. This last step is especially intensive in term of disk usage and I would recommand to have at least the total size of initial reads files as free space on your working disk. Temporary files should be removed at the end of the process.

At the end of the process, a `.cmp` file is generated for each reads file. An empty `cmp` file means that no difference has been observed between the original and the uncompressed file. 

It is really not recommanded, but you can also ignore the safety check using this command line `snakemake --cores 16 --until compress`.
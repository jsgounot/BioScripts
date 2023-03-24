## Bam scripts collection

Personal bam scripts collection

### compare_bam.py

Compare two bam files and extract unique reads ID from each of them. Allows statistic and reads grep.

```
$ python compare_bam.py --help
usage: compare_bam.py [-h] {bamcompare,readstat,readextract} ...

positional arguments:
  {bamcompare,readstat,readextract}
                        sub-command help
    bamcompare          Compare bam and extract reads IDS
    readstat            Extract statistic based on bamparser IDs
    readextract         Extract reads based on bamparser IDs

optional arguments:
  -h, --help            show this help message and exit
```

**Important note:** `bamcompare` uses some filter to select rids, check `python compare_bam.py bamcompare --help`

Reads similarity between two bams (`head` = only 10k)

```bash
python compare_bam.py bamcompare --bam1 file1.bam --bam2 file2.bam --head | python compare_bam.py readstat --rids -
```

Extract reads found in both bam (with grep: small number of rids)

```bash
python compare_bam.py bamcompare --bam1 file1.bam --bam2 file2.bam | python compare_bam.py readscompare --kind intersection --rids - | zgrep -A 3 -f - /path/to/fastq.gz | gzip > test.fastq.gz
```

Extract reads found in bam1 only

```bash
python compare_bam.py bamcompare --bam1 file1.bam --bam2 file2.bam | python compare_bam.py readscompare --kind bam1 --rids - | zgrep -A 3 -f - /path/to/fastq.gz | gzip > test.fastq.gz
```

Because grep can be slow with a lot of rids, I provide a small utility to mimic grep `fastq.gz` as well:

```bash
python compare_bam.py bamcompare --bam1 file1.bam --bam2 file2.bam | python compare_bam.py readscompare --kind intersection --rids - | python compare_bam.py fastaextract --rids - --fastagz /path/to/fastq.gz
```

Because extracting rids can be slow and you might want to do multiple operation after, one way is to save the result first:

```bash
python compare_bam.py bamcompare --bam1 file1.bam --bam2 file2.bam > bamcompare.rids.json
cat bamcompare.rids.json | python compare_bam.py readscompare --kind intersection --rids - | python compare_bam.py fastaextract --rids - --fastagz /path/to/fastq.gz
```


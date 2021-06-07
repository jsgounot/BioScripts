# GTDBTk analyses

A collections of scripts used to work with [GTDBTk](https://github.com/Ecogenomics/GTDBTk) results.

# Dependancies

Python

* click
* biopython

# Example

```
python gtdbtk.consolidate.py data/gtdbtk.ar122.summary.tsv data/gtdbtk.bac120.summary.tsv pplacer_taxonomy gtdbtk.pplacer.ranks.tsv
python gtdbtk.sankey.py gtdbtk.pplacer.ranks.tsv
```


# Classic mash distances run

Snakemake pipeline to run [mash](https://github.com/marbl/Mash) between two set of sketches and extract the best (smallest) distance for each query.
By default, a maximum mash distance of 0.3 is used, you can change this directly in the snakefile.

# Dependancies 

* Mash
* Python pandas

# Sketching

Before running `mash distance`, you need to sketch your fasta sequences. You can easily do it like this:

```
ls -d /path/to/your/fasta/files/*.fa > sketchlist.txt
mash sketch -l sketchlist.txt -p 10 -k 21 -s 10000 -o yoursketchfile.msh
```

With `k` being k-mer size, `s` the sketch size and `p` the number of threads used. 
Note that there is a *blur* limit of the total size of sequences you can sketch (depending also of the sketch size). This is why it's actually better to split big database in multiple sketches.

# How it works

You will need to change the path of the query and reference sketches directly inside the snakefile script.
For both queries and references sketches set, a json file is produced where is stored an unique identifier as an integer.
Mash distance results between one query and a subject will be stored as `workingdir/queryid_targetid.tsv.gz`.
The json is saved in the working dir and load / updated at each execution, meaning that you can update your sketch list with new sequences and previously generated pair results will not be regenerated.
However, you might need to remove the previously generated `mash_results.tsv` file.

# Outputs

The best mash result for each queries is given in the `mash_results.tsv` file. All mash distances are available in `workingdir/*.dist.tsv.gz` files. 

Note that for the `mash_results.tsv` file, some queries can have duplicates if the query minimum mash distance is shared between two or more reference sequences.


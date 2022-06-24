
# Kraken2 DBBUILDER

A simple script to help to quickly produce a [kraken](http://ccb.jhu.edu/software/kraken2/) database from scratch. This script is basicaly following recommandation from the [documentation](https://github.com/DerrickWood/kraken2/wiki/Manual#custom-databases) with basic options and will also create in the same directory the [bracken](https://ccb.jhu.edu/software/bracken/) database. 

#### Before going further

Kraken classic database use NCBI taxonomic classifications with 2 files: `nodes.dmp` and `names.dmp` (see [here](https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_readme.txt) for official files description, files can be truncated). This means you need to associate your genomes to either an existing NCBI taxonomic ID or create a new one. All taxonomic ID and relationship between them are saved in the `nodes.dmp` file, while the `names.dmp` files store each node metadata. 

1. If your goal is just to include a couple of missing species / strain into a pre-existing database, please look at the kraken2 `--add-to-library` option. You will need to add an existing or novel taxid of your species into the fasta header too. 
2. If your goal is to create a new database from scratch based on a MAGs dataset, use this script.

### Dependancies

Python

* click
* biopython

Tested on :

* kraken2 v.2.1.1
* bracken v.2.6.0

`conda create -n kraken -c bioconda kraken2 bracken click biopython`

### How to

```
$ python make_db.py --help
Usage: make_db.py [OPTIONS] FASTAS OUTDIR

Options:
  --krakenb TEXT         kraken2-build path
  --krakeni TEXT         kraken2-inspect path
  --brackenb TEXT        bracken-build path
  -r, --gtdbtk_res TEXT  GTDBtk result files
  --nodes TEXT           GTDB archeal metadata file
  --names TEXT           GTDB bacterial metadata file
  --ext TEXT             Fasta files extension (.fa, .fasta), usefull if you
                         want to link to DRep results
  --drep_cdb TEXT        DRep CDB result for novel SPECIES
  --no-prune             Don't prune the tree to keep only used nodes
  --threads INTEGER
  --help                 Show this message and exit.
```

#### Command lines an input files

Depending of your input files and what you want, there are different ways to use this tool

##### You don't care about the taxonomic tree

`python make_db.py '*.fasta' outdir`

This will create a completly fake tree where each node (species,genus,...,domain) will be tagged with a unique identifier (novel*i*) and with each sequence being completly independant of other sequences. This is <u>not</u> recommanded but can be interesting if you only want to check the classification rate though.

##### You only have dreplication data (with drep)

You can provide your [DRep](https://github.com/MrOlm/drep) information (CBD file). This will only group sequences from the same cluster into the same species and will not provide meaningful information for parent nodes (genus and higher). 

`python make_db.py '*.fasta' outdir --drep-cdb Cdb.csv`

##### You only have GTDBTk results

GTDBTk classification can be used to guide the generation of the tree. This will generate a reticulated full tree with meaningful taxon name. The taxonomic ids will be completly random though.

`python make_db.py '*.fasta' outdir -r gtdbtk.ar122.summary.tsv -r gtdbtk.bac120.summary.tsv `

You will most likely need to use the `--ext` option as well to match your files basename to GTDBtk identifiers.

**Note**: In case of a novel nodes (for example `s__` for a novel species), the script will create a unique identifiers like `novel1`, `novel2`, etc. This can lead to potential issues in the taxonomy and can be fixed if you did a dereplication analysis of you sample. Let's say you have 3 novel genomes with 2 you know come from the same species and 1 which is the only representative of this novel species. All are coming from the same genus. By default the script is unable to know if genomes come from the same or different species, therefore 3 novel species will be created. To fix this, either you modify GTDBTk file to reflect the change using a unique species identifier to add to the taxonomic rank, such as `s__` becoming `s__novel_species1` , or you can use the `--drep-cbd` option.

##### You have GTDBTk results that you want to map over an existing taxonomic tree

You might want to use an existing taxonomic files (`names.dmp` and `nodes.dmp`) from another database **which was also generated with GTDBTk results**, like the [UHGG database](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/uhgg_kraken2-db/taxonomy/), to map your GTDBTk results. It means that the `names.dmp` files must contain names like `g__Turicibacter`. This will create a tree with similar taxon id than the one provided in the taxonomic files.

`python make_db.py '*.fasta' outdir -r gtdbtk.ar122.summary.tsv -r gtdbtk.bac120.summary.tsv --nodes nodes.dmp --names names.dmp`

By default, only nodes which are found in your sequences will be conserved, but you can keep them all with `no-prune`. This is not recommended though.

##### You can combine everything

`python make_db.py '*.fasta' outdir -r gtdbtk.ar122.summary.tsv -r gtdbtk.bac120.summary.tsv --nodes nodes.dmp --names names.dmp --drep-cdb Cdb.csv`

This will take the GTDBTk results, map them against the known tree and will group novel genomes from the same cluster into the same taxonomic node. This is by far the cleanest solution you can provide. 

#### Memory usage

Please be aware that kraken2 database building may need a lot of memory, especially if you use thousand of sequences.
# ANI

A simple and fast python script for blastn pairwise average nucleotide identity (ANI) calculation. Use [pyblast](https://github.com/jsgounot/PyBlast) and [google fire](https://github.com/google/python-fire).

## Basic command :

`python ani.py 'path/to/your/*.fasta' --ncore 4 --bname`

Otherwise :

`python ani.py --help`

## Notes 

Related publications or tools :
* [The original paper](https://pubmed.ncbi.nlm.nih.gov/17220447/) covering the ANI calculation (blastall)
* [jspeciesWS](http://jspecies.ribohost.com/jspeciesws/) : A graphical interface (blast+)
* [pyani](https://github.com/widdowquinn/pyani) : Another python library (blast+)
* [chjp/ani](https://github.com/chjp/ANI) : Similar but in perl (blastall)

The original calculation of ANI is based on blastall command. This python script uses blast+ with the following arguments `-evalue 1e-15 -max_target_seqs 1 -xdrop_gap_final 150 -dust no`. Results might differ from other tools depending of the blast version. 
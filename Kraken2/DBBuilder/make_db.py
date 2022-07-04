# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-11-14 17:56:19
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-07-04 10:27:07

import glob, os, shutil
from Bio import SeqIO
from multiprocessing import Pool
from collections import defaultdict

import gtdbtk2ncbi
import click

@click.command()
@click.argument('fastas', type=str)
@click.argument('outdir', type=str)
@click.option('--krakenb', default="kraken2-build", type=str, help="kraken2-build path")
@click.option('--krakeni', default="kraken2-inspect", type=str, help="kraken2-inspect path")
@click.option('--brackenb', default="bracken-build", type=str, help="bracken-build path")
@click.option('--gtdbtk_res', '-r', type=str, multiple=True, help="GTDBtk result files")
@click.option('--nodes', type=str, default='', help="GTDB archeal metadata file")
@click.option('--names', type=str, default='', help="GTDB bacterial metadata file")
@click.option('--ext', type=str, default='', help="Fasta files extension (.fa, .fasta), usefull if you want to link to gtdbtk files")
@click.option('--drep_cdb', type=str, default='', help="DRep CBD result for novel SPECIES")
@click.option('--add-strain-level', is_flag=True, help='Add a taxonomic ID for each genome as individual strain')
@click.option('--no-prune', is_flag=True, help='Don\'t prune the tree to keep only used nodes')
@click.option('--threads', default=1, type=int)
def run(fastas, outdir, krakenb, krakeni, brackenb, gtdbtk_res, nodes, names, ext, drep_cdb, add_strain_level, no_prune, threads) :
    if os.path.isdir(outdir) :
        raise Exception("Directory already exists : %s" %(outdir))

    fnames = glob.glob(fastas)
    print ("%i files to use" %(len(fnames)))
    assert fnames

    bnames = {os.path.basename(fname) for fname in fnames}
    if len(bnames) != len(fnames):
        raise Exception('Some files share the same basename, which can lead to error during database creation.')

    # create of nodes and names
    taxdir = os.path.join(outdir, "taxonomy")
    os.makedirs(taxdir)

    print ("Generate NCBI tree")
    if not gtdbtk_res: 
        emptyfname = os.path.join(taxdir, 'empty.gtdbtkres.tsv')
        gtdbtk2ncbi.create_empty(fnames, emptyfname)
        gtdbtk_res = [emptyfname]
    
    taxids = gtdbtk2ncbi.main(gtdbtk_res, taxdir, nodes, names, ext, drep_cdb, no_prune, slevel=add_strain_level)

    # add to library
    # need to create a custom fasta file with specific header
    print ("Produce tmp fasta files")
    tmp_fdir = os.path.join(outdir, "tmp")
    os.makedirs(tmp_fdir)
    pool = Pool(processes=threads)
    res = []

    for idx, fname in enumerate(fnames):
        try: tid = taxids[os.path.basename(fname)]
        except KeyError: raise Exception(f'Basename not found {os.path.basename(fname)}, maybe you should add an extension with option --ext?')
        res.append(pool.apply_async(make_tmp_files, args=(fname, tid, idx, tmp_fdir)))

    pool.close()

    for e in res:
        # https://stackoverflow.com/a/28660669/5016055
        e.get()

    pool.join()

    print ("Add libraries")
    cmdline = "find %s -name '*.fa' -print0 | xargs -P %i -0 -I{} -n1 %s --add-to-library {} --db %s" %(
        tmp_fdir, threads, krakenb, outdir)

    print (cmdline)
    os.system(cmdline)

    print ("Remove temporary files")
    shutil.rmtree(tmp_fdir)
    
    print ("Build the database")
    cmdline = "%s --build --threads %i --db %s" %(krakenb, threads, outdir)
    print (cmdline)
    os.system(cmdline)

    print ("Build kraken inspect")
    outfile = os.path.join(outdir, "inspection.txt.gz")
    cmdline = "%s --db %s | gzip > %s" %(krakeni, outdir, outfile)
    print (cmdline)
    os.system(cmdline)

    # run bracken build
    print ("Build bracken database")
    cmdline = "%s -d %s -t %i" %(brackenb, outdir, threads)
    print (cmdline)
    os.system(cmdline)

def names_maxtaxid(names) :
    with open(names) as  f :
        return max(int(line.split()[0])
            for line in f)

def make_tmp_files(fname, taxid, idx, outdir) :
    print ("Create temp file : %s" %(fname))
        
    fdata = list(SeqIO.parse(fname, "fasta"))
    for record in fdata :
        record.id = record.name = record.id + "|kraken:taxid|" + str(taxid)
    
    # Write, process and remove
    tmp_fasta = os.path.join(outdir, "tmp.fasta.%i.%i.fa" %(taxid, idx))
    SeqIO.write(fdata, tmp_fasta, "fasta")

if __name__ == '__main__':
    run()
# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-11-14 17:56:19
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-06-24 10:00:07

import glob, os, shutil
from Bio import SeqIO
from multiprocessing import Pool
from collections import defaultdict

import gtdbtk2ncbi
import click

@click.command()
@click.argument('fastas', type=str)
@click.argument('outdir', type=str)
@click.option('--krakenuniqb', default="krakenuniq-build", type=str, help="krakenuniq-build executable")
@click.option('--gtdbtk_res', '-r', type=str, multiple=True, help="GTDBtk result files")
@click.option('--nodes', type=str, default='', help="GTDB archeal metadata file")
@click.option('--names', type=str, default='', help="GTDB bacterial metadata file")
@click.option('--ext', type=str, default='', help="Fasta files extension (.fa, .fasta), usefull if you want to link to gtdbtk files")
@click.option('--drep_cdb', type=str, default='', help="DRep CBD result for novel SPECIES")
@click.option('--no-prune', is_flag=True, help='Don\'t prune the tree to keep only used nodes')
@click.option('--threads', default=1, type=int)
def run(fastas, outdir, krakenuniqb, gtdbtk_res, nodes, names, ext, drep_cdb, no_prune, threads) :
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
    
    taxids = gtdbtk2ncbi.main(gtdbtk_res, taxdir, nodes, names, ext, drep_cdb, no_prune, subset=bnames)

    # add to library
    # need to create a custom fasta file with specific header
    print ('Link and map fasta')
    
    libdir = os.path.join(outdir, "library")
    os.makedirs(libdir)
    pool = Pool(processes=threads)
    res = []

    for fname in fnames:
        try: tid = taxids[os.path.basename(fname)]
        except KeyError: raise Exception(f'Basename not found {os.path.basename(fname)}, maybe you should add an extension with option --ext?')
        res.append(pool.apply_async(link_fasta, args=(fname, tid, libdir)))

    pool.close()

    for e in res:
        # https://stackoverflow.com/a/28660669/5016055
        e.get()

    pool.join()

    print ('Build database')
    cmdline = f'{krakenuniqb} --db {outdir} --threads {threads} --taxids-for-genomes --taxids-for-sequences'
    print (cmdline)
    os.system(cmdline)

def link_fasta(fname, tid, libdir):
    outfile = os.path.join(libdir, os.path.basename(fname))
    os.symlink(fname, outfile)

    fdata = SeqIO.parse(fname, 'fasta')
    with open(outfile + '.map', 'w') as f:
        for record in fdata:
            f.write(f'{record.id}\t{tid}\t{os.path.basename(fname)}\n')

if __name__ == '__main__':
    run()
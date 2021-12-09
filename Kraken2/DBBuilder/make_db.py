# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-11-14 17:56:19
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-12-09 13:00:39

import click
import glob, os, sys, shutil, csv

from Bio import SeqIO
from multiprocessing import Pool
from collections import defaultdict

bname = os.path.basename

@click.command()
@click.argument('fastas', type=str)
@click.argument('outdir', type=str)
@click.option('--rnode', default=2, type=int, help="mother taxid node")
@click.option('--taxidsfile', default='', type=str, help="fname to taxid map, see readme")
@click.option('--ignore_missing_taxids', default=False, type=bool, is_flag=True, help="ignore missing taxids")
@click.option('--taxonomy', default="", type=str, help="taxonomy directory")
@click.option('--krakenb', default="kraken2-build", type=str, help="kraken2-build path")
@click.option('--krakeni', default="kraken2-inspect", type=str, help="kraken2-inspect path")
@click.option('--brackenb', default="bracken-build", type=str, help="bracken-build path")
@click.option('--startid', default=0, type=int, help="The starting value for your taxonomic ids")
@click.option('--threads', default=1, type=int)
def run(fastas, outdir, rnode, taxidsfile, ignore_missing_taxids, taxonomy, krakenb, krakeni, brackenb, startid, threads) :
    
    if os.path.isdir(outdir) :
        raise Exception("Directory already exists : %s" %(outdir))

    fnames = glob.glob(fastas)
    print ("%i files to use" %(len(fnames)))
    assert fnames

    # look for taxonomic dir
    taxonomy = taxonomy or "./taxonomy"
    if not os.path.isdir("taxonomy") :
        raise Exception("taxonomy directory not found : %s" %(taxonomy))

    # look for nodes and names
    nodes = os.path.join(taxonomy, "nodes.dmp")
    if not os.path.isfile(nodes) : 
        raise Exception("Nodes file not found in taxonomy directory : %s" %(nodes))

    names = os.path.join(taxonomy, "names.dmp")
    if not os.path.isfile(names) : 
        raise Exception("Names file not found in taxonomy directory : %s" %(names))

    if not startid :
        print ("start taxid (startid) not given, search from names file ...")
        maxtaxid = names_maxtaxid(names)
        startid = maxtaxid + 1
        print ("found maxtaxid %i, startid sets to %i" %(maxtaxid, startid))

    # create of nodes and names
    taxdir = os.path.join(outdir, "taxonomy")
    os.makedirs(taxdir)

    # copy files
    print ("Copy nodes and names ...")
    shutil.copyfile(names, os.path.join(taxdir, "names.dmp"))
    shutil.copyfile(nodes, os.path.join(taxdir, "nodes.dmp"))
    print ("done")

    print ("Update names and nodes files")
    if taxidsfile: taxids = read_taxidsfile(taxidsfile, fnames, startid, ignore_missing_taxids)
    else: taxids = {idx: [fname] for idx, fname in enumerate(fnames, start=startid)}

    # names
    outfile = os.path.join(taxdir, "names.dmp")
    with open(outfile, "a") as f :
        for tid, fnames in taxids.items():
            name = bname(fnames[0])
            f.write("%s\t|\t%s\t|\t\t|\tscientific name\t|\n" %(tid, name))

    # nodes
    outfile = os.path.join(taxdir, "nodes.dmp")
    with open(outfile, "a") as f :
        for tid in taxids:
            f.write("%i\t|\t%i\t|\tspecies\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|\n" %(tid, rnode))

    # add to library
    # need to create a custom fasta file with specific header
    print ("Produce tmp fasta files")
    tmp_fdir = os.path.join(outdir, "tmp")
    os.makedirs(tmp_fdir)
    pool = Pool(processes=threads)

    for tid, fnames in taxids.items():
        for idx, fname in enumerate(fnames):
            pool.apply_async(make_tmp_files, args=(fname, tid, idx, tmp_fdir))

    pool.close()
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

def read_taxidsfile(taxidsfile, fnames, startid, ignore_missing_taxids):
    fnames = set(fnames)
    taxids = defaultdict(list)

    with open(taxidsfile) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            if len(row) != 2: continue
            fname, tid = row

            if fname not in fnames:
                print (fname, 'not found in fnames, ignore')
                continue

            tid = int(tid) + startid
            taxids[tid].append(fname)

    missing = set(fnames) - {fname for fnames in taxids.values() for fname in fnames}

    if missing:
        print (len(missing), 'missing fnames ..., first one: ', sorted(missing)[0])
        if not ignore_missing_taxids:
            raise Exception('Missing taxids raised an error')

    return taxids


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
import click
import glob, os, sys, shutil

from Bio import SeqIO
from multiprocessing import Pool

@click.command()
@click.argument('fastas', type=str)
@click.argument('outdir', type=str)
@click.argument('rnode', type=int)
@click.option('--taxonomy', default="", type=str, help="taxonomy directory")
@click.option('--krakenb', default="kraken2-build", type=str, help="kraken2-build path")
@click.option('--krakeni', default="kraken2-inspect", type=str, help="kraken2-inspect path")
@click.option('--brackenb', default="bracken-build", type=str, help="bracken-build path")
@click.option('--startid', default=0, type=int, help="The starting value for your taxonomic ids")
@click.option('--threads', default=1, type=int)
def run(fastas, outdir, rnode, taxonomy, krakenb, krakeni, brackenb, startid, threads) :
    
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
    taxids = {fname : idx for idx, fname in enumerate(fnames, start=startid)}
    shorts = {fname : os.path.basename(fname) for fname in fnames}
    separator = "\t|\t"

    # names
    outfile = os.path.join(taxdir, "names.dmp")
    with open(outfile, "a") as f :
        for fname in fnames :
            tid = taxids[fname]
            name = shorts[fname]
            f.write("%s\t|\t%s\t|\t\t|\tscientific name\t|\n" %(tid, name))

    # nodes
    outfile = os.path.join(taxdir, "nodes.dmp")
    with open(outfile, "a") as f :
        for fname in fnames :
            print (fname, tid)
            tid = taxids[fname]
            f.write("%s\t|\t2\t|\tspecies\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|\n" %(tid))

    # add to library
    # need to create a custom fasta file with specific header
    print ("Produce tmp fasta files")
    tmp_fdir = os.path.join(outdir, "tmp")
    os.makedirs(tmp_fdir)
    pool = Pool(processes=threads)
    for fname in fnames :
        pool.apply_async(make_tmp_files, args=(fname, taxids[fname], tmp_fdir))

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
    print ("Build braken database")
    cmdline = "%s -d %s -t %i" %(brackenb, outdir, threads)
    print (cmdline)
    os.system(cmdline)


def names_maxtaxid(names) :
    with open(names) as  f :
        return max(int(line.split()[0])
            for line in f)

def make_tmp_files(fname, taxid, outdir) :
    print ("Create temp file : %s" %(fname))
        
    fdata = list(SeqIO.parse(fname, "fasta"))
    for record in fdata :
        record.id = record.name = record.id + "|kraken:taxid|" + str(taxid)
    
    # Write, process and remove
    tmp_fasta = os.path.join(outdir, "tmp.fasta.%i.fa" %(taxid))
    SeqIO.write(fdata, tmp_fasta, "fasta")

if __name__ == '__main__':
    run()

# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-06-23 11:38:15
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-06-24 09:45:21

import os, csv
import click
from collections import Counter
from dataclasses import dataclass

@dataclass
class Node:
    tax_id: int
    parent_tax_id: int
    rank: str
    novel: bool = False

    def __post_init__(self):
        self.tax_id = int(self.tax_id)
        self.parent_tax_id = int(self.parent_tax_id)

    def get_row(self, sep='\t', end='\n'):
        return sep.join((str(self.tax_id), str(self.parent_tax_id), self.rank)) + end

@dataclass
class Name:
    tax_id: int
    name_txt: str
    unique_name: str
    name_class: str

    def __post_init__(self):
        self.tax_id = int(self.tax_id)

    def get_row(self, sep='\t', end='\n'):
        return sep.join((str(self.tax_id), self.name_txt, self.unique_name, self.name_class)) + end

@dataclass
class Seq:
    seqname: str
    tax_id: int
    
    def get_row(self, sep='\t', end='\n'):
        return sep.join((self. seqname, str(self.tax_id))) + end

RANKS = 'sgfocpdr'
RANKSSET = set(RANKS)
RANKSDEF = {
    's': 'species', 'g': 'genus', 'f': 'family', 'o': 'order',
    'c': 'class', 'p': 'phylum', 'd': 'division'
}

NOVEL_IDX = 1

# ------------------------------------------------------------------------------------------------

def check_file(fname):
    if fname and not os.path.isfile(fname):
        raise IOError(f'File not found: {fname}')

def check_dir(dname):
    if not os.path.isdir(dname):
        raise IOError(f'outdir does not exist: {dname}')

def parse_classification(classification):
    classification = classification.split(';')
    return dict(element.split('__') for element in classification)

def parse_drep(fname):
    data = {}
    with open(fname) as f:
        for line in csv.DictReader(f, delimiter=','):
            data[line['genome']] = line['secondary_cluster']
    return data

def parse_gtdbtk_summary(* fnames):
    for fname in fnames:
        with open(fname) as f:
            for line in csv.DictReader(f, delimiter='\t'):
                yield line['user_genome'], line['classification']

def parse_ncbi(fname):
    data = []
    with open(fname) as f:
        for line in f:
            line = line.strip().split('\t|\t')
            if line[-1].endswith('|'): line[-1] = line[-1][:-1].strip()
            yield line

def rec_check(classification, txt2tid, nodes, names, rankidx=0):
    '''
    Recursively check from species to division if the value provided
    is found in current node. If not, a new node is created with a taxid
    and current node is linked to parent node. Should be called with the leaf
    node (species), which is defined with rankdidx == 0 as specified by RANKS.
    '''

    rank = RANKS[rankidx]
    if rank == 'r' : return 1 # root
    
    value = classification[rank]
    if not value: 
        global NOVEL_IDX
        value = 'novel' + str(NOVEL_IDX)
        NOVEL_IDX += 1
    
    key = rank + '__' + value
    taxid = txt2tid.get(key, None)

    if not taxid:
        parent = rec_check(classification, txt2tid, nodes, names, rankidx+1)
        taxid = nodes[-1].tax_id + 1 # last node (sorted) taxid + 1
        nodes.append(Node(taxid, parent, RANKSDEF[RANKS[rankidx]], True))
        names.append(Name(taxid, key, '', 'scientific name'))
        txt2tid[key] = taxid
        print (f'- Element not found: {key}, node create with taxid {taxid} and parent taxid {parent}')

    return taxid

def fetch_taxids(taxids, nodes, res=None):
    '''
    Return all taxid, including parent nodes
    from a list of taxids
    '''

    if res is None: res = set()
    ptid = {node.tax_id: node.parent_tax_id for node in nodes}

    def rec_taxid(taxid, ptid, res):
        if taxid in res: return
        res.add(taxid)
        if taxid in ptid: 
            rec_taxid(ptid[taxid], ptid, res)

    for taxid in taxids:
        rec_taxid(taxid, ptid, res)

    return res

def write_dataclasses(data, outfile, sep, end='\n'):
    with open(outfile, 'w') as f:
        for element in data:
            f.write(element.get_row(sep, end))

# ------------------------------------------------------------------------------------------------

def main(gtdbtk_res, outdir, nodes, names, ext, drep_cdb, no_prune, subset=None) :
    check_file(nodes)
    check_file(names)
    check_dir(outdir)

    if names:
        names = [Name(* row[:4]) for row in parse_ncbi(names)]
    else:
        names = [Name(1, 'root', '', 'scientific name')]

    if nodes:
        nodes = [Node(* row[:3]) for row in parse_ncbi(nodes)]
    else:
        nodes = [Node(1, 1, 'root')]

    if drep_cdb:
        drep_cdb = parse_drep(drep_cdb)
    else:
        drep_cdb = None

    nodes = sorted(nodes, key=lambda element: element.tax_id)
    txt2tid = {name.name_txt: name.tax_id for name in names}
    seqtids = []

    for name, classification in parse_gtdbtk_summary(* gtdbtk_res):
        name = name + ext

        if subset and name not in subset:
            continue

        classification = parse_classification(classification)

        if not classification['s'] and drep_cdb is not None:
            try:
                classification['s'] = 'drep_' + drep_cdb[name]
            except KeyError:
                raise Exception(f'Unable to find genome `{name}` in DRep CDB file, maybe the exception is missing?')

        assert len(set(classification) - RANKSSET) == 0
        taxid = rec_check(classification, txt2tid, nodes, names)
        seqtids.append(Seq(name, taxid))


    if subset:
        missing = set(subset) - {seq.seqname for seq in seqtids}
        if missing:
            raise Exception(f'{len(missing)} sequences were not found in the GATK files, first: {next(iter(missing))}')

    nodes = sorted(nodes, key=lambda element: element.parent_tax_id)
    used_taxids = fetch_taxids((seq.tax_id for seq in seqtids), nodes)

    counter_all = Counter(node.rank for node in nodes)
    counter_novel = Counter(node.rank for node in nodes if node.novel)
    counter_used = Counter(node.rank for node in nodes if node.tax_id in used_taxids)

    print ('Statistic ...')
    print (f'{len(seqtids)} sequences provided')
    for rank in counter_all:
        print (f'{rank} - All: {counter_all[rank]} - Novel: {counter_novel[rank]} - Used: {counter_used[rank]}')

    if no_prune == False:
        nodes = (node for node in nodes if node.tax_id in used_taxids)
        names = (name for name in names if name.tax_id in used_taxids)

    write_dataclasses(nodes, os.path.join(outdir, 'nodes.dmp'), sep='\t|\t')
    write_dataclasses(names, os.path.join(outdir, 'names.dmp'), sep='\t|\t', end='\t|\n')
    write_dataclasses(seqtids, os.path.join(outdir, 'seq2taxid.tsv'), sep='\t')

    return {seq.seqname: seq.tax_id for seq in seqtids}

def create_empty(fnames, outfile):
    '''
    Create an empty gtdbtk_res based on file names
    '''
    s = 'd__;p__;c__;o__;f__;g__;s__'
    with open(outfile, 'w', newline='') as f:
        fieldnames = ['user_genome', 'classification']
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for fname in fnames:
            name = os.path.basename(fname)
            writer.writerow({'user_genome': name, 'classification': s})
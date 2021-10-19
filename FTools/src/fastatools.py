# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2020-12-08 11:27:43
# @Last modified by:   jsgounot
# @Last Modified time: 2021-10-19 16:43:47

import os, sys

import click

import ametric as ametricc
import compare as comparec
import extract as extractc
import generate as generatec
import info as infoc

@click.group()
def cli():
    pass

@cli.command(help="Display assembly metrics")
@click.argument("files", nargs=-1, type=str)
@click.option('--ref', type=str, help="Reference fasta path (LG stats)", default="")
@click.option('--refsize', type=int, help="Override --ref", default=0)
@click.option('--minsize', type=int, help="Sequence minimum size to consider", default=0)
@click.option('--maxsize', type=int, help="Sequence maximum size to consider", default=0)
@click.option('--nvalue', type=int, help="N value to consider (N50, N90). Default 50", default=50)
@click.option('--nsep', help="Add separator to numbers", is_flag=True)
@click.option('--fullname', help="Replace basename with full name", is_flag=True)
@click.option('--cpath', help="Try to complete the relative path with the absolute path", is_flag=True)
@click.option('--outfile', type=str, help="Path of the TSV table (default stdout)", default="")
@click.option('--force', help="Force analyze of unlikely fasta file", is_flag=True)
@click.option('--ncore', type=int, help="Number of core to use. Default 1", default=1)
def ametric(* args, ** kwargs) :
    ametricc.run(* args, ** kwargs)

@cli.command(help="Compare two fasta files")
@click.argument("fasta1", type=str)
@click.argument("fasta2", type=str)
@click.option('--minsize', type=int, help="Sequence minimum size to consider", default=0)
@click.option('--maxsize', type=int, help="Sequence maximum size to consider", default=0)
@click.option('--sort', type=str, help="Column name used to sort results", default="ID")
@click.option('--names', type=str, multiple=True, help="Extract only provided header names")
@click.option('--force', help="Force analyze of unlikely fasta file", is_flag=True)
def compare(* args, ** kwargs) :
    comparec.run(* args, ** kwargs)

@cli.command(help="Extract fasta sequences")
@click.argument("files", nargs=-1, type=str)
@click.option('--minsize', type=int, help="Sequence minimum size to consider", default=0)
@click.option('--maxsize', type=int, help="Sequence maximum size to consider", default=0)
@click.option('--start', type=int, help="Start position", default=-1)
@click.option('--end', type=int, help="End position", default=-1)
@click.option('--names', type=str, multiple=True, help="Extract only provided header names - Use : --names name1 name2 name3")
@click.option('--revcomp', help="Reverse complete DNA sequence", is_flag=True)
@click.option('--ncount', type=int, default=0, help="Sample only a number n of sequences from the output")
@click.option('--rand', help="Randomly selected sequences, can be coupled with ncount to generate random sampling", is_flag=True)
@click.option('--correct', help="Correct sequence by replacing ? to N", is_flag=True)
@click.option('--sort', type=click.Choice(['name', 'size']), help="Sequence sort order to apply. Can be either 'name' or 'size'")
@click.option('--reverse', help="Reverse the sorting order", is_flag=True)
@click.option('--force', help="Force analyze of unlikely fasta file", is_flag=True)
def extract(* args, ** kwargs) :
    extractc.run(* args, ** kwargs)

@cli.command(help="Generate random fasta sequences")
@click.option('--basename', type=str, help="Sequences basename", default="RandomSequence:")
@click.option('--start', type=int, help="Starting sequences indice", default=1)
@click.option('--number', type=int, help="Number of sequences to generate", default=1)
@click.option('--size', type=int, help="Size of generated sequence", default=100)
@click.option('--upperSize', type=int, help="If provided, sequence size will randomly range between (size, uppersize)", default=0)
@click.option('--bases', type=str, help="Bases to use", default="ATGG")
@click.option('--weights', type=str, help="Weights associated to base. Formatting : 1/1/1/1", default='1/1/1/1')
@click.option('--seed', type=int, help="Random seed", default=0)
def generate(* args, ** kwargs) :
    generatec.run(* args, ** kwargs)

@cli.command(help="Display sequences informations")
@click.argument("files", nargs=-1, type=str)
@click.option('--minsize', type=int, help="Sequence minimum size to consider", default=0)
@click.option('--maxsize', type=int, help="Sequence maximum size to consider", default=0)
@click.option('--sort', type=str, help="Column name used to sort results", default="ID")
@click.option('--head', type=int, help="N first rows to display", default=0)
@click.option('--names', type=str, multiple=True, help="Extract only provided header names")
@click.option('--force', help="Force analyze of unlikely fasta file", is_flag=True)
@click.option('--basecount', help="Return basecount instead", is_flag=True)
@click.option('--nsep', help="Add separator to numbers", is_flag=True)
@click.option('--asone', help="Produce a single dataframe", is_flag=True)
@click.option('--fullname', help="Replace basename with full name, imply --asone", is_flag=True)
@click.option('--outfile', type=str, help="Path of the TSV table (default stdout), imply --asone", default="")
def info(* args, ** kwargs) :
    infoc.run(* args, ** kwargs)

cli = click.CommandCollection(sources=[cli],
    help="Tools collection to work with fasta files")

if __name__ == '__main__':
    cli()
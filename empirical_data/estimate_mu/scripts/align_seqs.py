#!/usr/bin/env python3
import argparse
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
from Bio.Seq import Seq
import numpy as np
import subprocess
import shlex
import random
from collections import Counter
from datetime import date
import os
import sys


def align_seqs(seqs_path):
    align_cmd = \
        f"mafft --auto --reorder --anysymbol --nomemsave --adjustdirection --preservecase --thread 4 {seqs_path}" 
    aligned_path = os.path.splitext(seqs_path)[0] + '_aligned.fasta'
    subprocess.run(shlex.split(align_cmd), stdout=open(aligned_path, 'w'))
    return(aligned_path)


def align_seqs_to_ref(seqs_path, reference_name):
    # aligns relative
    seqs = list(SeqIO.parse(seqs_path, 'fasta'))
    ref_seq = [i for i in seqs if reference_name in i.description]
    if len(ref_seq) > 1:
        raise Exception('more than one match to reference name')
    seqs = [i for i in seqs if reference_name not in i.description]
    with open(os.path.splitext(seqs_path)[0] + '_noref.fasta', 'w') as out:
        SeqIO.write(seqs, out, 'fasta')
    with open(os.path.splitext(seqs_path)[0] + '_ref.fasta', 'w') as out:
        SeqIO.write(ref_seq, out, 'fasta')
    align_cmd = \
        f"mafft --auto --keeplength --preservecase --addfragments  \
            {os.path.splitext(seqs_path)[0] + '_noref.fasta'} \
            {os.path.splitext(seqs_path)[0] + '_ref.fasta'}" 
    aligned_path = os.path.splitext(seqs_path)[0] + '_aligned_ref.fasta'
    subprocess.run(shlex.split(align_cmd), stdout=open(aligned_path, 'w'))
    return(aligned_path)


def degap_seqs(seqs_path, reference_name):
    seqs = list(SeqIO.parse(seqs_path, 'fasta'))
    ref = [i for i in seqs if reference_name in i.name]
    if len(ref) > 1:
        raise Exception('multiple sequences match reference name')
    else:
        aligned_ref = ref[0]
    aligned_ref_seq = list(aligned_ref.seq)
    aligned_ref_gaps = \
        set([i for i, item in enumerate(aligned_ref_seq) if item == '-'])
    for item in seqs:
        item.seq = item.seq.upper()
        item_seq = list(item.seq)
        item_seq = [item for i, item in enumerate(item_seq) 
            if i not in aligned_ref_gaps]
        item.seq = Seq(''.join(item_seq))
    degapped_path = os.path.splitext(seqs_path)[0] + '_degap.fasta'
    with open(degapped_path, 'w') as out:
        SeqIO.write(seqs, out, 'fasta')
    return(degapped_path)


# lots of this from augur
# https://github.com/nextstrain/augur/blob/master/augur/filter.py#L221
def filter_seqs(aligned_path=None, min_length=None):
    aligned_seqs = list(SeqIO.parse(aligned_path, 'fasta'))
    # sequences to force keep even if they fail filtering
    # 2) filter by minimum length
    # as is done in the nextstrain pipeline, this is the sum of ACTG: 
    # https://github.com/nextstrain/augur/blob/efcc8be3bcf72ab4beb237afcc5b4b002bbaa806/augur/filter.py#L221
    if min_length:
        prev_n = len(aligned_seqs)
        aligned_seqs = [item for item in aligned_seqs if 
            sum(map(lambda x: item.seq.count(x), 
                ["a", "t", "g", "c", "A", "T", "G", "C"])) > min_length]
        print(f'{len(aligned_seqs) - prev_n} sequences removed by length filter', file=sys.stderr)
    filtered_path = os.path.splitext(aligned_path)[0]+'_filtered.fasta'
    with open(filtered_path, 'w') as out:
        SeqIO.write(aligned_seqs, out, 'fasta')
    return(filtered_path)


def mask_seqs(seq_path, to_mask):
    seqs = list(SeqIO.parse(seq_path, 'fasta'))
    for item in seqs:
        seq_list = list(item.seq)
        if to_mask['head']:
            seq_list[0:to_mask['head']] = list('N' * to_mask['head'])
        if to_mask['tail']:
            seq_list[-to_mask['tail']:] = list('N' * to_mask['tail'])
        if to_mask['sites']:
            for site in to_mask['sites']:
                seq_list[site-1] = 'N'
        item.seq = Seq(''.join(seq_list))
    outname = os.path.splitext(seq_path)[0] + '_masked.fasta'
    with open(outname, 'w') as out:
        SeqIO.write(seqs, out, 'fasta')
    return(outname)


def run():
    parser = argparse.ArgumentParser()
    # input files
    # todo outdir currently doesn't get used by anything
    parser.add_argument('--outdir', 
        default='.')
    parser.add_argument('--sequences', 
            default=None,
            help='path to fasta file to align')
    parser.add_argument('--referenceName', 
            default=None,
            help='name of the reference sequences')
    parser.add_argument('--minLength',
        default=0,
        type=int,
        help='minimum number of ACTG characters')
    parser.add_argument('--maskHead',
        default=0,
        type=int,
        help='number of characters to mask at the begining')
    parser.add_argument('--maskTail',
        default=0,
        type=int,
        help='number of characters to mask at end')
    parser.add_argument('--maskSites',
        default=[],
        type=int,
        nargs='+',
        help='specific sites to mask')
    parser.add_argument('--alignType',
        default='pairwise',
        help='one of: "pairwise" or "reference" to either do pairwise alignment or align to reference')
    args = parser.parse_args()
    to_mask = {} 
    to_mask['head'] = args.maskHead 
    to_mask['tail'] = args.maskTail
    to_mask['sites'] = args.maskSites
    #args.alignType = args.alignType.lower()
    if args.alignType not in ['pairwise', 'reference']:
        raise Exception('alignType must be one of pairwise or reference')
    if args.alignType == 'pairwise':
        aligned_seqs_path = align_seqs(args.sequences)
        degapped_seqs_path = \
        degap_seqs(aligned_seqs_path, args.referenceName)
    elif args.alignType == 'reference':
        aligned_seqs_path = align_seqs_to_ref(args.sequences, args.referenceName)
        degapped_seqs_path = aligned_seqs_path
    # quality filtering
    # sequences must be aligned for this filtering to work
    filtered_seqs_path = \
        filter_seqs(aligned_path=degapped_seqs_path, min_length=args.minLength)
    if args.maskHead or args.maskTail or args.maskSites:
        filtered_seqs_path = \
            mask_seqs(filtered_seqs_path, to_mask)
    print(filtered_seqs_path, file=sys.stdout)
    # removes reference so that we have a reference free alignment
    filtered_seqs = list(SeqIO.parse(filtered_seqs_path, 'fasta'))
    filtered_seqs = [i for i in filtered_seqs if args.referenceName not in i.description]
    with open(filtered_seqs_path.split('.')[0]+'_noref.fasta', 'w') as out:
        SeqIO.write(filtered_seqs, out, 'fasta')


if __name__ == "__main__":
    run()

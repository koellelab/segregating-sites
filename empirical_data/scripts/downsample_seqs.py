import argparse
from Bio import SeqIO
import numpy as np
import pandas as pd
import datetime
from utils import convert_to_datetime
#from scripts.utils import convert_to_datetime
import sys 
from collections import Counter


def run():
    parser = argparse.ArgumentParser()
    # input files
    # //// FASTA FILE PARSING ////
    parser.add_argument('--sequences', 
        default=None)
    parser.add_argument('--seqNameSep',
        default='|')
    parser.add_argument('--seqNameField',
        default=1, type=int, 
        help='which field in the sequence name corresponds to the metadata name column')
    # //// METADATA FILE PARSING ////
    parser.add_argument('--metadata', 
        default=None, help='metadata file path, assumes no header')
    parser.add_argument('--metadataDelim',
        default='\t', 
        help='metadata delimiter')
    parser.add_argument('--metadataIDCol',
        default=1,
        type=int,
        help='which column in the metadata file contains sequence names (corresponds to seqNameField)')
    parser.add_argument('--targetN', 
        default=100, type=int, help='total number of sequences in alignment to target when using weights')
    parser.add_argument('--goodChars', 
        help='which characters to allow in fasta file',
        default=['A', 'C', 'G', 'T', '-', 'N', 'R', 'Y', 'S', 'W', 
        'K', 'M', 'D', 'H', 'B', 'V', '?'],
        nargs='+')
    parser.add_argument('--outName', 
        help='outname to save data to if None, defaults to modified input fasta name',
        default=None)
    args = parser.parse_args()
    #args.sequences = 'data/gisaid_hcov-19_2021_01_26_aligned_ref_filtered_masked.fasta'
    #args.metadata = './data/metadata_aligned.tsv'
    #args.include = './config/include.tsv'
    #args.exclude = './config/exclude.csv'
    #args.maxDate = '03-31-2020'
    #args.regionWeights = 'data/country_case_weights.tsv'
    #args.targetN = 6000
    #args.interRegionWeights = 'data/gisaid_hcov-19_2021_01_26_EHC_aligned_ref_filtered_masked_min_dist_weights.tsv'
        # read in metadata
    metadata=pd.read_csv(args.metadata, 
        sep=args.metadataDelim, header=None, low_memory=False)
    # read in sequences 
    seqs = list(SeqIO.parse(args.sequences, 'fasta'))
    # //// 1) removes sequences from sequence file with bad characters ////
    # todo faster
    # todo move up
    good_chars_str = ''.join(set(args.goodChars))
    len_seqs = len(seqs)
    seqs = [i for i in seqs if not i.seq.strip(good_chars_str)]
    print(f'{len_seqs - len(seqs)} sequences removed for non nucleotide characters', 
        file=sys.stderr)
    seqs_names = set([i.description.split(args.seqNameSep)[args.seqNameField] for i in seqs])
    
    # //// 2) remove metadata entires not in the fasta file ////
    len_metadata = metadata.shape[0]
    metadata = metadata[metadata[args.metadataIDCol].isin(seqs_names)]
    print(f'{len_metadata - metadata.shape[0]} metadata entries removed because they are not in the fasta file', 
        file=sys.stderr)    
    
    # //// 6) finally, downsample metadata ///
    # args.SamplesPerWeek and args.weights are currently incompatible
    sampled_metadata = metadata.sample(args.targetN)
    sampled_seq_names = set(sampled_metadata[args.metadataIDCol])
    print(f'{len(sampled_seq_names)} samples selected for downsampling', file=sys.stderr)
    # want to remove any sequences that are duplicated in the fasta
    sampled_seqs = []
    included_seqs_descriptions = set([])
    for idx, i in enumerate(seqs):
        if i.description.split(args.seqNameSep)[args.seqNameField] in sampled_seq_names:
            if i.description.split(args.seqNameSep)[args.seqNameField] not in included_seqs_descriptions:
                sampled_seqs.append(i)
                included_seqs_descriptions.add(i.description.split(args.seqNameSep)[args.seqNameField])
    print(f'{len(sampled_seq_names) - len(sampled_seqs)} samples missing in fasta file', file=sys.stderr)
    print(f'{len(sampled_seqs)} samples selected by downsampling', file=sys.stderr)
    # remove spaces from sequence names for IQ-Tree
    # todo confirm this works
    for i in sampled_seqs:
        i.description = i.description.replace(' ', '_')
        i.name = i.description
        i.id = i.description
    if args.outName: 
        print(args.outName+'.fasta')
        with open(args.outName+'.fasta', 'w') as out:
            SeqIO.write(sampled_seqs, out, 'fasta')
    else:
        print(args.sequences.replace('.fasta', f'_{len(sampled_seqs)}.fasta'), file=sys.stdout)
        with open(args.sequences.replace('.fasta', f'_{len(sampled_seqs)}.fasta'), 'w') as out:
            SeqIO.write(sampled_seqs, out, 'fasta')



if __name__ == "__main__":
    run()

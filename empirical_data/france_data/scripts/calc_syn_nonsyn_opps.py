import numpy as np
import pandas as pd
import argparse
from collections import defaultdict, Counter

def n_syn_nonsyn(codon, codon_dict):
    nucs = ['A', 'C', 'T', 'G']
    # todo faster in numpy
    codon_list = list(codon)
    alt_codons = []
    for pos in [0,1,2]:
        for key in nucs:
            new_codon = codon_list.copy()
            new_codon[pos] = key
            if new_codon != codon_list:
                alt_codons.append(''.join(new_codon))
    codon_dict_copy = codon_dict.copy()
    syn = [codon_dict[i] == codon_dict_copy[codon] for i in alt_codons]
    n_total = len(alt_codons)
    n_syn = sum(syn)
    n_nonsyn = len(alt_codons) - n_syn
    return(pd.DataFrame([codon, n_total,n_syn, n_nonsyn], index=['codon', 'n_total', 'n_syn', 'n_nonsyn']).T)


def import_genome(gfile):
    with open(gfile, 'r') as fp:
        ref = [i.strip() for i in fp.readlines()]
    # get nucleotide sequence
    ref_seq = np.array(list(''.join([''.join(k.split(' ')[1:]) for
        k in ref[[idx for idx, i in enumerate(ref) if 
        i[:6] == 'ORIGIN'][0]:]]).upper()))
    cds = [ref[idx:idx+2] for idx, i in enumerate(ref) if i[:3] == 'CDS']
    cds = [[i[1].split('=')[1].replace('"',''), i[0].split()[1]] for i in cds]
    cds = [[i[0], i[1].replace('join(', '').replace(')', '').split(',')] for i in cds]
    cds = [[i[0], [[int(j) for j in k.split('..')] for k in i[1]]] for i in cds]
    # +1 because GenBank coordinates are INCLUSIVE
    cds = [[i[0], np.hstack([np.arange(k[0], k[1]+1) for k in i[1]])] for i in cds]
    cds = [[i[0], pd.DataFrame(i[1]).rename(columns={0: 'pos'})] for i in cds]
     # add nucleotides and create dict
    cds_dict = {}
    for i in cds:
        # gff file is 1-indexed
        i[1]['nuc'] = ref_seq[i[1]['pos'] - 1]
         # sometimes the same ORF is present multiple times
        # take the longest if so
        for i in cds:
            if i[0] not in cds_dict.keys():
                cds_dict[i[0]] = i[1]
            elif i[1].shape[0] > cds_dict[i[0]].shape[0]:
                cds_dict[i[0]] = i[1]
    return(cds_dict)


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genome', default=None)
    args = parser.parse_args()

    #args.genome='data/NC_045512.2.gb'

    # DEFINE CODON DICTIONARY
    # TODO READ FROM FILE 
    codon_dict=defaultdict(lambda: "NaN")
    codon_dict.update(
        {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 
        'TAT': 'Y', 'TAC': 'Y', 'TAA': 'X', 'TAG': 'X', 
        'TGT': 'C', 'TGC': 'C', 'TGA': 'X', 'TGG': 'W', 
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 
        '---': '-'})

    # for each codon count # of syn and nonsyn mutations
    mut_types = pd.concat([n_syn_nonsyn(codon, codon_dict) for 
        codon in codon_dict.keys() if 
        codon != '---'])

    # read in CDS from genome
    cds_dict = import_genome(args.genome)

    # chunk into codons and count # of codons
    codon_counts = pd.concat([
        pd.DataFrame.from_dict(
            dict(Counter([''.join(i) for i in 
                np.split(val['nuc'].values, int(val['nuc'].values.shape[0]/3))])), orient='index') for 
        key, val in cds_dict.items()]).reset_index().groupby('index').agg(sum).reset_index()

    # merge
    codon_counts = codon_counts.merge(mut_types, left_on='index', right_on='codon')
    # sum everythign up
    n_syn = sum(codon_counts[0] * codon_counts['n_syn'])
    n_nonsyn = sum(codon_counts[0] * codon_counts['n_nonsyn'])
    with open('.'.join(args.genome.split('.')[:-1]) + '_n_syn_n_nonsyn.tsv', 'w') as fp:
        fp.write('\tn_syn\tn_nonsyn\n')
        fp.write(args.genome.split('/')[-1].split('.')[0] + 
            '\t' + str(n_syn) + '\t' + str(n_nonsyn))



if __name__ == "__main__":
    run()




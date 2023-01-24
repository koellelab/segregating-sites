import numpy as np
import pandas as pd
import argparse
from collections import defaultdict, Counter
import json

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


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def format_seqs_arr(s, n_seqs, nuc_dict={97:[97], 99:[99], 116:[116], 103:[103], 110: [97, 99, 116, 103, 110]}):
    n_seqs = int(n_seqs)
    size = int(len(s)/n_seqs)
    seqs_arr = \
        np.frombuffer(s.lower().encode(), dtype=np.int8)
    seqs_arr = np.copy(seqs_arr)
    #[(97, 'A'), (114, 'R'), (119, 'W'), (109, 'M'), (100, 'D'), (104, 'H'), (118, 'V'), 
    #(110, 'N'), (99, 'C'), (121, 'Y'), (115, 'S'), (109, 'M'), (98, 'B'), (104, 'H'), 
    #(118, 'V'), (110, 'N'), (117, 'U'), (121, 'Y'), (119, 'W'), (107, 'K'), (98, 'B'), 
    #(100, 'D'), (104, 'H'), (110, 'N'), (103, 'G'), (114, 'R'), (115, 'S'), (107, 'K'), 
    #(98, 'B'), (100, 'D'), (118, 'V'), (110, 'N'), (116, 'T')]
    all_nucs = np.hstack(list(nuc_dict.values()))
    seqs_arr[~np.in1d(seqs_arr, all_nucs)] = 110
    seqs_arr = \
        seqs_arr.reshape((n_seqs, int(seqs_arr.shape[0]/n_seqs)))
    return(seqs_arr)


def import_fasta(fasta_path, nuc_dict={97:[97], 99:[99], 116:[116], 103:[103], 110: [97, 99, 116, 103, 110]}):
    s_names = []
    all_s = ''
    fh = open(fasta_path, 'rt')
    with fh as fasta:
        for h,s in read_fasta(fasta):
            s_names.append(h)
            all_s += s
    fh.close()
    s_arr = format_seqs_arr(all_s, len(s_names), nuc_dict=nuc_dict)
    return(s_arr, s_names)


def import_cds(gfile):
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


# somethign is off with default nuc dict
def get_seg_sites(seqs, nuc_dict={97:np.array([97]), 99:np.array([99]), 116:np.array([116]), 103:np.array([103]), 110: np.array([97, 99, 116, 103, 110])}):
    if seqs.shape[0] >= 1:
        u,inv = np.unique(seqs,return_inverse = True)
        nucs = np.array([nuc_dict[x] for x in u], dtype=object)[inv].reshape(seqs.shape)
        l = max(max(nuc_dict.keys()), max([max(item) for item in nuc_dict.values()]))
        if not isinstance(nucs[0][0], np.ndarray):
            bins = np.array([np.bincount(i, minlength=l).max() for i in nucs.astype(np.int64).T])
        else:
            bins = np.array([np.bincount(np.concatenate(i), minlength=l).max() for i in nucs.T])
        seg_sites = (np.where(bins < seqs.shape[0])[0])
        return(seg_sites)
    else:
        return(np.nan)



def translate_seq(seq):
    # translate to protein
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
    # chunk to codons
    codons = [''.join(k).upper() for k in np.split(seq, int(seq.shape[0]/3))]
    aa = ''.join([codon_dict[i] for i in codons])
    return(aa)


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genome', default=None)
    parser.add_argument('--seqs', default=None)
    parser.add_argument('--nucDict', 
        default='scripts/nuc_dict_all.json',
        help='which nucleotide dictionary to use')
    args = parser.parse_args()

    #args.genome='data/NC_045512.2.gb'
    #args.seqs = 'data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs.fasta'

    nuc_dict = \
            {np.frombuffer(key.lower().encode(), dtype=np.int8)[0]: 
                np.frombuffer(''.join(value).lower().encode(), dtype=np.int8) 
                for key, value in json.load(open(args.nucDict, 'r')).items()}


    s_arr, s_names = import_fasta(args.seqs)

    cds_dict = import_cds(args.genome)
    cds_df = pd.DataFrame({'gene': [], 'pos': [], 'nuc': []})
    for key, val in cds_dict.items():
        val['gene'] = key
        cds_df = pd.concat([cds_df, val])


    prot = cds_df.groupby('gene').apply(lambda k: translate_seq(k['nuc'].values))

    # get all segregating sites
    # IMPORT FULL NUCLEOTIDE DICTIONARY
    seg_sites = get_seg_sites(s_arr, nuc_dict)

    # for each seg site get list of nucleotides
    seg_sites_nucs = [np.unique(i) for i in s_arr[:,seg_sites].T]
    # add all possible ambiguity interpretations but exclude Ns 
    seg_sites_unambig = [np.array(list(np.unique(np.hstack([nuc_dict[i] if i != 110 else np.array([]) for 
        i in j ])).astype(np.int8).tobytes().decode().upper())) for j in seg_sites_nucs]

    # this is sloppy but it works
    nonsyn={}
    for pos, nucs in zip(seg_sites, seg_sites_unambig):
        this_pos = []
        for nuc in nucs:
            alt_cds = cds_df.copy()
            # -1 because cds coordinates are 1 indexed
            alt_cds.iloc[np.where(alt_cds['pos']-1==pos)[0],2] = nuc
            alt_prot = alt_cds.groupby('gene').apply(lambda k: translate_seq(k['nuc'].values))
            this_pos.extend(alt_prot!=prot)
        nonsyn[pos] = any(this_pos)

    sum(nonsyn.values())
    n = len(nonsyn.values())
    n_nonsyn = sum(nonsyn.values())


    with open('.'.join(args.seqs.split('.')[:-1])+'_syn_nonsyn.tsv', 'w') as fp:
        fp.write('\tn_syn\tn_nonsyn\n')
        fp.write(args.seqs.split("/")[0]+'\t'+str(n-n_nonsyn)+'\t'+str(n_nonsyn))






if __name__ == "__main__":
    run()





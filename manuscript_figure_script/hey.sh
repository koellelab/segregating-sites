#!/bin/bash

cd empirical_data/france_data

FranceSeqs=data/gisaid_hcov-19_2021_04_29_16.fasta
echo $(awk -F'\t' 'FNR==NR{hash[$2]=$3; next}{print $0"\t"hash[$1]}' \
            ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_match_largest_seqs.tsv \
             data/france_random_global_ref_figure.tsv)

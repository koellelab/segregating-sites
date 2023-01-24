#!/bin/bash

# NEW FRANCE DATA
# france data includes complete, high coverage, collection date complete seqeucnes from France 
# sampled thorugh march 17th uploaded through April 29th 2021
# global data includes complete, high coverage, collection date complete sequences 
# sampled through march 17th uploaded through April 29th 2021

# per GISAID rules we cannot reshare data
# list of accession IDs saved in:
# data/gisaid_hcov-19_2021_04_29_16_ids.tsv
# datagisaid_hcov-19_2021_04_29_16_match_largest_seqs_ids.tsv
# data/france_random_global_ref_ids.tsv
FranceSeqs=data/gisaid_hcov-19_2021_04_29_16.fasta
ReferenceSeq=data/EPI_ISL_402125.fasta
ReferenceGenome=data/NC_045512.2.gb
globalSeqs=data/gisaid_hcov-19_2021_06_12_01.fasta


# --------------------------------------------------------------------------------------#
# //// 1. ALIGN FRENCH SEQUENCE DATA ////
# --------------------------------------------------------------------------------------#
# need to add a country column to the GISAID metadata
# by splitting their Location column
awk -F'\t' '{split($4, a, " / "); print $2","$3","a[2]}' \
    ${FranceSeqs%.fasta}.tsv | tail -n +2 \
    > ${FranceSeqs%.fasta}_date_country.csv

# get lists of IDS used
grep ">" $FranceSeqs | cut -d'|' -f 2 > ${FranceSeqs%.fasta}_ids.tsv 

# align sequences to reference
# masking standard sites from early de Maio et al work 
# NOTE: not a pairwise alignment
# add reference
cat $ReferenceSeq $FranceSeqs | sed 's/ //g' > ${FranceSeqs%.fasta}_ref.fasta

alignSeqs=$(python3 scripts/align_seqs.py \
    --sequences ${FranceSeqs%.fasta}_ref.fasta \
    --referenceName EPI_ISL_402125 \
    --minLength 2800 \
    --maskHead 55 \
    --maskTail 100 \
    --maskSites 11083 15324 21575 \
    --alignType reference)

# --------------------------------------------------------------------------------------#
# //// 2. GET CONNECTED COMPONENTS ////
# --------------------------------------------------------------------------------------#
# get pairwise distances
python3 scripts/calc_pairwise_distances.py \
    --seqs ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref.fasta \
    --nucDict scripts/nuc_dict_all.json

# get the connected components from those pairwise distances
python3 scripts/get_connected_components.py \
    --dist  ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref.dist \
    --metadata ${FranceSeqs%.fasta}.tsv \
    --ref $ReferenceSeq \
    --seqs ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref.fasta \
    --nucDict scripts/nuc_dict_all.json

# --------------------------------------------------------------------------------------#
# //// 3. GET ALL SEQUENCES THAT DESCEND FROM THE LARGEST CONNECTED COMPONENT ////
# --------------------------------------------------------------------------------------#
# get most common CC
largest=$(awk '{print $2}' ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc.tsv | \
    uniq -c | sort -k1 -n - | tail -n 1 | awk '{print $2}')

# get all CCs which descend from the largest
python3 scripts/get_match_clades.py \
    --ccSNPs ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_snp.tsv \
    --cc ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc.tsv \
    --dists ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref.dist \
    --toMatch $largest \
    --maxDist  7 \
    --nucDict scripts/nuc_dict_all.json \
    --ref $ReferenceSeq 

# get sequences which descend 
python3 scripts/get_seqs.py \
    --outName ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_0 \
    --seqs ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref.fasta \
    --getSeqs ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest.tsv \
    --outName ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs

# identify synonymous and nonsynonymous mutations
python3 scripts/calc_syn_nonsyn_sites.py \
    --genome ${ReferenceGenome} \
    --seqs ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs.fasta

# calcualte syn nonsyn opportunities
python3 scripts/calc_syn_nonsyn_opps.py \
    --genome ${ReferenceGenome}

# get metadata just for those sequences
awk -F'\t' 'FNR==NR {hash[$1]; next} $2 in hash' \
    <(grep ">" ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs.fasta |\
    sed 's/\>//g' | awk -F'|' '{print $2}') \
    ${FranceSeqs%.fasta}.tsv \
    > ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_match_largest_seqs.tsv

# get list of sequences used
grep ">" ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs.fasta \
    | cut -d'|' -f 2 > ${FranceSeqs%.fasta}_match_largest_seqs_ids.tsv


# --------------------------------------------------------------------------------------#
# //// 4. CALCULATE SEGREGATING SITES OF FRANCE DATA ////
# --------------------------------------------------------------------------------------#
# everythign just equals itself, exlcuding N
python3 scripts/calc_segregating_sites.py \
    --metadata ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_match_largest_seqs.tsv \
    --seqs ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs.fasta \
    --metadataDelim \\t \
    --metadataNameCol 1 \
    --metadataDateCol 2 \
    --nucDict scripts/nuc_dict_unambig.json

mv ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s.tsv \
   ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s_unambig.tsv

# everything just equals itself, including N
python3 scripts/calc_segregating_sites.py \
    --metadata ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_match_largest_seqs.tsv \
    --seqs ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs.fasta \
    --metadataDelim \\t \
    --metadataNameCol 1 \
    --metadataDateCol 2 \
    --nucDict scripts/nuc_dict_none.json

mv ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s.tsv \
   ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s_none.tsv

# properly account for ambiguous nucleotides
python3 scripts/calc_segregating_sites.py \
    --metadata ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_match_largest_seqs.tsv \
    --seqs ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs.fasta \
    --metadataDelim \\t \
    --metadataNameCol 1 \
    --metadataDateCol 2

# --------------------------------------------------------------------------------------#
# //// 5. CALCULATE THE DISTANCE OF EACH OF THESE SEQUENCES TO WUHAN/HU-1 ////
# --------------------------------------------------------------------------------------#
python3 scripts/calc_distance.py \
    --querySeqs ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs.fasta  \
    --targetSeq data/EPI_ISL_402125.fasta \
    --nucDict scripts/nuc_dict_all.json

# --------------------------------------------------------------------------------------#
# //// 6. INCORPORATE EXOGENEOUS SEQUENCE DATA ////
# --------------------------------------------------------------------------------------#
# regrettably, need to align all global sequences for later steps
python3 scripts/align_seqs.py \
    --sequences ${globalSeqs} \
    --referenceName EPI_ISL_402125 \
    --minLength 2800 \
    --maskHead 55 \
    --maskTail 100 \
    --maskSites 11083 15324 21575 \
    --alignType reference

# get lists of IDS used
grep ">" ${globalSeqs} | cut -d'|' -f 2 > ${globalSeqs%.fasta}_ids.tsv 

# add to france sequences
cat ${globalSeqs%.fasta}_ids.tsv  ${FranceSeqs%.fasta}_ids.tsv > france_global_seqs.tsv

# downsample 100 random global sequences
python3 scripts/downsample_seqs.py \
    --sequences $globalSeqs \
    --metadata <(awk -F'\t' '{gsub (" ", "", $4); split($4,a,"/"); if (a[2] != "France") print $1"\t"$2"\t"$3"\t"a[2]}' ${globalSeqs%.fasta}.tsv | grep -v EPI_ISL_402125) \
    --targetN 100 \
    --outName data/random_global_seqs

# add reference
cat data/random_global_seqs.fasta data/EPI_ISL_402125.fasta > data/random_global_seqs_ref.fasta

# get the metadata for these sequences
awk -F'\t' 'FNR==NR {hash[$1]; next} $2 in hash' \
    <(grep ">"  data/random_global_seqs_ref.fasta |\
    sed 's/\>//g' | awk -F'|' '{print $2}') \
    ${globalSeqs%.fasta}.tsv \
    | awk -F'\t'  '{gsub (" ", "", $4); split($4,a,"/");  print $1"\t"$2"\t"$3"\t"a[2]}' \
    > data/random_global_seqs_ref.tsv

# align the downsampled sequences
python3 scripts/align_seqs.py \
    --sequences data/random_global_seqs_ref.fasta \
    --referenceName EPI_ISL_402125 \
    --minLength 2800 \
    --maskHead 55 \
    --maskTail 100 \
    --maskSites 11083 15324 21575 \
    --alignType reference

# add to france sequences
cat ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref.fasta data/random_global_seqs_ref_aligned_ref_filtered_masked.fasta \
    > data/france_random_global_ref_aligned.fasta

# concatenate the metadata
cat <(tail -n +2 ${FranceSeqs%.fasta}.tsv | awk -F '\t' '{gsub (" ", "", $4); split($4,a,"/");  print $1"\t"$2"\t"$3"\t"a[2]}') data/random_global_seqs_ref.tsv \
    > data/france_random_global_ref.tsv

# get list of sequences used
grep ">" data/france_random_global_ref_aligned.fasta \
    | cut -d'|' -f 2 > data/france_random_global_ref_ids.tsv


# --------------------------------------------------------------------------------------#
# //// 7. BUILD ML AND TIME ALIGNED TREE  ////
# --------------------------------------------------------------------------------------#
iqtree2 -s data/france_random_global_ref_aligned.fasta  --polytomy -T 4

python3 scripts/refine_tree.py \
    --sequences data/france_random_global_ref_aligned.fasta \
    --tree data/france_random_global_ref_aligned.fasta.treefile \
    --resolve_polytomies True \
    --metadata data/france_random_global_ref.tsv \
    --root_name EPI_ISL_402125 \
    --metadataNameCol 1 \
    --metadataDateCol 2  \
    --metadataSep "\t"

cat \
    <(awk -F'|' '{print $2"\tincluded"}' ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest.tsv) \
    <(awk -F'\t' 'FNR==NR {hash[$1]; next} !($1 in hash) {print $0"\tnot_included_france"}' \
        <(awk -F'|' '{print $2"\tincluded"}' ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest.tsv) \
        <(tail -n +2 ${FranceSeqs%.fasta}.tsv | awk -F'\t' '{print $2}')) \
    <(awk -F'\t' '{print $2"\tnot_included_exog"}' data/random_global_seqs_ref.tsv) \
    > data/france_random_global_ref_figure.tsv


# --------------------------------------------------------------------------------------#
# //// 9. GET EARLIEST BASAL GENOTYPE ////
# --------------------------------------------------------------------------------------#
# get all sequences which have a given set of SNPs
# MAY HAVE OTHER SNPS AS WELL
python3 scripts/get_match_seqs.py \
    --SNPs <(awk -F'\t' -v largeset=$largest '{if ($1==largest) print $0}' ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_snp.tsv | cut -f2-) \
    --seqs ${globalSeqs%.fasta}_aligned_ref_filtered_masked.fasta \
    --ref $ReferenceSeq \
    --nucDict scripts/nuc_dict_all.json \
    --outName globalSeqs_match_basal

# --------------------------------------------------------------------------------------#
# //// 14. PLOT RESULTS ////
# --------------------------------------------------------------------------------------#
python3 scripts/plot_tree.py \
    --metadata <(awk -F'\t' 'FNR==NR{hash[$2]=$3; next}{print $0"\t"hash[$1]}' \
        ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_match_largest_seqs.tsv \
        data/france_random_global_ref_figure.tsv) \
    --trees data/france_random_global_ref_aligned_refine.nwk data/france_random_global_ref_aligned_refine_time.nwk \
    --outName data/france_trees \
    --distDat ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_EPI_ISL_402125.tsv


python3 scripts/plot_empirical_data.py \
    --sDat ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s.tsv \
    --metadata ${FranceSeqs%.fasta}_date_country.csv \
    --outName ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s



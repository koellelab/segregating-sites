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
globalSeqs=data/gisaid_hcov-19_2021_06_12_01.fasta
globalSeqs2=data/gisaid_hcov-19_2020_03_31_complete_hc_date_EHC_metadata_aligned_ref_filtered_masked.fasta


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

# get metadata just for those sequences
awk -F'\t' 'FNR==NR {hash[$1]; next} $2 in hash' \
    <(grep ">" ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs.fasta |\
    sed 's/\>//g' | awk -F'|' '{print $2}') \
    ${FranceSeqs%.fasta}.tsv \
    > ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_match_largest_seqs.tsv

# get list of sequences used
grep ">" ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs.fasta \
    | cut -d'|' -f 2 > data/gisaid_hcov-19_2021_04_29_16_match_largest_seqs_ids.tsv


# --------------------------------------------------------------------------------------#
# //// 4. CALCULATE SEGREGATING SITES OF FRANCE DATA ////
# --------------------------------------------------------------------------------------#
# calculate segregating sites 
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
# need to downsample 100 random global sequences
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


python3 scripts/refine_tree.py \
    --sequences data/gtr_test/france_random_global_ref_aligned.fasta \
    --tree data/gtr_test/france_random_global_ref_aligned.fasta.treefile \
    --resolve_polytomies True \
    --metadata data/france_random_global_ref.tsv \
    --root_name EPI_ISL_402125 \
    --metadataNameCol 1 \
    --metadataDateCol 2  \
    --metadataSep "\t"


# --------------------------------------------------------------------------------------#
# //// 9. GET EARLIEST BASAL GENOTYPE ////
# --------------------------------------------------------------------------------------#
python3 scripts/get_match_seqs.py \
    --SNPs <(awk -F'\t' -v largeset=$largest '{if ($1==largest) print $0}' ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_snp.tsv | cut -f2-) \
    --seqs $globalSeqs2 \
    --ref $ReferenseSeq \
    --nucDict scripts/nuc_dict_all.json \
    


clades.py \
    --ccSNPs ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_snp.tsv \
    --cc ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc.tsv \
    --dists ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref.dist \
    --toMatch $largest \
    --maxDist  0 \
    --nucDict scripts/nuc_dict_all.json \
    --ref $ReferenceSeq 



# --------------------------------------------------------------------------------------#
# //// 14. PLOT RESULTS ////
# --------------------------------------------------------------------------------------#
python3 scripts/plot_empirical_data.py \
# --------------------------------------------------------------------------------------#
# //// 9. GET DATA FOR MU ESTIMATION ////
# --------------------------------------------------------------------------------------#
# data files for estimation of mu
# all from original studies 
# will send if you email mmart59@emory.edu and ask
# but don't want to republish data that doesn't belong to us
popaData1='data/scitranslmed.abe2555_Data_file_S1.tsv'
popaData4='data/scitranslmed.abe2555_Data_file_S4.tsv'
braunData2='data/biorxiv.2021.04.30.440988_s2.tsv'
braunData4='data/biorxiv.2021.04.30.440988_S4.tsv'
jamesData3='data/biorxiv.2020.11.15.20231993_3.tsv'
jamesData1='data/biorxiv.2020.11.15.20231993_S1.tsv'
# for privacy reasons we cannot release these data
lythgoeData1='data/science.abg082_household_samples.csv'
# gisaid nextmeta file, just to associates name w/ accession
lythgoeData2='data/science.abg082_metadata.tsv'
# file name for downloaded sequence data
seqData='data/gisaid_hcov-19_2021_06_22_20.fasta'

# format and concatenate data
awk -F'\t' 'FNR==NR { map[$14] = $13; next } { $2 = map[$2]; $10 = map[$10]; print $2"\t"$10}' \
    $popaData1 $popaData4 | tail -n +2 > data/scitranslmed.abe2555_pair_accessions.tsv

awk -F'\t' 'FNR==NR { map[$1] = $6; next } {  $1 = map[$1]; $2 = map[$2]; print $1"\t"$2}' \
    $braunData2 $braunData4 > data/biorxiv.2021.04.30.440988_accessions.tsv

awk -F'\t' 'FNR==NR { gsub(/\r/,"",$12); map[$12] = $4; next } {$1 = map[$1]; $2 = map[$2]; print $1"\t"$2}' \
    <(grep "CH1" $jamesData1) \
    <(awk -F'\t' '{print $1}' $jamesData3 | sed $'s/ to /\t/g' | tail -n +3) \
    > data/biorxiv.2020.11.15.20231993_accessions.tsv

# for the lythgoe data, we have sequence names, not gisaid accession numbers, so need to convert
tail -n +2 $lythgoeData1 | \
    awk -F',' '{split($10,a,"_"); split($11,b,"_"); if (length(a) > 1) print a[1]"\t"b[1]}' \
    > data/science.abg082_pair_names.tsv

awk -F'\t' '{print $1"\n"$2}' data/science.abg082_pair_names.tsv \
    > data/science.abg082_all_names.tsv

awk -F'\t' 'FNR==NR {split($1,a,"/"); hash[a[2]] = $2; next} print hash[$1]}' \
    $lythgoeData2 \
    data/science.abg082_pair_names.tsv

awk -F'\t' 'FNR==NR {split($1, a, "/"); hash[a[2]]=$3; next} {if ($1 in hash && $2 in hash) print hash[$1]"\t"hash[$2]}' \
    <(grep "United Kingdom" $lythgoeData2 | grep 2020) \
    data/science.abg082_pair_names.tsv \
    > data/science.abg082_pair_accessions.tsv

cat \
    <(awk -F'\t' '{print "popa\t"$0}' data/scitranslmed.abe2555_pair_accessions.tsv) \
    <(awk -F'\t' '{print "braun\t"$0}' data/biorxiv.2021.04.30.440988_accessions.tsv) \
    <(awk -F'\t' '{print "james\t"$0}' data/biorxiv.2020.11.15.20231993_accessions.tsv) \
    <(awk -F'\t' '{print "lythgoe\t"$0}' data/science.abg082_pair_accessions.tsv) \
    > data/combined_metadata.tsv

# get a list of all accessions
awk -F'\t' '{print $2"\n"$3}' data/combined_metadata.tsv | \
    sort | uniq \
    > data/all_accessions.tsv

# download accessions from gisaid

# combine with ref
cat $seqData $refSeq > ${seqData%.fasta}_ref.fasta

# --------------------------------------------------------------------------------------#
# //// 10. ALIGN SEQUENCES FOR MU ESTIMATION ////
# --------------------------------------------------------------------------------------#
# align accessions to Wuhan/Hu-1 using MAFFT
python3 scripts/align_seqs.py \
    --sequences  ${seqData%.fasta}_ref.fasta \
    --maskHead 55 \
    --maskTail 100 \
    --alignType reference \
    --referenceName EPI_ISL_402125

# --------------------------------------------------------------------------------------#
# //// 11. GET PAIRWISE DISTANCES ////
# --------------------------------------------------------------------------------------#

python3 scripts/pair_distances.py \
    --seqs ${seqData%.fasta}_ref_aligned_ref_filtered_masked.fasta \
    --pairDat data/combined_metadata.tsv \
    --nucDict scripts/nuc_dict_all.json


# --------------------------------------------------------------------------------------#
# //// 12. FORMAT DATA TABLE ////
# --------------------------------------------------------------------------------------#
cat \
    <(echo "Study\tDonor\tRecipient\t# SNPs") \
    <(sed 's/popa/Popa et al. 2020/g' data/combined_metadata_dist.tsv | \
        sed 's/braun/Braun et al. 2021/g' | \
        sed 's/james/James et al. 2020/g' | \
        sed 's/lythgoe/Lythgoe et al. 2021/g' | \
        awk -F'\t' '{if ($1 == "Lythgoe et al. 2021") print $1"\tNA\tNA\t"$4; else print $0}') \
    > data/combined_metadata_dist_format.tsv


# --------------------------------------------------------------------------------------#
# //// 13. ESTIMATE MU ////
# --------------------------------------------------------------------------------------#
python3 scripts/calc_mu.py \
    --distDat data/combined_metadata_dist.tsv


    --sDat ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s.tsv \
    --muDat data/combined_metadata_dist.tsv \
    --muEst data/combined_metadata_dist_mle.tsv \
    --metadata ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_match_largest_seqs.tsv \
    --outName figures/france_s_dat

python3 scripts/plot_tree.py \
    --metadata <(awk -F'\t' 'FNR==NR{hash[$2]=$3; next}{print $0"\t"hash[$1]}' \
        ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_match_largest_seqs.tsv \
        data/france_random_global_ref_figure.tsv) \
    --trees data/france_random_global_ref_aligned_refine.nwk data/france_random_global_ref_aligned_refine_time.nwk \
    --outName figures/france_trees \
    --distDat ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_EPI_ISL_402125.tsv


head ${FranceSeqs%.fasta}_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_EPI_ISL_402125.tsv




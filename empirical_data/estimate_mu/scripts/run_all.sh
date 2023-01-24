
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


# --------------------------------------------------------------------------------------#
# //// 1. ORGANIZE DATA ////
# --------------------------------------------------------------------------------------#
# Popa
awk -F'\t' 'FNR==NR { map[$14] = $13; next } { $2 = map[$2]; $10 = map[$10]; print $2"\t"$10}' \
    $popaData1 $popaData4 | tail -n +2 > data/scitranslmed.abe2555_pair_accessions.tsv
# Braun
awk -F'\t' 'FNR==NR { map[$1] = $6; next } {  $1 = map[$1]; $2 = map[$2]; print $1"\t"$2}' \
    $braunData2 $braunData4 > data/biorxiv.2021.04.30.440988_accessions.tsv
# James
awk -F'\t' 'FNR==NR { gsub(/\r/,"",$12); map[$12] = $4; next } {$1 = map[$1]; $2 = map[$2]; print $1"\t"$2}' \
    <(grep "CH1" $jamesData1) \
    <(awk -F'\t' '{print $1}' $jamesData3 | sed $'s/ to /\t/g' | tail -n +3) \
    > data/biorxiv.2020.11.15.20231993_accessions.tsv

# Lythgoe
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

# get list of all sequences
awk -F'\t' '{print $2"\n"$3}' data/combined_metadata_dist_format.tsv | grep -v 'NA' > all_seqs.tsv
# --------------------------------------------------------------------------------------#
# //// 2. ALIGN SEQUENCES FOR MU ESTIMATION ////
# --------------------------------------------------------------------------------------#
# align accessions to Wuhan/Hu-1 using MAFFT
python3 scripts/align_seqs.py \
    --sequences  ${seqData%.fasta}_ref.fasta \
    --maskHead 55 \
    --maskTail 100 \
    --alignType reference \
    --referenceName EPI_ISL_402125

# --------------------------------------------------------------------------------------#
# //// 3. GET PAIRWISE DISTANCES ////
# --------------------------------------------------------------------------------------#

python3 scripts/pair_distances.py \
    --seqs ${seqData%.fasta}_ref_aligned_ref_filtered_masked.fasta \
    --pairDat data/combined_metadata.tsv \
    --nucDict scripts/nuc_dict_all.json


# --------------------------------------------------------------------------------------#
# //// 4. FORMAT DATA TABLE ////
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
# //// 5. ESTIMATE MU ////
# --------------------------------------------------------------------------------------#
python3 scripts/calc_mu.py \
    --distDat data/combined_metadata_dist_format.tsv





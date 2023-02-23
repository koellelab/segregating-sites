FranceSeqs=data/gisaid_hcov-19_2021_04_29_16.fasta


python3 scripts/plot_empirical_data.py \
    --sDat france_data/data/gisaid_hcov-19_2021_04_29_16_ref_aligned_ref_filtered_masked_noref_cc_match_largest_seqs_s.tsv \
    --muDat estimate_mu/data/combined_metadata_dist_format.tsv \
    --muEst estimate_mu/data/combined_metadata_dist_format_mle.tsv

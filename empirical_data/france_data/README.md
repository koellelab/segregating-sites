# Summary
This workflow analyzes sequences data from France presented in Fig 6 and S7 of [Epidemiological inference for emerging viruses using segregating sites](https://doi.org/10.1101/2021.07.07.451508). 
These scripts were primarily written by [Michael Martin](https://github.com/m-a-martin) who can be reached on GitHub or at [mmart108@jhmi.edu](mailto:mmart108@jhmi.edu). 

## Data
These analyses assume that the following data files are present in the `./data` directory. Unfortunately, due to GISAID rules we cannot share these data directly, but list of accession numbers for each file have been provided. 

1. `data/gisaid_hcov-19_2021_04_29_16.fasta` This data file includes complete, high coverage, collection date complete seqeucnes from France sampled thorugh march 17th uploaded through April 29th 2021. 
2. 1. `data/gisaid_hcov-19_2021_04_29_16.tsv` Metadata file associated with France sequence data. 
3. `data/gisaid_hcov-19_2021_06_12_01.fasta` This data file includes complete, high coverage, collection date complete sequences sampled through march 17th uploaded through April 29th 2021. 
4. `data/gisaid_hcov-19_2021_06_12_01.fasta` Metadata file associated with global sequence data. 
5. `data/EPI_ISL_402125.fasta` Wuhan/Hu-1 reference genome

## Workflow
This workflow completes 8 main steps, many of which include sub-steps. 

1. Align France sequence data. 
	1. Add country column to GISAID metadata using awk.
	2. Get list of IDs which are used, these are provided on the GitHub instead of the FASTA due to GISAID sharing rules. 
	3. Concatenate France sequence data and reference sequence. 
	4. Finally, align France sequences to the reference sequence masking first 55 nucleotides, last 100 nucleotides, and site 11083, 15324, and 21575 using the `scripts/align_seqs.py` script. 
2. Get connected components of the France data. 
	1. Calcualte the pairwise distances between all France sequences using the `scripts/calc_pairwise_distances.py` file. 
	2. Get the connected components using a cutoff of 1 SNP based on the pairwise distances using the `scripts/get_connected_components.py` file. 
3. Get sequences that are part of or descend from the largest connected component. 
	1. Identify the largest connected component using `awk`. 
	2. Get all connected components which share the SNPs present in the largest connected components and match the reference sequence at all other positions using the `scripts/get_match_clades.py` file. Excludes any connected components which differ from the largest by >7 SNPs. 
	3. Get the sequences present in the connected components identifed in #2 using teh `scripts/get_seqs.py` file and their associated metadata using `awk`. 
	4. Get a list of the sequences in #3, this list is provided on GitHub instead of the FASTA due to the GISAID sharing rules. 
4. Use the `scripts/calc_segregating_sites.py` file to calcualte the segregating sites trajectory for the France data identified in #3 using a time window of 4 days. 
5. Use the `scripts/calc_distance.py` file to calculate the nucleotide distance between each sequence identified in #3 and Wuhan/Hu-1. 
6. Incorporate global sequence data for phylogenetic tree plot. 
	1. Downsample the global seqs file to 100 random sequences using the `scripts/downsample_seqs.py` script. 
	2. Add the reference to those downsampled sequences with `cat`. 
	3. Get metadata for those downsampled sequences using `awk`. 
	4. Align downsampled global sequences to reference sequence using the `scripts/align_seqs.py` script masking the first 55 and last 100 nucleotides and sites 11083, 15324, 21575. 
	5. Combine downsampled global sequences with France sequences using `cat`. 
	6. Get a list of combined sequences which is provided on GitHub instead of the FASTA due to GISAID sharing rules. 
7. Build ML divergence and time aligned trees. 
	1. Use IQtree2 to infer a ML phylogeny of the combined France and global sequences, identifying the best fit nucleotide substituiton model with ModelTest and collapsing near zero branches. 
	2. Use the `scripts/refine_tree.py` script which internally uses TreeTime to infer a time aligned ML tree rooting at EPI_ISL_402125. 
8. Plot the results. 
	1. Using the `scripts/plot_empirical_data.py` script plot the segregating site trajectory of the largest CC and its descendants in France. 
	2. Using the `scripts/plot_tree.py` script plot the phylogenetic trees of all France sequences and the downsampled global sequences.  

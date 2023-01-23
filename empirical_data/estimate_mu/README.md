# Summary
This workflow analyzes consensus sequence from epidemiologically linked transmission pairs to estimate the P(mutation | transmission) as shown in Fig 6C of [Epidemiological inference for emerging viruses using segregating sites](https://doi.org/10.1101/2021.07.07.451508). 
These scripts were primarily written by [Michael Martin](https://github.com/m-a-martin) who can be reached on GitHub or at [mmart108@jhmi.edu](mailto:mmart108@jhmi.edu). 

## Data
This analysis combines data from four studies: 
1. Data File S1 and Data File S4 of [Genomic epidemiology of superspreading events in Austria reveals mutational dynamics and transmission properties of SARS-CoV-2](https://doi.org/10.1126/scitranslmed.abe2555). 
2. Supplementary Table 2 and Supplementary Table 4 of [Acute SARS-CoV-2 infections harbor limited within-host diversity and transmit via tight transmission bottlenecks](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1009849#sec025). 
3. Table 3 and Supplementary Table 1 of [Transmission dynamics of SARS-CoV-2 within-host diversity in two major hospital outbreaks in South Africa](https://academic.oup.com/ve/article/7/1/veab041/6248115?login=false#supplementary-data). 
4. Unreleased source data from [SARS-CoV-2 within-host diversity and transmission](https://www.science.org/doi/10.1126/science.abg0821)which was provided by the authors. Due to privacy reasons we cannot share these data. 

## Workflow
This workflow completes the following steps: 
1. Accession numbers involved in all transmission pairs for each study are organized and combined into a single metadata file and the corresponding sequences downloaded from GISAID. 
2. Sequences are aligned to Wuhan/Hu-1 using MAFFT and masking the first 55 and last 100 nucleotides using the `scripts/align_seqs.py` script. 
3. Nucleotide distances between all transmission pairs are calculated using the `scripts/pair_distances.py` script. 
4. Distances are used to calculate the average number of mutations per transmission assuiming a poisson distribution using the `scripts/calc_mu.py` file. 
5. The observed and estimated probability distributions are plotted using the `scripts/plot_mu.py` script. 


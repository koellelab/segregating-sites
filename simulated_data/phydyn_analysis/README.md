# Summary
This workflow runs the BEAST2 PhyDyn analysis presented in Fig 3B, 5B, and 5D of [Epidemiological inference for emerging viruses using segregating sites](https://doi.org/10.1101/2021.07.07.451508). 
These scripts were primarily written by [Michael Martin](https://github.com/m-a-martin) who can be reached on GitHub or at [mmart108@jhmi.edu](mailto:mmart108@jhmi.edu). 


These analyses are conducted on two sets of simulated data from the same simulated epidemic (`data/SEIR_simple_simtaj.npz`) with different sampling strategies. The first sampling strategy (`data/SEIR_simple_prop500_0928_seed123.npz`) includes 500 viral sequences randomly sampled from throughout the simulated epidemic. The second sampling strategy (`data/SEIR_seed1234_unif10_sampled_during_32-52_220428.npz`) includes 50 sequences sampled between day 32 and day 52 of the simulated epidemic. Finally, we simulate a third epidemic (`data/SEIR_[mu=4e-1]_seed221110/SEIR_[mu=4e-1]_seed221110_simtaj.npz`) using a slightly higher substitution rate and sample sequences between day 32 and 52 (`data/SEIR_[mu=4e-1]_seed221110/SEIR_[mu=4e-1]_seed221110_unif10_32-52.npz`).

For each, the following steps are conducted:

1. `scripts/generate_fasta.py` is used to simulate a fasta alignment from the simulated mutation sparse array assuming a JC69 substitution model. 
2. `scripts/generate_xml.py` is ued to generate a BEAST XML file from the simulated alignment. 
3. BEAST2 with PhyDyn is run on the BEAST XML. 
4. `python3 scripts/density_plot.py` is used to plot the results. 


Precise code to replicate each step is available in `scripts/run_all.sh`.
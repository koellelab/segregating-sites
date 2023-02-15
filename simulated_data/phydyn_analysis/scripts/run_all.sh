R0=1.6

#### 1. 500 SAMPLED PROPORTIONALLY THROUGH EPIDEMIC ####
mu=0.2
# generate fasta
python3 scripts/generate_fasta.py \
	--simulationData data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_simtaj.npz \
	--sampledData data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_sampling.npz

# make metadata file 
grep ">" data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_sampling.fasta | \
	sed 's/\>//g' | \
	awk '{split($1, a, "_"); print $1"\t"a[5]/366}' > \
	data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_sampling_metadata.tsv

# clock rate is in units of subs/site
# but size of alignment changes with shorter simulations
# so, we need to convert
# set max of ucld mean to:
# 10 * subs/transmission * R0 transmissions/infection * 366 days/yr / 5 days per transmission
first_seq=$(tail -n +2 "data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_sampling.fasta" | head -n 1)
genome_size=${#first_seq}

python3 scripts/generate_xml.py \
	--seqs data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_sampling.fasta \
	--alnName simpleSEIR_R0=16e-1_mu2e-1_seed201234_prop500_win4 \
	--metadata data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/seed230201_sampling_metadata.tsv \
	--ucldMeanMax $(((10 * mu * R0 * 366 /5)/genome_size))


# run beast2
beast2 -threads 4  -working data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_prop500_win4/simpleSEIR_R0=16e-1_mu2e-1_seed201234_prop500_win4.xml

# generate figure
python3 scripts/density_plot.py \
	--metadata data/simpleSEIR_R0=24e-1_mu2e-1/seed1234_prop500_win2/seed230127_sampling_metadata.tsv \
	--logFile data/simpleSEIR_R0=24e-1_mu2e-1/seed1234_prop500_win2/simpleSEIR_R0=24e-1_mu2e-1_seed20230127_prop500_win2.log

#### 2. SAMPLED UNIFORMLY BETWEEN DAYS 32 AND 52 ####
mu=0.2
# generate fasta
python3 scripts/generate_fasta.py \
	--simulationData data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_simtaj.npz \
	--sampledData data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif10_win4_32-52/seed230201_sampling.npz

# make metadata file 
grep ">" data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif10_win4_32-52/seed230201_sampling.fasta  | \
	sed 's/\>//g' | \
	awk '{split($1, a, "_"); print $1"\t"a[5]/366}' > \
	data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif10_win4_32-52/seed230201_sampling.tsv


# clock rate is in units of subs/site
# but size of alignment changes with shorter simulations
# so, we need to convert
# set max of ucld mean to:
# 10 * subs/transmission * R0 transmissions/infection * 366 days/yr / 5 days per transmission
first_seq=$(tail -n +2 "data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif10_win4_32-52/seed230201_sampling.fasta" | head -n 1)
genome_size=${#first_seq}

python3 scripts/generate_xml.py \
	--seqs data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif10_win4_32-52/seed230201_sampling.fasta \
	--alnName simpleSEIR_R0=16e-1_mu2e-1_seed1234_unif10_win4_32-52 \
	--ucldMeanMax $(((10 * mu * R0 * 366 /5)/genome_size)) \
	--metadata data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif10_win4_32-52/seed230201_sampling.tsv

# run beast
beast2 -threads 4 -working data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif10_win4_32-52/simpleSEIR_R0=16e-1_mu2e-1_seed1234_unif10_win4_32-52.xml
# make figure
python3 scripts/density_plot.py \
	--metadata data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif10_win4_32-52/seed230201_sampling.tsv \
	--logFile data/simpleSEIR_R0=16e-1_mu2e-1/seed1234_unif10_win4_32-52/simpleSEIR_R0=16e-1_mu2e-1_seed1234_unif10_win4_32-52.log

#### 3. SAMPLED UNIFORMLY BETWEEN DAYS 32 AND 52, MU=0.4 ####
mu=0.4
# generate fasta
python3 scripts/generate_fasta.py \
	--simulationData "data/simpleSEIR_R0=16e-1_mu4e-1/seed221110_simtaj.npz" \
	--sampledData "data/simpleSEIR_R0=16e-1_mu4e-1/seed221110_unif10_win4_32-52/seed230201_sampling.npz"


# make metadata file 
grep ">" data/simpleSEIR_R0=16e-1_mu4e-1/seed221110_unif10_win4_32-52/seed230201_sampling.fasta | \
	sed 's/\>//g' | \
	awk '{split($1, a, "_"); print $1"\t"a[5]/366}' > \
	data/simpleSEIR_R0=16e-1_mu4e-1/seed221110_unif10_win4_32-52/seed230201_sampling.tsv


# clock rate is in units of subs/site
# but size of alignment changes with shorter simulations
# so, we need to convert
# set max of ucld mean to:
# 10 * subs/transmission * R0 transmissions/infection * 366 days/yr / 5 days per transmission
first_seq=$(tail -n +2 "data/simpleSEIR_R0=16e-1_mu4e-1/seed221110_unif10_win4_32-52/seed230201_sampling.fasta" | head -n 1)
genome_size=${#first_seq}

python3 scripts/generate_xml.py \
	--seqs "data/simpleSEIR_R0=16e-1_mu4e-1/seed221110_unif10_win4_32-52/seed230201_sampling.fasta" \
	--alnName "simpleSEIR_R0=16e-1_mu4e-1_seed221110_unif10_win4_32-52" \
	--ucldMeanMax $(((10 * mu * R0 * 366 /5)/genome_size)) \
	--metadata data/simpleSEIR_R0=16e-1_mu4e-1/seed221110_unif10_win4_32-52/seed230201_sampling.tsv

# run beast
beast2 -threads 4 -working data/simpleSEIR_R0=16e-1_mu4e-1/seed221110_unif10_win4_32-52/simpleSEIR_R0=16e-1_mu4e-1_seed221110_unif10_win4_32-52.xml


# generate figure
python3 scripts/density_plot.py \
	--metadata data/simpleSEIR_R0=16e-1_mu4e-1/seed20230201_unif10_win4_32-52/seed20230201_sampling.tsv \
	--logFile data/simpleSEIR_R0=16e-1_mu4e-1/seed20230201_unif10_win4_32-52/simpleSEIR_R0=16e-1_mu4e-1_seed20230201_unif10_win4_32-52.log





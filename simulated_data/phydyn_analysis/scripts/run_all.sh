#### 1. 500 SAMPLED PROPORTIONALLY THROUGH EPIDEMIC ####
# generate fasta
python3 generate_fasta.py \
	--simulationData data/SEIR_simple_simtaj.npz \
	--sampledData data/SEIR_simple_prop500_0928_seed123.npz

# make metadata file 
grep ">" data/SEIR_simple_prop500_0928_seed123.fasta | \
	sed 's/\>//g' | \
	awk '{split($1, a, "_"); print $1"\t"a[5]/366}' > \
	data/SEIR_simple_prop500_0928_seed123.tsv

# generate beast xml
python3 generate_xml.py \
	--seqs data/SEIR_simple_prop500_0928_seed123.fasta \
	--alnName data/SEIR_simple_prop500_0928_seed123 \
	--metadata data/SEIR_simple_prop500_0928_seed123.tsv

# run beast2
beast2 -threads 4  data/SEIR_simple_prop500_0928_seed123.xml

# generate figure
python3 scripts/density_plot.py \
	--metadata data/SEIR_simple_prop500_0928_seed123.tsv \
	--logFile data/SEIR_simple_prop500_0928_seed123.log

#### 2. SAMPLED UNIFORMLY BETWEEN DAYS 32 AND 52 ####
# generate fasta
python3 generate_fasta.py \
	--simulationData SEIR_simple_simtaj.npz \
	--sampledData SEIR_seed1234_unif10_sampled_during_32-52_220428.npz

# make metadata file 
grep ">" SEIR_seed1234_unif10_sampled_during_32-52_220428.fasta | \
	sed 's/\>//g' | \
	awk '{split($1, a, "_"); print $1"\t"a[5]/366}' > \
	SEIR_seed1234_unif10_sampled_during_32-52_220428.tsv

# make xml
# clock rate is in units of subs/site
# but size of alignment changes with shorter simulations
# so, we need to convert
# assume 0.2 subs/transmission
# assume R0 = 1.6 trans/5 days
# know 365 days/year
# thus 0.2*1.6*366/5 = 23.424 substitution/year
#know genome size of 125
# 23.424/125 = 0.19 subs/site/year
# set max of prior to double this 
python3 generate_xml.py \
	--seqs SEIR_seed1234_unif10_sampled_during_32-52_220428.fasta \
	--alnName SEIR_seed1234_unif10_sampled_during_32-52_220428 \
	--ucldMeanMax 0.38

# run beast
beast2 -threads 4  data/SEIR_simple_prop500_0928_seed123.xml

# make figure
python3 scripts/density_plot.py \
	--metadata <(
		grep ">" "data/SEIR_seed1234_unif10_sampled_during_32-52_220428.fasta" | \
			awk '{split($1, a, "_"); print $0"\t"a[5]/366}' | \
			sed 's/\>//g') \
	--logFile data/SEIR_seed1234_unif10_sampled_during_32-52_220428.log


#### 3. SAMPLED UNIFORMLY BETWEEN DAYS 32 AND 52, MU=0.4 ####
# generate fasta
python3 scripts/generate_fasta.py \
	--simulationData "data/SEIR_[mu=4e-1]_seed221110/SEIR_[mu=4e-1]_seed221110_simtaj.npz" \
	--sampledData "data/SEIR_[mu=4e-1]_seed221110/SEIR_[mu=4e-1]_seed221110_unif10_32-52.npz"

# generate xml
# clock rate is in units of subs/site
# but size of alignment changes with shorter simulations
# so, we need to convert
# assume 0.4 subs/transmission
# assume R0 = 1.6 trans/5 days
# know 365 days/year
# thus 0.4*1.6*366/5 = 46.848 substitution/year
# know genome size of 79
# 46.848/79 = 0.59 subs/site/year
# set max of prior to double this 
python3 scripts/generate_xml.py \
	--seqs "data/SEIR_[mu=4e-1]_seed221110/SEIR_[mu=4e-1]_seed221110_unif10_32-52.fasta" \
	--alnName "SEIR_mu04_seed221110_unif10_32-52" \
	--ucldMeanMax 1.19

# run beast
beast2 -threads 4  data/SEIR_[mu=4e-1]_seed221110/SEIR_simple_prop500_0928_seed123.xml

# generate figure
python3 scripts/density_plot.py \
	--metadata <(
		grep ">" "data/SEIR_[mu=4e-1]_seed221110/SEIR_[mu=4e-1]_seed221110_unif10_32-52.fasta" | \
			awk '{split($1, a, "_"); print $0"\t"a[5]/366}' | \
			sed 's/\>//g') \
	--logFile "data/SEIR_[mu=4e-1]_seed221110/SEIR_mu04_seed221110_unif10_32-52.log"





# first need to generate fasta alignments

python3 generate_fasta.py \
	--simulationData SEIR_simple_simtaj.npz \
	--sampledData SEIR_simple_prop500_0928_seed123.npz


python3 generate_xml.py \
	--seqs SEIR_simple_prop500_0928_seed123.fasta \
	--alnName SEIR_simple_prop500_0928_seed123 \
	--metadata SEIR_simple_prop500_0928_seed123.tsv

python3 generate_fasta.py \
	--simulationData SEIR_simple_simtaj.npz \
	--sampledData SEIR_seed1234_unif10_sampled_during_32-52_220428.npz

# make metadata file 

 
grep ">" SEIR_seed1234_unif10_sampled_during_32-52_220428.fasta | \
	sed 's/\>//g' | \
	awk '{split($1, a, "_"); print $1"\t"a[5]/366}' > \
	SEIR_seed1234_unif10_sampled_during_32-52_220428.tsv
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



python3 generate_fasta.py \
	--simulationData SEIR_simple_simtaj.npz \
	--sampledData SEIR_seed1234_unif10_sampled_during_32-60_220430.npz
	
# clock rate is in units of subs/site
# but size of alignment changes with shorter simulations
# so, we need to convert
# assume 0.2 subs/transmission
# assume 1 trans/5 days
# know 365 days/year
# thus 0.2*1.6*366/5 = 23.424 substitution/year
#know genome size of 219
# 23.424/219 = 0.11 subs/site/year
# set max of prior to double this 
python3 generate_xml.py \
	--seqs SEIR_seed1234_unif10_sampled_during_32-60_220430.fasta \
	--alnName SEIR_seed1234_unif10_sampled_during_32-60_220430 \
	--ucldMeanMax 0.22

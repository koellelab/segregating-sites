import argparse
import numpy as np


def generate_sampled_fasta(sampled_individuals, genotype_info):
	genotype_info = \
		[genotype_info[:,0][genotype_info[:,1][i]:genotype_info[:,1][i+1]] 
			for i in np.arange(genotype_info[:,1].max())]
	sampled_genotype_info = \
		[genotype_info[i]
			for i in sampled_individuals[:,1].astype(int)]
	mutation_locations = np.unique(np.concatenate(sampled_genotype_info))
	# generate random reference geome
	# asume all frequencies and all transition probabilities are equal
	aln_length = mutation_locations.max() + 1
	rng = np.random.default_rng()
	nucleotides = np.array(['A', 'C', 'T', 'G'])
	ref_genome = rng.choice(nucleotides, size=aln_length)
	# for each polymorphic site, pick an alt nucleotide that is different form the ref nucleotide
	mutation_ref_nucs = ref_genome[mutation_locations]
	mutation_alt_nucs = \
		np.tile(nucleotides,mutation_locations.shape[0]).\
			reshape((mutation_locations.shape[0], nucleotides.shape[0]))
	nucleotides_idx = {'A': 0, 'C': 1, 'T': 2, 'G': 3}
	u,inv = np.unique(mutation_ref_nucs,return_inverse = True)
	# dummy character so we can sort and remove
	mutation_alt_nucs[np.arange(mutation_alt_nucs.shape[0]), 
		np.array([nucleotides_idx[x] for x in u])[inv].reshape(mutation_ref_nucs.shape)] = 'Z'
	# sorts each row and removes last column (the dummy Z)
	mutation_alt_nucs = np.sort(mutation_alt_nucs)[:,:3]
	alt_nucs = np.zeros(aln_length, dtype='<U1')
	alt_nucs[mutation_locations] = \
		mutation_alt_nucs[np.arange(mutation_alt_nucs.shape[0]), rng.choice(3, size=mutation_alt_nucs.shape[0])]
	# now generate aln
	aln = \
		np.tile(ref_genome,sampled_individuals.shape[0]).reshape(
			(sampled_individuals.shape[0], ref_genome.shape[0]))
	aln_names = []
	for sampled_individual_idx, sampled_individual in enumerate(sampled_individuals):
		aln[sampled_individual_idx, sampled_genotype_info[sampled_individual_idx]] = \
			alt_nucs[ sampled_genotype_info[sampled_individual_idx]]
		aln_names.append(f'sample_{sampled_individual_idx}_genotype_{sampled_individual[1].astype(int)}_{sampled_individual[0]}')
	return(aln_names, aln)


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--sampledData')
	parser.add_argument('--simulationData')
	args = parser.parse_args()
	#args.simulationData = 'SEIR_simple_simtaj.npz'
	#args.sampledData = 'SEIR_seed1234_unif10_sampled_during_32-52_220428.npz'
	#args.sampledData = 'SEIR_simple_prop500_0928_seed123.npz'
	simulation_data = np.load(args.simulationData)
		

	#args.sampledData = 'SEIR_simple_prop500_0928_seed123.npz'
	sampled_data = np.load(args.sampledData)

	sampled_individuals = sampled_data['sampled_individuals']
	genotype_info = simulation_data['genotype_info']
	aln_names, aln = generate_sampled_fasta(sampled_individuals, genotype_info)
	with open('.'.join(args.sampledData.split('.')[0:-1])+'.fasta', 'w') as fp:
		for idx, i in enumerate(aln_names):
			_ = fp.write('>'+i+'\n'+''.join(aln[idx,:])+'\n')




if __name__ == "__main__":
	run()



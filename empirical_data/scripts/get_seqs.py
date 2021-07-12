#!/usr/bin/env python3
import pandas as pd
from Bio import SeqIO
import argparse


def get_seqs(seqs=None, sequenceField=1, sequenceSep='|',
  getSeqs=None, getSeqsCol=0, getSeqsSep='|', getSeqsField=1, getSeqsDelim='\t',
  exclude=False, formatNames=None, outName='filtered_seqs'):
	print(getSeqsCol)
	get_seqs = set(pd.read_csv(getSeqs, sep=getSeqsDelim, header=None)[getSeqsCol].str.split(getSeqsSep, expand=True)[getSeqsField])
	orig_seqs = list(SeqIO.parse(seqs, 'fasta'))
	if formatNames:
		get_seqs = set([i.replace(formatNames[0], formatNames[1]) for i in get_seqs])
		for i in orig_seqs:
			i.description = i.description.replace(formatNames[0], formatNames[1])
			i.name = i.description
			i.id = i.description
	if exclude == True:
		print(len(orig_seqs))
		filtered_seqs = \
			[item for item in orig_seqs if 
				item.description.split(sequenceSep)[sequenceField] not in get_seqs]
		print(len(filtered_seqs))
	else:
		print(len(orig_seqs))
		filtered_seqs = \
			[item for item in orig_seqs if 
				item.description.split(sequenceSep)[sequenceField] in get_seqs]
		print(f'{len(filtered_seqs)} of the desired sequences found in alignment')
	with open(f'{outName}.fasta', 'w') as out:
		SeqIO.write(filtered_seqs, out, 'fasta')


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--outName', default='seqs',
	                    help='name of output files')
	parser.add_argument('--seqs', default=None,
	                    help='fasta file with sequences')
	parser.add_argument('--sequenceField', default=1, type=int,
						help='field in fasta header with the sequence ID to match to input file')
	parser.add_argument('--sequenceSep', default='|', 
						help='seperator for fasta header')
	parser.add_argument('--getSeqs', default=None,
		help='list of sequence IDs to get from fasta (assume no header, csv)')
	parser.add_argument('--getSeqsCol', default=0, type=int,
		help='column in sequence list')
	parser.add_argument('--getSeqsSep', default='|',
		help='how to split sequence names')
	parser.add_argument('--getSeqsField', default=1, type=int,
		help='column in sequence list')
	parser.add_argument('--getSeqsDelim', default='\t')
	parser.add_argument('--excludeSeqs', 
		help='whether to exlude sequences in list instead of include them',
		action='store_true', dest="exclude")
	parser.add_argument('--formatNames', 
		help='string to replace in names',
		nargs=2)
	parser.set_defaults(exclude=False)
	args = parser.parse_args()
	get_seqs(**args.__dict__)


if __name__ == "__main__":
	run()


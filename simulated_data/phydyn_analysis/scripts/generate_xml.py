#!/usr/bin/env python3
from Bio import SeqIO
import pandas as pd
import argparse



def numeric_from_datetime(dt):
    from calendar import isleap
    import datetime
    if dt is None:
        dt = datetime.datetime.now()
    days = 366 if isleap(dt.year) else 365
    # removing the 0.5 adjustment in the standard treetime code 
    # to align with tajimas D inference method
    res =  dt.year + (dt.timetuple().tm_yday) / days
    return(res)


def generate_sequence_block(sequence_list):
    sequence_block = ''
    for item in sequence_list:
        sequence_block += '<sequence id="seq_{0}" spec="Sequence" taxon="{0}" \
        totalcount="4" value="{1}" />\n'.format(item.description, item.seq)
    return(sequence_block)


def generate_time_block(sequence_list, date_dict):
    time_block = ''
    for item in sequence_list:
        time_block += f'{item.description}={date_dict[item.description]},'
    # removes trailing ','
    time_block = ''.join(list(time_block)[0:-1])
    return(time_block)



def generate_xml():
    parser = argparse.ArgumentParser()
    # input files
    parser.add_argument('--xmlTemplate', 
        help='path to xml template file', 
        default='config/template.xml')
    parser.add_argument('--seqs', 
        help='path to aln')
    parser.add_argument('--seqNameSep',
        help='character to split sequence name on to match metadata',
        default='|')
    parser.add_argument('--seqNameField',
        help='field of split sequence name to match metadata',
        default=0, type=int)
    parser.add_argument('--seqNameAppend',
        default='_I')
    parser.add_argument('--metadata',
        help='path to metadata')
    parser.add_argument('--metadataDelim',
        help='metadata file delimiter',
        default='\t')
    parser.add_argument('--metadataNameCol',
        help='column in metadata with sequence names',
        type=int,
        default=0)
    parser.add_argument('--metadataDateCol',
        help='column in metadata with date',
        type=int,
        default=1)
    parser.add_argument('--metadataDateFmt',
        help='date format',
        default='numeric')
    parser.add_argument('--alnName',
        help='name of discrete trait',
        default='sars_cov_2')
    parser.add_argument('--ucldMeanMax',
        default=1E-2, type=float)
    args = parser.parse_args()
    print(args.ucldMeanMax)
    # default values used for testing
    #args.seqs = 'SEIR_simple_prop60_sampled_during_early45_1218.fasta'
    #args.metadata = 'SEIR_simple_prop150_sampled_during_28-56_1219.tsv'
    #args.alnName = 'SEIR_simple_prop60_sampled_during_early45_1218'


    seqs = list(SeqIO.parse(args.seqs, 'fasta'))
    n_seqs = len(seqs)
    print(f'{n_seqs} sequences in input fasta')


    if args.metadata:
        metadata = \
            pd.read_csv(args.metadata, sep=args.metadataDelim, header=None)
        metadata_names = set(metadata[args.metadataNameCol])
        seqs = [i for i in seqs if i.description.split(args.seqNameSep)[args.seqNameField] in metadata_names]
        for seq in seqs:
            seq.description = seq.description + args.seqNameAppend
        print(f'{len(seqs)} sequences found in metadata')
        # creates a new column with trait labels 
        # build date dict
        if args.metadataDateFmt != 'numeric':
            metadata['numeric_date'] = \
                pd.to_datetime(metadata[args.metadataDateCol], format=args.metadataDateFmt).apply(numeric_from_datetime)
        else:
            metadata['numeric_date'] = \
                metadata[args.metadataDateCol]
        if args.seqNameAppend:
            metadata[args.metadataNameCol] = metadata[args.metadataNameCol].apply(lambda k: k+args.seqNameAppend)
        date_dict = \
            {row[args.metadataNameCol]: row['numeric_date'] for idx, row in metadata.iterrows()}

    else:
        for seq in seqs:
            seq.description = seq.description + args.seqNameAppend
        date_dict = {i.description: float(i.name.split('_')[-1])/366 for i in seqs}


    sequence_block = generate_sequence_block(seqs)
    time_block = generate_time_block(seqs, date_dict)
    n_dimensions = str(int(2*len(seqs) - 2))
    with open(args.xmlTemplate, 'r') as infile:
            template = infile.read()

            
    print(args.alnName)
    template = template.replace('<!-- ALN_NAME -->', args.alnName)
    template = template.replace('<!-- SEQUENCE_BLOCK -->', sequence_block)
    template = template.replace('<!-- TIME_BLOCK -->', time_block)
    template = template.replace('<!-- N_DIMENSIONS -->', n_dimensions)
    template = template.replace('<!-- UCLD_MEAN_MAX -->', str(args.ucldMeanMax))
    template = template.replace('<!-- UCLD_MEAN_INIT -->', str(args.ucldMeanMax*0.2))
    with open('/'.join(args.seqs.split('/')[:-1])+'/' + args.alnName+'.xml', 'w') as outfile:
        outfile.write(template)


if __name__ == '__main__':
    generate_xml()

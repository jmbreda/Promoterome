import pandas as pd
import numpy as np
import h5py
import re
import sys
sys.path.insert(0, '/home/jbreda/Genome_sequence')
import genome_sequence
import argparse


def parse_argument():
    parser = argparse.ArgumentParser(description='Get promoter plus window around tss.')
    parser.add_argument('--genome'
                        ,required=True
                        ,type=str
                        ,choices=['mm10','mm39','hg19','hg38']
                        ,help="Genome organism and version")
    parser.add_argument('--infile'
                        ,required=True
                        ,type=str
                        ,help="Input gff table with promoter")
    parser.add_argument('--outfile_hdf5'
                        ,required=True
                        ,type=str
                        ,help="Output datastructure with promoter window sequences per promoter")

    return parser.parse_args()

if __name__ == '__main__':

    args = parse_argument()

    # load promoter table
    promoterome = pd.read_csv(args.infile,sep='\t')

    # replace chrN by N for compatibility with genome_sequence
    promoterome.loc[:,'chr'] = promoterome.chr.apply( lambda c: re.sub('chr','',c) )

    N_prom = promoterome.shape[0]
    win_size = promoterome.at[0,'end'] - promoterome.at[0,'start']
    N_nuc = 4
    # save promoter one hot encoded sequence in hdf5
    Prom_seq = np.zeros((N_prom,win_size,N_nuc),dtype=np.uint8)
    for p in promoterome.index:
        coord = promoterome.loc[p,['chr','start','end','strand']].values
        Prom_seq[p,:,:] = genome_sequence.OneHotEncoding(genome_sequence.SamCoord_2_Seq(args.genome,coord))

    with h5py.File(args.outfile_hdf5, 'w') as hf:
        prom = hf.create_dataset('sequence',data=Prom_seq)
        prom.attrs['genome'] = args.genome
        prom.attrs['window_size'] = win_size
        prom.attrs['promoterome_bed_file'] = args.infile

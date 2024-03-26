import pandas as pd
import os
import argparse

def parse_argument():
    parser = argparse.ArgumentParser(description='Get promoter plus window around tss.')
    parser.add_argument('--genome'
                        ,required=True
                        ,type=str
                        ,choices=['mm10','hg38']
                        ,help="Genome organism and version")
    parser.add_argument('--window_kb'
                        ,default=2
                        ,type=int
                        ,help="window size in kb")
    parser.add_argument('--promoterome'
                        ,required=True
                        ,type=str
                        ,help="Input gff table with promoter")
    parser.add_argument('--blacklist'
                        ,required=True
                        ,type=str
                        ,help="Input bed file with blacklist")
    parser.add_argument('--outfile_table'
                        ,required=True
                        ,type=str
                        ,help="Output table with promoter windows")

    return parser.parse_args()

if __name__ == '__main__':

    args = parse_argument()

    # load promoter table
    #if args.genome == 'mm10':
    #    infile_promoterome = '/bigdata/jbreda/genome/mm10/mm10_promoters_v2.gff.gz'
    #elif args.genome == 'hg38':
    #    infile_promoterome = '/bigdata/jbreda/genome/hg38/hg38_promoters_v1.gff.gz'

    # blacklist
    #infile_blacklist = f'../Blacklist/lists/{args.genome}-blacklist.v2.bed.gz'

    # intersect promoterome with blacklist
    promoterome_blacklist_annot = f'resources/{args.genome}/promoterome_blacklist_annotation.gff'
    os.system(f'bedtools intersect -a {args.promoterome} -b {args.blacklist} -loj > {promoterome_blacklist_annot}')
    
    # read promoterome with blacklist annotation
    promoterome = pd.read_csv(promoterome_blacklist_annot,sep='\t',header=None,comment='#',usecols=[0,3,4,6,8,12])
    promoterome.columns = ['chr','start','end','strand','name','black_listed']

    # remove 'M' chr
    promoterome.drop(index=promoterome[promoterome.chr=='chrM'].index,inplace=True)

    # get gene name and id (chr_strand_start_end)
    promoterome.loc[:,'gene'] = promoterome.name.apply( lambda p: p.split('|')[1] )
    if args.genome == 'mm10':
        promoterome.loc[:,'id'] = promoterome.name.apply( lambda p: p.split('|')[0].split(';')[0].replace('ID=mm10_v2_','') )
    elif args.genome == 'hg38':
        promoterome.loc[:,'id'] = promoterome.name.apply( lambda p: p.split('|')[0].split(';')[0].replace('ID=hg38_v1_','') )
    promoterome.drop(columns='name',inplace=True)
    promoterome.reset_index(drop=True,inplace=True)

    # get position of promoter
    middle = ((promoterome.loc[:,'end'] + promoterome.loc[:,'start'])/2).astype(int)
    promoterome.loc[:,'start'] = middle - args.window_kb*1000
    promoterome.loc[ promoterome.start<0,'start'] = 0
    promoterome.loc[:,'end'] = middle + args.window_kb*1000

    promoterome.loc[promoterome.black_listed != 'High Signal Region',:].to_csv(args.outfile_table,sep='\t',index=False)

            

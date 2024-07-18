import numpy as np
import pandas as pd
import argparse

def parse_argument():
    parser = argparse.ArgumentParser(description='Cluster gene with close by promoters.')
    parser.add_argument('--infile_promoterome'
                        ,required=True
                        ,type=str
                        ,help="Input bed table with promoter")
    parser.add_argument('--threshold'
                        ,default=200
                        ,type=int
                        ,help="distance threshold in bp")
    parser.add_argument('--outfile_promoterome'
                        ,required=True
                        ,type=str
                        ,help="Output table with clustered promoter windows")

    return parser.parse_args()

def find_clusters(pos,idx, threshold):
    
    assert np.all(pos == sorted(pos)), "Points must be sorted"
    n = len(pos)

    # Initialize variables
    groups = []
    current_group = []
    
    # Iterate through the sorted points
    for i in range(n):
        if not current_group:
            # Start the first group with the first point
            current_group.append(idx[i])
        else:
            # Check if the current point can be added to the current group
            if pos[i] - pos[i-1] <= threshold:
                current_group.append(idx[i])
            else:
                # If not, finalize the current group and start a new one
                groups.append(current_group)
                current_group = [idx[i]]
    # Finally, add the last group
    if current_group:
        groups.append(current_group)

    return groups


if __name__ == '__main__':

    args = parse_argument()

    Promoterome = pd.read_csv(args.infile_promoterome, sep='\t')
    CHR = ['chr' + str(i) for i in range(1, 20)] + ['chrX', 'chrY']
    STRAND = ['+', '-']

    # get the index of promoter clusters
    groups = []
    for chr in CHR:
        for strand in STRAND:
            # find duplicated genes names in promoterome
            Promoterome_chr = Promoterome[(Promoterome['chr'] == chr) & (Promoterome['strand'] == strand)]
            duplicated_genes = Promoterome_chr[Promoterome_chr.duplicated(subset='gene', keep=False)]

            if np.any( duplicated_genes.loc[:,'start'].values[1:] - duplicated_genes.loc[:,'start'].values[:-1] <= 0 ):
                print('Warning: duplicated genes are not sorted by start position')
                print(chr,strand)
                break

            for g in duplicated_genes.gene.unique():
                # find duplicated genes indices
                idx = Promoterome_chr[Promoterome_chr['gene'] == g].index
                pos = np.sum(Promoterome_chr.loc[idx, ['start','end']].values,1) / 2
                groups.extend(find_clusters(pos,idx,args.threshold))

    # merge promoters in each group
    for group in groups:
        if len(group) > 1:
            Promoterome.loc[group,'gene'].values

            # merge all promoters in group to the 1st one
            start = int(np.round(Promoterome.loc[group,'start'].mean()))
            end = int(np.round(Promoterome.loc[group,'end'].mean()))
            id = '|'.join( Promoterome.loc[group,'id'].values )
            Promoterome.loc[group[0],'start'] = start
            Promoterome.loc[group[0],'end'] = end
            Promoterome.loc[group[0],'id'] = id
            # drop the other promoters
            Promoterome.drop(group[1:], inplace=True)

    # check if windows are unique
    assert len(np.unique(np.diff(Promoterome.loc[:,['start','end']].values,1))), "Promoterome windows are not unique"

    # save the clustered promoterome
    Promoterome.to_csv(args.outfile_promoterome, sep='\t', index=False)


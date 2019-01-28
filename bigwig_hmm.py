"""
This is a wrapper script to run HMMs (pomegranate or hmmlearn)
with a few bells and whistles
v1.1
"""
import bbi
import numpy as np
import pandas as pd
import bioframe
from pomegranate import *

def get_chroms(genome='hg38', ignoreXYMT=True):
    "Get list of chroms to analyze"
    chromsizes = bioframe.fetch_chromsizes(genome)
    chr_list = list(chromsizes.index)
    if ignoreXYMT == True:
        chr_list = [i for i in chr_list if i not in ('chrM', 'chrY','chrY')]
    return chr_list

def create_df(filepath,chroms):
    "Create dataframe from bigwig"
    df=pd.DataFrame(columns=['chrom','start','end','value'])
    for item in chroms:
        ivals = list(bbi.fetch_intervals(filepath,item,0,-1))
        df_new=pd.DataFrame(ivals, columns=['chrom','start','end','value'])
        df=df.append(df_new, ignore_index=True)
    return df

def hmm(df,num_states=2):
    "HMM program"
    model = HiddenMarkovModel.from_samples(NormalDistribution,X=[df['value'].values], n_components=num_states)
    states=model.viterbi(df['value'].values)
    listofstates = [i[0] for i in states[1]]
    listofstates.pop(0) #first value isn't a state
    df['state']=listofstates
    return df

def sparse(df):
    "Merge neighboring bins with same state"
    state_list=[]
    chr_list=[]
    start_list=[]
    end_list=[]
    state_list.append(df['state'].iloc[0])
    chr_list.append(df['chrom'].iloc[0])
    start_list.append(df['start'].iloc[0])
    previous_state=df['state'].iloc[0]

    for index, row in df.iloc[1:].iterrows():
        if index != 0 and df['start'].iloc[index] == 0:
            end_list.append(df['end'].iloc[(index-1)])
            chr_list.append(df['chrom'].iloc[index])
            start_list.append(df['start'].iloc[index])
            state_list.append(df['state'].iloc[index])
            continue
        if df['state'].iloc[index] == df['state'].iloc[(index-1)]:
            continue
        else:
            end_list.append(df['end'].iloc[(index-1)])
            chr_list.append(df['chrom'].iloc[index])
            start_list.append(df['start'].iloc[index])
            state_list.append(df['state'].iloc[index])
    end_list.append(df['end'].iloc[(index)])

    keys=['chrom','start','end','state']
    values=[chr_list,start_list,end_list,state_list]
    dictionary = dict(zip(keys, values))
    df_sparse=pd.DataFrame.from_dict(dictionary)
    return df_sparse #<-need to sort this out with write to file

def sort_states(df):
    "Re-sort states with LAD regions labelled State 0"
    df = df[(df['chrom'] == 'chr11') & (df['start'] >= 40000000) & (df['start'] >= 42000000)]

def write_to_file(df,filename,num_states=2,):
    df['score'] = '0'
    df['strand'] = '.'
    filename=filename+str(num_states)+'_state_HMM.bed'
    df.loc[df['state'] == 0, 'RGB'] = '0,255,0'
    df.loc[df['state'] == 1, 'RGB'] = '0,255,255'
    df.loc[df['state'] == 2, 'RGB'] = '255,0,0'
    cols_to_keep=['chrom','start','end','state','score','strand','start','end','RGB']
    df.to_csv(filename, sep='\t', header=False, columns=cols_to_keep, index=False)
    return print('All Finished!')

def main():
    chroms=get_chroms()
    print(chroms)
    df=create_df(filepath='DKO_ladseq.log2ratio.bw',chroms=chroms)
    print('made df!')
    print(df.head())
    df=hmm(df)
    print('finished hmm!')
    print(df.head())
    df_sparse=sparse(df)
    print('finished sparse_df!')
    print(df_sparse.head())
    write_to_file(df_sparse,filename='test')
    print("Finished writing to file!!")


if __name__ == '__main__':
    main()

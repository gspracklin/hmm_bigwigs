"""
This is a wrapper script to run HMMs (pomegranate or hmmlearn)
with a few bells and whistles
v1.3
"""
#!/bin/env python
import argparse
import bbi
import numpy as np
import pandas as pd
import bioframe
from pomegranate import HiddenMarkovModel, NormalDistribution
import matplotlib.pyplot as plt


def get_chroms(genome, ignoreXYMT=True):
    "Get list of chroms to analyze"
    print("Using chroms from " + genome)
    chromsizes = bioframe.fetch_chromsizes(genome, filter_chroms=True)
    chr_list = list(chromsizes.index)
    if ignoreXYMT == True:
        chr_list = [i for i in chr_list if i not in ("chrM", "chrX", "chrY")]
    return chr_list


def create_df(inputfile, chroms):
    "Create dataframe from bigwig"
    df = pd.DataFrame(columns=["chrom", "start", "end", "value"])
    for item in chroms:
        ivals = list(bbi.fetch_intervals(inputfile, item, 0, -1))
        df_new = pd.DataFrame(ivals, columns=["chrom", "start", "end", "value"])
        df = df.append(df_new, ignore_index=True)
    return df


def hmm(df, num_states):
    "HMM program"
    # df['value']=df['value'].replace(0,np.nan) #this removes unmappable areas of chr
    # df_dropna=df.dropna(subset=['value']) #this removes unmappable areas of chr (NaN is otherwise considered 0)
    vals = df["value"].values
    model = HiddenMarkovModel.from_samples(
        NormalDistribution, X=[vals], n_components=num_states
    )
    states = model.predict(vals)

    # Rename states to increase with mean signal
    order = np.argsort(df['value'].groupby(states).mean())
    states = [order[s] for s in states]
    df["state"] = states
    df['state'][np.isnan(df['value'])] = np.nan
    return df


def sparse(df):
    "Merge neighboring bins with same state"
    chr_list = []
    start_list = []
    state_list = []
    end_list = []

    for item in df["chrom"].unique():
        chrom_df = df[df["chrom"] == item].reset_index()

        chr_list.append((chrom_df["chrom"].iloc[0]))
        start_list.append((chrom_df["start"].iloc[0]))
        state_list.append((chrom_df["state"].iloc[0]))
        for index, row in chrom_df[1:].iterrows():
            if chrom_df["state"].iloc[index] == chrom_df["state"].iloc[(index - 1)]:
                continue
            else:
                end_list.append(chrom_df["end"].iloc[(index - 1)])
                chr_list.append(chrom_df["chrom"].iloc[index])
                start_list.append(chrom_df["start"].iloc[index])
                state_list.append(chrom_df["state"].iloc[index])
        if len(start_list) != len(end_list):
            end_list.append(chrom_df["end"].iloc[(index)])

    keys = ["chrom", "start", "end", "state"]
    values = [chr_list, start_list, end_list, state_list]
    dictionary = dict(zip(keys, values))
    df_sparse = pd.DataFrame.from_dict(dictionary)
    return df_sparse.dropna()


def merge_different_hmmstates(df, cLAD, open):
    "merge strong and weak HMM states into 2 "
    import pandas as pd

    chr_list = []
    start_list = []
    end_list = []
    weak = int(3 - (cLAD + open))

    for item in df["chrom"].unique():
        chrom_df = df[df["chrom"] == item]
        start = 1
        for index, row in chrom_df.iterrows():
            if start == 1:
                if df["state"].iloc[index] == cLAD or df["state"].iloc[index] == weak:
                    chr_list.append(df["chrom"].iloc[index])
                    start_list.append(df["start"].iloc[index])
                    start = 0
                    continue
                else:
                    continue
            elif df["state"].iloc[index] == open:
                end_list.append(df["end"].iloc[(index - 1)])
                start = 1
            else:
                continue
        if start == 0:
            end_list.append(df["end"].iloc[(index)])

        if len(chr_list) != len(start_list) or len(start_list) != len(end_list):
            print("Wrong Lengths!")
            break

    keys = ["chrom", "start", "end", "state"]
    values = [chr_list, start_list, end_list]
    dictionary = dict(zip(keys, values))
    df_merge = pd.DataFrame.from_dict(dictionary)
    return df_merge

def write_to_file(df, outputfile, num_states, cmap='coolwarm'):
    states = list(range(num_states))
    cmap = plt.get_cmap(cmap)
    colors = {s:cmap(s/states[-1]) for s in states}
    df["score"] = "0"
    df["strand"] = "."
    filename = outputfile + "_" + str(num_states) + "_state_HMM_colored.bed"
    df['RGB'] = df["state"].apply(lambda x: ','.join([str(int(round(c*255))) for c in colors[x][:-1]]))
    cols_to_keep = [
        "chrom",
        "start",
        "end",
        "state",
        "score",
        "strand",
        "start",
        "start",
        "RGB",
    ]
    df.to_csv(filename, sep="\t", header=False, columns=cols_to_keep, index=False)

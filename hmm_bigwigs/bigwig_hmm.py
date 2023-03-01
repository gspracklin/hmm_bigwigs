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


def create_df(inputfile, view):
    "Create dataframe from bigwig"
    df_list = []
    for i, row in view.iterrows():
        ivals = bbi.fetch_intervals(inputfile, row.chrom, row.start, row.end)
        df_new = pd.DataFrame(ivals, columns=["chrom", "start", "end", "value"])
        df_list.append(df_new)
    return pd.concat(df_list).reset_index(drop=True)


def hmm(df, num_states):
    "HMM program"
    # df["value"] = df["value"].replace(0, np.nan)  # this removes unmappable areas of chr
    df = df.dropna(
        subset=["value"]
    )  # this removes unmappable areas of chr (NaN is otherwise considered 0)
    if df.shape[0] == 0:
        df = df.assign(state=None)
        return df
    vals = df["value"].values
    model = HiddenMarkovModel.from_samples(
        NormalDistribution, X=[vals], n_components=num_states
    )
    states = model.predict(vals)

    # Rename states to increase with mean signal
    order = np.argsort(df["value"].groupby(states).mean())
    states = [order[s] for s in states]
    df["state"] = states
    df["state"][np.isnan(df["value"])] = np.nan
    return df


def sparse(df):
    "Merge neighboring bins with same state"
    df_sparse = bioframe.merge(df, on=["state"]).dropna()
    return df_sparse


def merge_different_hmmstates(df, cLAD, open):
    "merge strong and weak HMM states into 2"
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


def write_to_file(df, outputfile=None, cmap="coolwarm"):
    states = list(sorted(df["state"].unique()))
    cmap = plt.get_cmap(cmap)
    colors = {s: cmap(s / states[-1]) for s in states}
    df["score"] = "0"
    df["strand"] = "."
    df["RGB"] = df["state"].apply(
        lambda x: ",".join([str(int(round(c * 255))) for c in colors[x][:-1]])
    )
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
    if outputfile:
        df.to_csv(outputfile, sep="\t", header=False, columns=cols_to_keep, index=False)
    else:
        print(
            df.to_csv(
                outputfile, sep="\t", header=False, columns=cols_to_keep, index=False
            )
        )

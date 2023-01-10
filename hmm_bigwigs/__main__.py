#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from hmm_bigwigs import *


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create bedfile of HMM states from bigwig"
    )
    parser.add_argument(
        "-i",
        "--inputfile",
        type=str,
        dest="inputfile",
        help="path to the bigwig file(s)",
        action="store",
        required=True,
    )
    parser.add_argument(
        "-g", "--genome", dest="genome", help="genome (i.e. hg38)", action="store"
    )
    parser.add_argument(
        "-n",
        "--num_states",
        type=int,
        dest="num_states",
        help="number of hidden states",
        action="store",
        default=2,
    )
    parser.add_argument(
        "--cmap",
        type=str,
        dest="cmap",
        help="Colormap to map states to colors",
        action="store",
        required=False,
        default="coolwarm",
    )
    parser.add_argument(
        "-o",
        "--outputfile",
        dest="outputfile",
        help="output filename (endswith:#_state_HMM.bed)",
        action="store",
        required=False,
    )
    parser.add_argument(
        "-save-split-files",
        dest="savesplit",
        help="Whether to save separate bed files split by state in addition to the output file",
        action="store",
        required=False,
    )
    args = parser.parse_args()
    if args.outputfile is None:
        args.outputfile = args.inputfile
    return args


def main():
    args = parse_args()
    # print("Starting HMM on " + args.inputfile)
    chroms = get_chroms(args.genome)
    df = create_df(inputfile=args.inputfile, chroms=chroms)
    df = hmm(df, args.num_states)
    # print("Finished hmm!")
    df_sparse = sparse(df)
    write_to_file(df_sparse, args.outputfile, args.num_states, cmap=args.cmap)
    # df_final=merge_different_hmmstates(df_sparse, cLAD=cLAD, open=open_state)
    # df_final.to_csv(args.outputfile+'_combined_state.bed', sep='\t', header=False, index=False)
    # print("write first file")
    if args.savesplit:
        df_sparse[df_sparse["state"] == 0].to_csv(
            args.outputfile + "_0_state.bed", sep="\t", header=False, index=False
        )
        df_sparse[df_sparse["state"] == 1].to_csv(
            args.outputfile + "_1_state.bed", sep="\t", header=False, index=False
        )
        df_sparse[df_sparse["state"] == 2].to_csv(
            args.outputfile + "_2_state.bed", sep="\t", header=False, index=False
        )
    # print("Finished writing to file")

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from natsort import natsort_keygen
import bbi
import bioframe as bf
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
        "--view", dest="view_file", help="file with genomic view", action="store"
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
        "--save-split-files",
        dest="savesplit",
        help="Whether to save separate bed files split by state in addition to the output file",
        action="store_true",
        required=False,
        default=False,
    )
    args = parser.parse_args()
    if args.outputfile is None:
        args.outputfile = args.inputfile
    return args


def main():
    args = parse_args()
    if args.view_file is not None:
        view = bf.make_viewframe(
            bf.read_table(args.view_file, schema="bed4", index_col=False)
        )
    elif args.genome is not None:
        view = bf.make_viewframe(bf.fetch_chromsizes(args.genome, as_bed=True))
    else:
        with bbi.open(args.inputfile) as f:
            view = (
                bf.make_viewframe(dict(f.chromsizes))
                .sort_values("chrom", key=natsort_keygen())
                .reset_index(drop=True)
            )
    df = create_df(inputfile=args.inputfile, view=view)
    df = bf.assign_view(df, view).dropna(subset="value")
    try:
        df = df.groupby("view_region").apply(hmm, num_states=args.num_states)
    except ValueError:
        print(df)
    df_sparse = sparse(df)
    write_to_file(df_sparse, args.outputfile, cmap=args.cmap)
    if args.savesplit:
        for state in range(args.num_states):
            df_sparse[df_sparse["state"] == state].to_csv(
                f"{args.outputfile}_{state}_state.bed",
                sep="\t",
                header=False,
                index=False,
            )

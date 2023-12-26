#!/usr/bin/env python3

import argparse
import pandas as pd
from functools import reduce


def parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--samples", nargs="+", type=str, help="List of samples", required=True
    )
    args = parser.parse_args()
    return args


def merge_dfs(dfs: list) -> pd.DataFrame:
    df = reduce(lambda left, right: pd.merge(left, right, on="Name"), dfs)
    df.set_index("Name", inplace=True)
    return df


def main():
    args = parser()
    dfs = []

    for sample in args.samples:
        name = sample.split("/")[-1].split(".genes")[0]
        df = pd.read_csv(sample, sep="\t", usecols=[0, 4])
        df.rename(columns={"NumReads": name}, inplace=True)
        dfs.append(df)

    mtx = merge_dfs(dfs)
    mtx.to_csv("counts.txt", sep="\t")

    return None


if __name__ == "__main__":
    main()

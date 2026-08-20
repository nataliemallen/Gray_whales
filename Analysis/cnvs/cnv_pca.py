#!/usr/bin/env python3

import os
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# paths

results_dir = "/scratch/gautschi/allen715/2026_whales/cnvs/results"

pop1_samples = [
    l.strip()
    for l in open("/scratch/gautschi/allen715/2026_whales/lists/east.txt")
]

pop2_samples = [
    l.strip()
    for l in open("/scratch/gautschi/allen715/2026_whales/lists/west_ds.txt")
]

all_samples = pop1_samples + pop2_samples

colors = {
    "ENP": "#3FA5D5",
    "WNP": "#E7AF23",
}


# functions

def parse_file(file_path, sample):
    df = pd.read_csv(file_path, sep="\t", header=None)

    df.columns = [
        "type",
        "region",
        "length",
        "rd",
        "p1",
        "p2",
        "p3",
        "p4",
        "a",
        "b",
        "c",
    ]

    df[["scaffold", "coords"]] = df["region"].str.split(":", expand=True)
    df[["start", "end"]] = df["coords"].str.split("-", expand=True)

    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["sample"] = sample

    return df[["type", "scaffold", "start", "end", "sample"]]


def run_pca(bin_size):

    print(f"\nProcessing {bin_size}...")

    calls = []

    for sample in all_samples:
        file_path = os.path.join(
            results_dir,
            f"{sample}.{bin_size}.calls.tsv"
        )

        if not os.path.exists(file_path):
            raise FileNotFoundError(file_path)

        calls.append(parse_file(file_path, sample))

    cnv = pd.concat(calls, ignore_index=True)

    # Define unique CNV regions
    regions = (
        cnv.groupby(["type", "scaffold", "start", "end"])
        .size()
        .reset_index(name="count")
    )
    regions["id"] = regions.index

    sample_region = []

    for _, r in regions.iterrows():

        mask = (
            (cnv["type"] == r["type"])
            & (cnv["scaffold"] == r["scaffold"])
            & (cnv["start"] == r["start"])
            & (cnv["end"] == r["end"])
        )

        for sample in cnv.loc[mask, "sample"].unique():
            sample_region.append(
                {
                    "sample": sample,
                    "region": r["id"],
                }
            )

    sr = pd.DataFrame(sample_region)

    matrix = pd.crosstab(sr["sample"], sr["region"])

    pca = PCA(n_components=2)
    coords = pca.fit_transform(matrix)

    explained = pca.explained_variance_ratio_ * 100

    print(
        f"  PC1: {explained[0]:.2f}%   "
        f"PC2: {explained[1]:.2f}%"
    )

    plt.figure(figsize=(8, 6))

    for samples, label in [
        (pop1_samples, "ENP"),
        (pop2_samples, "WNP"),
    ]:

        idx = matrix.index.isin(samples)

        plt.scatter(
            coords[idx, 0],
            coords[idx, 1],
            color=colors[label],
            label=label,
            s=40,
        )

    plt.xlabel(f"PC1 ({explained[0]:.2f}%)", fontsize=16)
    plt.ylabel(f"PC2 ({explained[1]:.2f}%)", fontsize=16)

    # Increase tick label size
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    # Legend with border (boxed key)
    plt.legend(
        frameon=True,
        facecolor="white",
        edgecolor="black",
        fontsize=14
    )

    plt.tight_layout()

    outfile = f"cnv_pca_{bin_size}.png"
    plt.savefig(outfile, dpi=300)
    plt.close()


# run both bin sizes

for bin_size in ["10kb", "100kb"]:
    run_pca(bin_size)

print("\nDone.")
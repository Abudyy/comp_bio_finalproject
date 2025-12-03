import pandas as pd
from intervaltree import IntervalTree

#Load CTCF motif positions
ctcf = pd.read_csv(
    r"C:\Users\saleh\561_project\data\CTCF_sites.bed",
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "tf", "score", "strand"]
)

# Load GM12878 regulatory regions
gm = pd.read_csv(
    r"C:\Users\saleh\561_project\data\GM12878_only.bed",
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "tf_reg", "score_reg", "cells"]
)

#Dictionary of interval trees for each chromosome
trees = {}

print("started Building interval trees:")

for chrom in gm["chrom"].unique():
    gm_chr = gm[gm["chrom"] == chrom]

    tree = IntervalTree()
    for _, row in gm_chr.iterrows():
        tree[row["start"]:row["end"]] = True

    trees[chrom] = tree

print("Interval trees built.")

#CTCF overlap
positive_rows = []

print("fiinding overlaps:")

for idx, row in ctcf.iterrows():
    chrom = row["chrom"]
    if chrom not in trees:
        continue

    start = int(row["start"])
    end = int(row["end"])

    # intervaltree lookup (FAST)
    if trees[chrom].overlaps(start, end):
        positive_rows.append(row)

positives = pd.DataFrame(positive_rows)

print("Positive CTCF sites in GM12878 regulatory regions:", len(positives))

positives.to_csv(
    r"C:\Users\saleh\561_project\data\CTCF_positive.bed",
    sep="\t",
    header=False,
    index=False
)

print("Saved:", r"C:\Users\saleh\561_project\data\CTCF_positive.bed")

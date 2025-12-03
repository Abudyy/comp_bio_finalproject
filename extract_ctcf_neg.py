import pandas as pd
from intervaltree import IntervalTree
import random



#All CTCF motifs (used to EXCLUDE from negatives)
ctcf_all = pd.read_csv(
    r"C:\Users\saleh\561_project\data\CTCF_sites.bed",
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "tf", "score", "strand"]
)

# GM12878 regulatory regions
gm = pd.read_csv(
    r"C:\Users\saleh\561_project\data\GM12878_only.bed",
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "tf_reg", "score_reg", "cells"]
)

# Number of positives (we'll match this)
n_positives = len(pd.read_csv(
    r"C:\Users\saleh\561_project\data\CTCF_positive.bed",
    sep="\t", header=None))

print("Will generate this many negatives:", n_positives)

# 

# motif length = typical CTCF window length (first row)
motif_length = int(ctcf_all["end"].iloc[0] - ctcf_all["start"].iloc[0])
print("Using motif length:", motif_length)

# 

trees = {}
print("Building exclusion interval trees for CTCF motifs...")

for chrom in ctcf_all["chrom"].unique():
    ctcf_chr = ctcf_all[ctcf_all["chrom"] == chrom]
    tree = IntervalTree()

    for _, row in ctcf_chr.iterrows():
        tree[row["start"]:row["end"]] = True

    trees[chrom] = tree

print("Trees built.")


# ---------- SAMPLE NEGATIVES ----------

negatives = []
attempts = 0
MAX_ATTEMPTS = 5_000_000  #limit

print("Sampling negatives...")

while len(negatives) < n_positives and attempts < MAX_ATTEMPTS:

    # pick a random GM12878 region
    row = gm.sample(1).iloc[0]
    chrom = row["chrom"]

    if chrom not in trees:
        attempts += 1
        continue

    region_start = int(row["start"])
    region_end   = int(row["end"])

    # region too small?
    if region_end - region_start <= motif_length:
        attempts += 1
        continue

   
    s = random.randint(region_start, region_end - motif_length)
    e = s + motif_length

   
    if trees[chrom].overlaps(s, e):
        attempts += 1
        continue

    # valid negative example
    negatives.append([chrom, s, e, "CTCF", 0.0, "+"])
    attempts += 1

    if len(negatives) % 5000 == 0:
        print("Negatives collected:", len(negatives))



neg_df = pd.DataFrame(negatives)
neg_df.to_csv(
    r"C:\Users\saleh\561_project\data\CTCF_negative.bed",
    sep="\t",
    header=False,
    index=False
)

print("DONE. Total negatives:", len(neg_df))
print("Saved to CTCF_negative.bed")

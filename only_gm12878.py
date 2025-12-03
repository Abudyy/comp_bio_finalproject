import pandas as pd

df = pd.read_csv(
    r"C:\Users\saleh\561_project\data\wgEncodeRegTfbsClusteredWithCellsV3.bed",
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "tf", "score", "cells"]
)

# keep only  rows where the cell list contains GM12878

#this reduces the rows from ~4 million to ~1 million

df_gm = df[df["cells"].str.contains("GM12878", na=False)]

print(df_gm.head())
print("Total rows for GM12878:", len(df_gm))

df_gm.to_csv("GM12878_only.bed", sep="\t", header=False, index=False)

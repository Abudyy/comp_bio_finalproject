import pandas as pd

#thr  factorbook motif positions
df = pd.read_csv(
    r"C:\Users\saleh\561_project\data\factorbookMotifPos.txt",
    sep="\t",
    header=None,
    names=["ignore", "chrom", "start", "end", "tf", "score", "strand"]
)

# we only want to keep CTCF
df_ctcf = df[df["tf"] == "CTCF"]

print(df_ctcf.head())
print("Total CTCF sites:", len(df_ctcf))


df_ctcf[["chrom", "start", "end", "tf", "score", "strand"]].to_csv(
    r"C:\Users\saleh\561_project\data\CTCF_sites.bed",
    sep="\t",
    header=False,
    index=False
)

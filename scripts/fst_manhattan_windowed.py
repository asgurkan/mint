import pandas as pd
import matplotlib.pyplot as plt

# input 
fst_file = snakemake.input.windowed_fst

# output
fst_image = snakemake.output.fst_manhattan_image

df = pd.read_csv(fst_file, sep="\t")

# positioning the scatters to the middle of chromosome intervals in the figure 
df["BIN_MID"] = (df["BIN_START"] + df["BIN_END"]) / 2

# sorted list of chromosomes
try:
    chroms = sorted(df["CHROM"].unique(), key=lambda x: int(x))
except ValueError:
    chroms = sorted(df["CHROM"].unique())

# indexing each chromosome for coloring
chrom_index = {chrom: i for i, chrom in enumerate(chroms)}
df["chrom_idx"] = df["CHROM"].map(chrom_index)

# max BIN_END per chromosome to scale each to [0..1]
max_end_per_chrom = df.groupby("CHROM")["BIN_END"].max().to_dict()

def scale_position(row):
    i = chrom_index[row["CHROM"]]
    max_end = max_end_per_chrom[row["CHROM"]]
    scaled_within_chrom = row["BIN_MID"] / max_end  # scale to 0-1
    return i + scaled_within_chrom  # offset by chromosome index

df["scaled_pos"] = df.apply(scale_position, axis=1)

plt.figure(figsize=(20, 8))

# Scatter plot
scatter = plt.scatter(
    df["scaled_pos"],
    df["WEIGHTED_FST"],
    c=df["chrom_idx"],
    cmap="tab20",
    s=2,
    alpha=0.8
)

plt.title("Windowed Weir and Cockerham FST Plot", fontsize=14)
plt.xlabel("Locus", fontsize=12)
plt.ylabel("FST", fontsize=12)

# one block per chromosome
xticks = []
xlabels = []
for i, chrom in enumerate(chroms):
    xticks.append(i + 0.5)  # midpoint of [i, i+1]
    xlabels.append(chrom)

plt.xticks(xticks, xlabels, fontsize=8)
plt.xlim(0, len(chroms))
plt.tight_layout()
plt.savefig(fst_image, dpi = 300)
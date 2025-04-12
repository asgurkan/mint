import pandas as pd
import matplotlib.pyplot as plt

# input 
fst_file_path = snakemake.input.windowed_fst

# output
figure_output_path = snakemake.output.fst_manhattan_image
fst_file_path = "/home/asgurkan/Documents/population_genomics/workspace/march_study/fst_calculation/FST_dav_mysv2.weir.fst"
figure_output_path = "/home/asgurkan/Documents/population_genomics/report/figures/fst_manhattan/fst_manhattan_locus.png"

df = pd.read_csv(fst_file_path, sep = "\t")
df.columns = ["CHROM", "POS", "FST_raw"]

df["CHROM"] = df["CHROM"].astype(str)

# sort by chrom
df = df.sort_values(["CHROM", "POS"]).copy()

# indexing and sorting chromosomes to give them unique color 
unique_chroms = df["CHROM"].unique()
try:
    chroms_sorted = sorted(unique_chroms, key=lambda x: int(x))
except ValueError:
    chroms_sorted = sorted(unique_chroms)

chrom_index = {chrom: i for i, chrom in enumerate(chroms_sorted)}
df["chrom_idx"] = df["CHROM"].map(chrom_index)

# scale chromosomes by max position
max_pos_per_chrom = df.groupby("CHROM")["POS"].max().to_dict()

def scale_position(row):
    chrom = row["CHROM"]
    i = chrom_index[chrom]
    max_pos = max_pos_per_chrom[chrom]
    scaled_within_chrom = row["POS"] / max_pos  
    return i + scaled_within_chrom 

df["scaled_pos"] = df.apply(scale_position, axis=1)

# plot
plt.figure(figsize=(20, 8)) 

# scatter plot
plt.scatter(
    df["scaled_pos"],   # x: scaled_positions
    df["FST_raw"],      # y: fst value
    c=df["chrom_idx"],  # c: chr index 
    cmap="tab20",       # colormap
    s=2,                # scatter size
    alpha=0.5           # transparency
)

plt.title("Weir and Cockerham FST Plot by Locus", fontsize=14)
plt.xlabel("Locus", fontsize=12)
plt.ylabel("FST", fontsize=12)

# Each chromosome occupies [i, i+1], label the midpoint (i+0.5)
xticks = []
xlabels = []
for i, chrom in enumerate(chroms_sorted):
    xticks.append(i + 0.5)
    xlabels.append(chrom)

plt.xticks(xticks, xlabels, fontsize=8)
plt.xlim(0, len(chroms_sorted))  # from 0..n_chroms

plt.tight_layout()
plt.savefig(figure_output_path, dpi=300)

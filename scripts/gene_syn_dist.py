import pandas as pd
import matplotlib.pyplot as plt
import os

# input
csv_path = snakemake.input.annotations_csv

# output
output_dir = snakemake.output.locus_wise_syn_dir

os.makedirs(output_dir, exist_ok=True)
df = pd.read_csv(csv_path)

# Clean gene names and IsSynonymous column
df = df[df["Gene"].notna()]
df = df[df["Gene"].str.strip() != "."]
df = df[df["IsSynonymous"].notna()]
df["IsSynonymous"] = pd.to_numeric(df["IsSynonymous"], errors="coerce")
df = df[df["IsSynonymous"].isin([0, 1, 9])]

# Universal color map
color_map = {
    "Non-Synonymous (0)": "#ad2407",  # red
    "Synonymous (1)": "#4575b4",      # blue
    "Unknown (9)": "#aaaaaa"         # gray
}

# --- Per-locus plots ---
for locus, group in df.groupby("Locus"):
    gene_syn_counts = (
        group.groupby("Gene")["IsSynonymous"]
        .value_counts()
        .unstack(fill_value=0)
    )

    gene_syn_counts.columns = [
        "Non-Synonymous (0)" if c == 0 else
        "Synonymous (1)" if c == 1 else
        "Unknown (9)" for c in gene_syn_counts.columns
    ]

    # Ensure all columns exist in the right order for consistent coloring
    for col in color_map:
        if col not in gene_syn_counts.columns:
            gene_syn_counts[col] = 0
    gene_syn_counts = gene_syn_counts[list(color_map.keys())]

    top_genes = gene_syn_counts.sum(axis=1).nlargest(20).index
    gene_syn_counts = gene_syn_counts.loc[top_genes]

    ax = gene_syn_counts.plot(
        kind="bar", stacked=True, figsize=(12, 6), color=color_map
    )
    plt.title(f"Synonymous vs Non-Synonymous Variants in {locus} (Top 20 Genes)")
    plt.xlabel("Gene")
    plt.ylabel("Variant Count")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.legend(loc="upper right")

    for bars in ax.containers:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_y() + height / 2,
                    f"{int(height)}",
                    ha="center", va="center", fontsize=8, color="black"
                )
    locus_clean = locus.replace(":", "_").replace("/", "_")
    plt.savefig(f"{output_dir}/{locus_clean}_synonymous_variant_plot.png")
    plt.close()

# --- Global plot ---
global_gene_syn_counts = (
    df.groupby("Gene")["IsSynonymous"]
    .value_counts()
    .unstack(fill_value=0)
)

global_gene_syn_counts.columns = [
    "Non-Synonymous (0)" if c == 0 else
    "Synonymous (1)" if c == 1 else
    "Unknown (9)" for c in global_gene_syn_counts.columns
]

for col in color_map:
    if col not in global_gene_syn_counts.columns:
        global_gene_syn_counts[col] = 0
global_gene_syn_counts = global_gene_syn_counts[list(color_map.keys())]

top_global_genes = global_gene_syn_counts.sum(axis=1).nlargest(20).index
global_gene_syn_counts = global_gene_syn_counts.loc[top_global_genes]

ax = global_gene_syn_counts.plot(
    kind="bar", stacked=True, figsize=(12, 6), color=color_map
)
plt.title("Global Synonymous vs Non-Synonymous Variants (Top 20 Genes)")
plt.xlabel("Gene")
plt.ylabel("Variant Count")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.legend(loc="upper right")

for bars in ax.containers:
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_y() + height / 2,
                f"{int(height)}",
                ha="center", va="center", fontsize=8, color="black"
            )

plt.savefig(f"{output_dir}/global_synonymous_variant_plot.png")
plt.close()

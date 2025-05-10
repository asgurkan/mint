import re
import matplotlib.pyplot as plt
import pandas as pd

# Input: list of CV error files from Snakemake
cv_files = snakemake.input.cv_errors
output_csv = snakemake.output.adx_error_table
output_fig = snakemake.output.adx_error_figure

k_values = []
cv_errors = []

for path in cv_files:
    with open(path, "r") as f:
        for line in f:
            match = re.search(r"CV error \(K=(\d+)\):\s*([\d.]+)", line)
            if match:
                k = int(match.group(1))
                error = float(match.group(2))
                k_values.append(k)
                cv_errors.append(error)
                break  # Only the first match needed

# Sort results by K
df = pd.DataFrame(sorted(zip(k_values, cv_errors)), columns=["K", "CV_Error"])

# Save to CSV
df.to_csv(output_csv, index=False)

# Plot
plt.figure(figsize=(8, 5))
plt.plot(df["K"], df["CV_Error"], marker='o', linestyle='-')
for _, row in df.iterrows():
    plt.text(row["K"], row["CV_Error"] + 0.002, f"{row['CV_Error']:.4f}", 
             ha='center', va='bottom', fontsize=9)

plt.title("ADMIXTURE CV Error vs. K")
plt.xlabel("K")
plt.ylabel("CV Error")
plt.xticks(df["K"])
plt.grid(True)
plt.tight_layout()
plt.savefig(output_fig, dpi=300)

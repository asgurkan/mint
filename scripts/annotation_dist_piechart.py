import pandas as pd 
import matplotlib.pyplot as plt 

# input 
annotations_csv = snakemake.input.annotations_csv

# output
piechart = snakemake.output.piechart

### piecharts 

data = pd.read_csv(annotations_csv)
label_mapping = {0: "Nonsynonymous", 1: "Synonymous", 9: "Unknown"}
genic_label_mapping = {0: "Intergenic", 1: "Genic", 9: "Unknown"}

syn_counts = data["IsSynonymous"].value_counts().rename(index=label_mapping)
genic_counts = data["IsGenic"].value_counts().rename(index=genic_label_mapping)

colors_syn = ["#689be8", "#fc8d62", "#89e8aa"]  
colors_genic = ["#8da0cb", "#e78ac3", "#a6d854"]  

fig, axs = plt.subplots(1, 2, figsize=(12, 6))

axs[0].pie(syn_counts, labels=syn_counts.index, autopct="%1.1f%%", colors=colors_syn, textprops={'fontsize': 8})
axs[0].set_title("Synonymous vs. Nonsynonymous Mutations", fontsize = 10)

axs[1].pie(genic_counts, labels=genic_counts.index, autopct="%1.1f%%", colors=colors_genic, textprops={'fontsize': 8})
axs[1].set_title("Genic vs. Intergenic Mutations", fontsize = 10)
plt.savefig(piechart)
plt.close()

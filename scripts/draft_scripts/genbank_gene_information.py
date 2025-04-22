from Bio import SeqIO
import pandas as pd

def create_gene_data_csv(genbank_file, output_csv):
    """
    Extract gene data from a multi-locus GenBank file and save to CSV.

    :param genbank_file: Path to GenBank (.gbff or .gb) file with multiple loci.
    :param output_csv: Path to output CSV.
    """
    gene_data = []

    # Parse all records (each locus/contig)
    for record in SeqIO.parse(genbank_file, "genbank"):
        locus = record.id
        for feature in record.features:
            if feature.type == "gene":
                start = int(feature.location.start)
                end = int(feature.location.end)
                gene_name = feature.qualifiers.get("gene", ["Unknown"])[0]

                gene_data.append({
                    "locus": locus,
                    "gene": gene_name,
                    "start": start,
                    "end": end
                })

    # Convert to DataFrame and save
    gene_df = pd.DataFrame(gene_data)
    gene_df.to_csv(output_csv, index=False)
    print(f"✅ Gene data saved to: {output_csv}")

# Example usage
genbank_file = "/home/asgurkan/Documents/population_genomics/GCF_014108235.1_mMyoMyo1.p_genomic.gbff"
output_csv = "/home/asgurkan/Documents/population_genomics/data/mysdav/gene_data_all_loci.csv"

create_gene_data_csv(genbank_file, output_csv)

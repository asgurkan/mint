import csv
from cyvcf2 import VCF

# input 
vcf_path = snakemake.input.annotated_vcf
csv_path = snakemake.output.annotations_csv

vcf = VCF(vcf_path)

with open(csv_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Locus", "Position", "Gene", "IsSynonymous", "IsTransition", "IsGenic", "IsPseudo"])

    for i, variant in enumerate(vcf):
        writer.writerow([
            variant.CHROM,
            variant.POS,
            variant.INFO.get("Gene", "."),
            variant.INFO.get("IsSynonymous", "."),
            variant.INFO.get("IsTransition", "."),
            variant.INFO.get("IsGenic", "."),
            variant.INFO.get("IsPseudo", ".")
        ])

        if i % 10000 == 0:
            print(f"Processed {i} variants")


from Bio import SeqIO

# input
gb_file = "/home/asgurkan/Documents/population_genomics/GCF_014108235.1_mMyoMyo1.p_genomic.gbff"
vcf_in = "data/mysdav/variants/filtered.vcf"
# output
vcf_out = "data/mysdav/variants/filtered_renamed.vcf"

loci = [record.id for record in SeqIO.parse(gb_file, "genbank")]

vcf_to_gb_map = {str(i + 1): loci[i] for i in range(len(loci))}

with open(vcf_in) as fin, open(vcf_out, "w") as fout:
    for line in fin:
        if line.startswith("##contig=<ID="):
            # Replace contig header lines
            for k, v in vcf_to_gb_map.items():
                if f"ID={k}" in line:
                    line = line.replace(f"ID={k}", f"ID={v}")
        elif line.startswith("#"):
            fout.write(line)
            continue
        else:
            cols = line.strip().split("\t")
            if cols[0] in vcf_to_gb_map:
                cols[0] = vcf_to_gb_map[cols[0]]
            line = "\t".join(cols) + "\n"
        fout.write(line)

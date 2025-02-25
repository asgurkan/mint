import pandas as pd

## Functions
samples = pd.read_csv("samples.tsv", sep="\t")

def get_samples():
    """Returns list of samples."""
    return list(samples['sample'].unique())

def get_projects():
    """Returns list of projects."""
    return list(samples['proje'].unique())

## Rules    
rule all:
    input:
        "data/pca/pca.eigenval"
    
rule vcf_filter:
    input: 
        raw_vcf = "data/raw/{sample}.vcf"
    output: 
        filtered_vcf = "data/filtered/{sample}/{sample}.vcf"
    conda: 
        "envs/bcftools.yaml"
    params: 
        filter_expr = 'N_ALT==1 && TYPE=="snp" && INFO/DP>5 && CHROM!="NC_029346.1" && QUAL>=30'
    shell:
        """
        bcftools filter -i '{params.filter_expr}' {input.raw_vcf} -o {output.filtered_vcf}
        """

rule bgzip_and_index:
    input:
        "data/filtered/{sample}/{sample}.vcf"
    output:
        vcf_gz="data/filtered/{sample}/{sample}.vcf.gz",
        index="data/filtered/{sample}/{sample}.vcf.gz.tbi"
    conda:
        "envs/bcftools.yaml"
    shell:
        "bgzip -c {input} > {output.vcf_gz} && tabix -p vcf {output.vcf_gz}"

rule merge_vcfs:
    input:
        vcfs = expand("data/filtered/{sample}/{sample}.vcf.gz", sample=get_samples())
    output:
        merged = "data/merged/merged.vcf"
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools merge {input.vcfs} -o {output.merged} -O v"

rule pop_variant_filter:
    input:
        merged_vcf = "data/merged/merged.vcf"
    output:
        filtered_vcf = "data/merged/population_filtered.vcf"
    conda:
        "envs/bcftools.yaml"
    params:
        filter_expr = 'MAF>=0.05 && COUNT(GT=="./.")/N_SAMPLES<0.5'
    shell:
        "bcftools filter -i '{params.filter_expr}' {input.merged_vcf} -o {output.filtered_vcf}"

rule ld_pruning:
    input:
        vcf = "data/merged/population_filtered.vcf"
    output:
        prune_in = "data/ld_pruning/temp.prune.in",
        prune_out = "data/ld_pruning/temp.prune.out",
        bed = "data/ld_pruning/temp.bed",
        bim = "data/ld_pruning/temp.bim",
        fam = "data/ld_pruning/temp.fam"
    conda:
        "envs/plink.yaml"
    shell:
        """
        mkdir -p data/ld_pruning
        plink --vcf {input.vcf} --make-bed --allow-extra-chr --out data/ld_pruning/temp
        plink --bfile data/ld_pruning/temp --allow-extra-chr --indep-pairwise 50kb 1 0.2 --out data/ld_pruning/temp
        """

rule extract_pruned_variants:
    input:
        bed = "data/ld_pruning/temp.bed",
        bim = "data/ld_pruning/temp.bim",
        fam = "data/ld_pruning/temp.fam",
        prune_in = "data/ld_pruning/temp.prune.in"
    output:
        bed = "data/ld_pruning/pruned.bed",
        bim = "data/ld_pruning/pruned.bim",
        fam = "data/ld_pruning/pruned.fam"
    conda:
        "envs/plink.yaml"
    shell:
        """
        plink --bfile data/ld_pruning/temp \
              --extract {input.prune_in} \
              --make-bed \
              --allow-extra-chr \
              --out data/ld_pruning/pruned
        """

rule pca:
    input:
        bed = "data/ld_pruning/temp.bed",
        bim = "data/ld_pruning/temp.bim",
        fam = "data/ld_pruning/temp.fam"
    output:
        eigenvec = "data/pca/pca.eigenvec",
        eigenval = "data/pca/pca.eigenval"
    conda:
        "envs/plink.yaml"
    shell:
        """
        mkdir -p data/pca
        plink --bfile data/ld_pruning/temp \
              --allow-extra-chr \
              --pca \
              --out data/pca/pca
        """

rule admixture:
    input:
        bed = "data/ld_pruning/pruned.bed",
        bim = "data/ld_pruning/pruned.bim",
        fam = "data/ld_pruning/pruned.fam"
    output:
        cv_error = "data/admixture/admixture_cv_error.txt"
    params:
        K = 3  
    conda:
        "envs/admixture.yaml"
    shell:
        """
        mkdir -p data/admixture
        admixture --cv {input.bed} {params.K} | tee {output.cv_error}
        """

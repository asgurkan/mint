import pandas as pd
import glob

configfile: "config.yml"

### Variables 
samples = pd.read_csv("sample.tsv", sep = "\t")

### Functions
def get_samples():
    """Returns list of samples."""
    return list(samples["sample"].unique())

def get_projects():
    """Returns list of projects."""
    return list(samples["project"].unique())

def get_raw_vcf():
    """ Returns raw vcf path. """
    return samples["path"].unique()

# population file paths variable creation and control
pop_files = sorted(glob.glob("population_files/*.txt"))
if len(pop_files) < 2:
    raise ValueError("At least two population files are required for FST analysis.")

# Rule all
rule all:
    input:
        expand("data/{project}_report/figures/pca", project = get_projects()) ,
        # expand("data/raw/{sample}/{sample}.vcf", sample=get_samples()),
        expand(["data/admixture/cv_error/admixture_cv_error_K{K}.txt"], K = config["admixture"]["k_value"]),
        expand("data/{project}_report/figures/fst/fst_manhattan_locus.png", project = get_projects()),
        "data/dxy_results.csv"

## Rules 
rule vcf_split_sample:
    input:
        raw_vcf = get_raw_vcf()[0]
    output:
        splitted_vcf = "data/raw/{sample}/{sample}.vcf"
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        mkdir -p data/raw/{wildcards.sample}
        bcftools view -c1 -s {wildcards.sample} -Oz -o {output.splitted_vcf} {input.raw_vcf} || \
        echo "Warning: No variants found for {wildcards.sample}. Creating empty file." && touch {output.splitted_vcf}
        """

# rule vcf_filter:
#     input:
#         raw_vcf = "data/raw/{sample}/{sample}.vcf"
#     output:
#         filtered_vcf = "data/filtered/{sample}/{sample}.vcf",
#     conda:
#         "envs/bcftools.yaml"
#     params:
#         filter_expr='N_ALT==1 && QUAL>=30 && CHROM!~"NC_029346.1" && CHROM!~"Un" && CHROM!~"MT"'

#     shell:
#         """
#         bcftools filter -i '{params.filter_expr}' {input.raw_vcf} -o {output.filtered_vcf}
#         """


rule bgzip_and_index:
    input:
        "data/raw/{sample}/{sample}.vcf"
    output:
        vcf_gz="data/filtered/{sample}/{sample}.vcf.gz",
        index="data/filtered/{sample}/{sample}.vcf.gz.tbi",
    conda:
        "envs/bcftools.yaml"
    shell:
        "bgzip -c {input} > {output.vcf_gz} && tabix -p vcf {output.vcf_gz}"


rule merge_vcfs:
    input:
        vcfs=expand("data/filtered/{sample}/{sample}.vcf.gz", sample=get_samples()),
    output:
        merged="data/merged/merged.vcf",
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools merge {input.vcfs} -o {output.merged} -O v"


rule pop_variant_filter:
    input:
        merged_vcf="data/merged/merged.vcf",
    output:
        filtered_vcf="data/merged/population_filtered.vcf",
    conda:
        "envs/bcftools.yaml"
    params:
        filter_expr='MAF>=0.05 && COUNT(GT=="./.")/N_SAMPLES<0.5',
    shell:
        "bcftools filter -i '{params.filter_expr}' {input.merged_vcf} -o {output.filtered_vcf}"


rule ld_pruning:
    input:
        # vcf="data/merged/population_filtered.vcf",
        vcf="data/merged/merged.vcf",
    output:
        prune_in="data/ld_pruning/temp.prune.in",
        prune_out="data/ld_pruning/temp.prune.out",
        bed="data/ld_pruning/temp.bed",
        bim="data/ld_pruning/temp.bim",
        fam="data/ld_pruning/temp.fam",
    conda:
        "envs/plink.yaml"
    shell:
        """
        mkdir -p data/ld_pruning
        plink --vcf {input.vcf} \
              --make-bed \
              --allow-extra-chr \
              --chr-set 92 \
              --double-id \
              --out data/ld_pruning/temp

        plink --bfile data/ld_pruning/temp \
              --allow-extra-chr \
              --chr-set 92 \
              --double-id \
              --indep-pairwise 50kb 1 0.2 \
              --out data/ld_pruning/temp
        """

rule extract_pruned_variants:
    input:
        bed="data/ld_pruning/temp.bed",
        bim="data/ld_pruning/temp.bim",
        fam="data/ld_pruning/temp.fam",
        prune_in="data/ld_pruning/temp.prune.in",
    output:
        bed="data/ld_pruning/pruned.bed",
        bim="data/ld_pruning/pruned.bim",
        fam="data/ld_pruning/pruned.fam",
    conda:
        "envs/plink.yaml"
    shell:
        """
        plink --bfile data/ld_pruning/temp \
              --extract {input.prune_in} \
              --make-bed \
              --chr-set 92 \
              --allow-extra-chr \
              --out data/ld_pruning/pruned
        """

rule pca:
    input:
        bed="data/ld_pruning/temp.bed",
        bim="data/ld_pruning/temp.bim",
        fam="data/ld_pruning/temp.fam",
    output:
        eigenvec="data/pca/pca.eigenvec",
        eigenval="data/pca/pca.eigenval",
    conda:
        "envs/plink.yaml"
    shell:
        """
        mkdir -p data/pca
        plink --bfile data/ld_pruning/temp \
              --allow-extra-chr \
              --chr-set 92 \
              --pca \
              --out data/pca/pca
        """

rule pca_graph:
    input: 
        eigenvec_path = "data/pca/pca.eigenvec",
        eigenval_path = "data/pca/pca.eigenval",
        sample_tsv_path = "sample.tsv"
    output:
        pca_output_dir = directory("data/{project}_report/figures/pca")
    conda:
        "envs/pca_graph.yaml"
    log:
        "logs/{project}_pca_plot.log"
    script:
        "scripts/pca_plot.R"

rule admixture:
    input:
        bed="data/ld_pruning/pruned.bed",
        bim="data/ld_pruning/pruned.bim",
        fam="data/ld_pruning/pruned.fam",
    output:
        cv_error="data/admixture/cv_error/admixture_cv_error_K{K}.txt",
        output=directory("data/admixture/K{K}"),
        P="data/admixture/K{K}/pruned.{K}.P",
        Q="data/admixture/K{K}/pruned.{K}.Q"
    params:
        K=config["admixture"]["k_value"]
    conda:
        "envs/admixture.yaml"
    shell:
        """
        admixture --cv {input.bed} {params.K} | tee {output.cv_error}
        mv pruned.{params.K}.P {output.P}
        mv pruned.{params.K}.Q {output.Q}
        """
        
rule fst:
    input:
        input_vcf = "/home/asgurkan/Documents/population_genomics/march_data/raw/mysdav_renamed.vcf"
    output:
        windowed_fst = "data/fst/{project}.windowed.weir.fst"
    params:
        window_size = config["fst"]["window-size"],
        window_step = config["fst"]["window-step"],
        weir_fst_pops = " ".join([f"--weir-fst-pop {pop_file}" for pop_file in pop_files]),
        prefix = "data/fst/{project}"
    conda:
        "envs/vcftools.yaml"
    shell:
     """
        vcftools \
            --vcf {input.input_vcf} {params.weir_fst_pops} \
            --fst-window-size {params.window_size} \
            --fst-window-step {params.window_step} \
            --out {params.prefix}
        """

rule fst_manhattan:
    input:
        windowed_fst = "data/fst/{project}.windowed.weir.fst"
    output:
        fst_manhattan_image = "data/{project}_report/figures/fst/fst_manhattan_locus.png"
    conda:
        "envs/py_visualization.yaml"
    script:
        "scripts/fst_manhattan_windowed.py"

rule pi_dxy_calculation:
    input:
        merged_vcf = "/home/asgurkan/Documents/population_genomics/data/merged/merged.vcf.gz",
        population1 = pop_files[0],
        population2 = pop_files[1]
    output:
        dxy_text_file = "data/dxy_results.csv",
        pi_text_file = "data/pi_results.csv"
    conda:
        "envs/pca_graph.yaml"
    log:
        dxy_log_file = "logs/dxy_log_file.txt",
        pi_log_file = "logs/pi_log_file.txt"
    script:
        "scripts/pi_dxy_calculation.R"
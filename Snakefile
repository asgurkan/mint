import pandas as pd
import glob

configfile: "config.yml"

### Variables 
samples = pd.read_csv("sample.tsv", sep = "\t")

### Functions
def get_samples():
    """Returns list of samples."""
    return list(samples["sample"].unique())

def get_project():
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
        expand("data/{project}/ld_pruning/temp.prune.in", project = get_project())
        #expand("data/{project}/variants/filtered.vcf.gz", project = get_project())


# Rules 
rule retrieve_raw_vcf:
    input: 
        raw_vcf = get_raw_vcf()
    output: 
        raw_vcf_project = "data/{project}/variants/raw.vcf"
    shell:
        "cp {input.raw_vcf} {output.raw_vcf_project}"

rule pop_variant_filter:
    input:
        raw_vcf = "data/{project}/variants/raw.vcf"
    output:
        filtered_vcf = "data/{project}/variants/filtered.vcf"
    conda:
        "envs/bcftools.yaml"
    params:
        filter_expr='MAF>=0.05 && COUNT(GT=="./.")/N_SAMPLES<0.5',
    shell:
        "bcftools filter -i '{params.filter_expr}' {input.raw_vcf} -o {output.filtered_vcf}"


rule gzip_merge:
    input:
        merged = "data/{project}/variants/filtered.vcf"
    output:
        merged_gz = "data/{project}/variants/filtered.vcf.gz"
    shell:
        "bgzip -k {input.merged} -o {output.merged_gz}"


rule ld_pruning:
    input:
        # vcf = "data/{project}/variants/filtered.vcf"      # check filtering exp. 
        vcf = "data/{project}/variants/raw.vcf"
    output:
        prune_in="data/{project}/ld_pruning/temp.prune.in",
        prune_out="data/{project}/ld_pruning/temp.prune.out",
        bed="data/{project}/ld_pruning/temp.bed",
        bim="data/{project}/ld_pruning/temp.bim",
        fam="data/{project}/ld_pruning/temp.fam"
    params: 
        temp_folder = "data/{project}/ld_pruning/temp"
    conda:
        "envs/plink.yaml"
    shell:
        """
        plink --vcf {input.vcf} \
              --make-bed \
              --allow-extra-chr \
              --chr-set 92 \
              --double-id \
              --out {params.temp_folder}

        plink --bfile {params.temp_folder}\
              --allow-extra-chr \
              --chr-set 92 \
              --double-id \
              --indep-pairwise 50kb 1 0.2 \
              --out {params.temp_folder}
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
        merged_vcf = "data/merged/merged.vcf.gz",
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
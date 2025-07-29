import pandas as pd
import glob

configfile: "config.yml"

# paths 
gb_file_path = "gb_file_path.txt"
sample_tsv_path = "sample.tsv"
K_VALUES = list(map(str, range(1, 6))) 
### Variables 
samples = pd.read_csv(sample_tsv_path, sep = "\t")

### Functions
def get_samples():
    """Returns list of samples."""
    return list(samples["sample"].unique())

def get_project():
    """Returns list of projects."""
    return list(samples["project"].unique())

def get_vcf_path():
    """Returns list of projects."""
    return list(samples["path"].unique())

def get_raw_vcf():
    """ Returns raw vcf path. """
    return samples["path"].unique()

def ref_gb():
    f = open("gb_file_path.txt", "r")
    gb_file_path = f.readline().rstrip("\n")
    return gb_file_path

# population file paths variable creation and control
pop_files = sorted(glob.glob("population_files/*.txt"))
if len(pop_files) < 2:
    raise ValueError("At least two population files are required for FST analysis.")

# Rule all
rule all:
    input:
        expand( "data/{project}/report_complete{K}.txt", sample=get_samples(), project = get_project(), K=K_VALUES),
        #expand("data/{project}/figures/fst/fst_manhattan_locus.png", project = get_project())
        #expand("data/{project}/figures/admixture/plots/admix_plot_K_{K}.png", project = get_project(), K=K_VALUES)
        #expand("data/{project}/report_complete{K}.txt", project = get_project(), K=K_VALUES) 

# Rules 

# rule fastp:
#     input:
#         r1 = "/home/asgurkan/Documents/population_genomics/workspace/files_from_bengisu/RP23-238-22020_S131_L003_R1_001.fastq.gz",
#         r2 = "/home/asgurkan/Documents/population_genomics/workspace/files_from_bengisu/RP23-238-22020_S131_L003_R2_001.fastq.gz"
#     output:
#         r1_trimmed = "data/raw/{sample}/clean/{sample}_R1.trimmed.fastq.gz",
#         r2_trimmed = "data/raw/{sample}/clean/{sample}_R2.trimmed.fastq.gz",
#         html = "data/raw/{sample}/qc/{sample}_fastp.html",
#         json = "data/raw/{sample}/qc/{sample}_fastp.json"
#     threads: 12
#     conda:
#         "envs/fastp.yaml"
#     shell:
#         """
#         fastp \
#             -i {input.r1} -I {input.r2} \
#             -o {output.r1_trimmed} -O {output.r2_trimmed} \
#             -h {output.html} -j {output.json} \
#             -w {threads} --detect_adapter_for_pe
#         """

# rule bwa_mem_pe:
#     input:
#         ref = "workspace/GCF_000327345.1_ASM32734v1_genomic.fna",
#         r1 = "data/raw/{sample}/clean/{sample}_R1.trimmed.fastq.gz",
#         r2 = "data/raw/{sample}/clean/{sample}_R2.trimmed.fastq.gz"
#     output:
#         sam = "data/{sample}/aligned/{sample}.pe.sam"
#     threads: 20
#     conda:
#         "envs/bwa.yaml"
#     shell:
#         """
#         bwa mem -t {threads} -L 1000,1000 {input.ref} {input.r1} {input.r2} > {output.sam}
#         """


# rule retrieve_raw_vcf:
#     input: 
#         raw_vcf = get_raw_vcf()
#     output: 
#         raw_vcf_project = "data/{project}/variants/raw.vcf"
#     shell:
#         "cp {input.raw_vcf} {output.raw_vcf_project}"

# rule pop_variant_filter:
#     input:
#         raw_vcf = "data/{project}/variants/raw.vcf"
#     output:
#         filtered_vcf = "data/{project}/variants/filtered.vcf"
#     conda:
#         "envs/bcftools.yaml"
#     params:
#         filter_expr='MAF>=0.05'
#         # filter_expr='MAF>=0.05 && COUNT(GT=="./.")/N_SAMPLES<0.5'
#     shell:
#         "bcftools filter -i '{params.filter_expr}' {input.raw_vcf} -o {output.filtered_vcf}"

# rule rename_vcf_locus:
#     input:
#         vcf = "data/{project}/variants/filtered.vcf",
#         ref_gb = ref_gb()
#     output:
#         annotated_vcf = "data/{project}/variants/filtered_renamed.vcf"
#     conda:
#         "envs/vcf_annotator.yaml"
#     script:
#         "scripts/rename_vcf.py"

# rule annotate_vcf:
#     input:
#         vcf = "data/{project}/variants/filtered_renamed.vcf",
#         ref_gb = ref_gb()
#     output:
#         annotated_vcf = "data/{project}/annotations/annotated.vcf"
#     conda:
#         "envs/vcf_annotator.yaml"
#     shell:
#         """
#         vcf-annotator --output {output.annotated_vcf} {input.vcf} {input.ref_gb} 
#         """

# rule vcf2csv:
#     input:
#         annotated_vcf = "data/{project}/annotations/annotated.vcf"
#     output:
#         annotations_csv = "data/{project}/annotations/annotations.csv"
#     conda:
#         "envs/vcf_annotator.yaml"
#     script:
#         "scripts/annotation2csv.py"

# rule annotation_dist_piechart:
#     input:
#         annotations_csv = "data/{project}/annotations/annotations.csv"
#     output:
#         piechart = "data/{project}/figures/annotation_dist_piechart.png"
#     conda:
#         "envs/py_visualization.yaml"
#     script:
#         "scripts/annotation_dist_piechart.py"

# rule gene_syn_dist:
#     input:
#         annotations_csv = "data/{project}/annotations/annotations.csv"
#     output:
#         locus_wise_syn_dir = directory("data/{project}/figures/locus_wise_syn_analysis")
#     conda:
#         "envs/py_visualization.yaml"
#     script:
#         "scripts/gene_syn_dist.py"

# rule gzip_merge:
#     input:
#         merged = "data/{project}/variants/filtered.vcf"
#     output:
#         merged_gz = "data/{project}/variants/filtered.vcf.gz"
#     shell:
#         "bgzip -k {input.merged} -o {output.merged_gz}"



############################################################
# rule heterozygosity:
#     input:
#         vcf = "results/filtered/{sample}.vcf"
#     output:
#         tsv = "results/heterozygosity/{sample}_heterozygosity.tsv",
#         log = "results/heterozygosity/{sample}_heterozygosity.log"
#     params:
#         script = "scripts/heterozygosity.sh"
#     shell:
#         """
#         bash {params.script} {input.vcf} {output.tsv} > /dev/null 2>&1
#         """
############################################################

rule ld_pruning:
    input:
        vcf = get_vcf_path()
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
        bed="data/{project}/ld_pruning/temp.bed",
        bim="data/{project}/ld_pruning/temp.bim",
        fam="data/{project}/ld_pruning/temp.fam",
        prune_in="data/{project}/ld_pruning/temp.prune.in",
    output:
        bed="data/{project}/ld_pruning/pruned.bed",
        bim="data/{project}/ld_pruning/pruned.bim",
        fam="data/{project}/ld_pruning/pruned.fam",
    conda:
        "envs/plink.yaml"
    params: 
        temp_folder = "data/{project}/ld_pruning/temp",
        pruned_folder = "data/{project}/ld_pruning/pruned"
    shell:
        """
        plink --bfile {params.temp_folder} \
              --extract {input.prune_in} \
              --make-bed \
              --chr-set 92 \
              --allow-extra-chr \
              --out {params.pruned_folder}
        """

rule pca:
    input:
        bed = "data/{project}/ld_pruning/pruned.bed",
        bim = "data/{project}/ld_pruning/pruned.bim",
        fam = "data/{project}/ld_pruning/pruned.fam",
    output:
        eigenvec = "data/{project}/pca/pca.eigenvec",
        eigenval = "data/{project}/pca/pca.eigenval",
    conda:
        "envs/plink.yaml"
    params: 
        pca_folder = "data/{project}/pca",
        temp_folder = "data/{project}/ld_pruning/pruned"
    shell:
        """
        mkdir -p {params.pca_folder}
        plink --bfile {params.temp_folder} \
              --allow-extra-chr \
              --chr-set 92 \
              --pca \
              --out {params.pca_folder}/pca
        """

rule pca_graph:
    input: 
        eigenvec_path = "data/{project}/pca/pca.eigenvec",
        eigenval_path = "data/{project}/pca/pca.eigenval",
        sample_tsv = sample_tsv_path
    output:
        pca_output_dir = directory("data/{project}/figures/pca")
    conda:
        "envs/pca_graph.yaml"
    log:
        "logs/{project}_pca_plot.log"
    script:
        "scripts/pca_plot.R"

rule admixture:
    input:
        bed = "data/{project}/ld_pruning/pruned.bed",
        bim = "data/{project}/ld_pruning/pruned.bim",
        fam = "data/{project}/ld_pruning/pruned.fam",
    output:
        cv_error = "data/{project}/admixture/K{K}/admixture_cv_error_K{K}.txt",
        output = directory("data/{project}/admixture/K{K}"),
        P = "data/{project}/admixture/K{K}/pruned.{K}.P",
        Q = "data/{project}/admixture/K{K}/pruned.{K}.Q"
    params:
        K = lambda wildcards: wildcards.K
    conda:
        "envs/admixture.yaml"
    shell:
        """
        admixture --cv {input.bed} {params.K} | tee {output.cv_error}
        mv pruned.{params.K}.P {output.P}
        mv pruned.{params.K}.Q {output.Q}
        """

rule admixture_visualization:
    input:
        cv_errors = expand("data/{project}/admixture/K{K}/admixture_cv_error_K{K}.txt",
                            K=K_VALUES, 
                            project=get_project())
    output:
        adx_error_table = "data/{project}/admixture/admixture_cv_errors.csv",
        adx_error_figure = "data/{project}/figures/admixture/error_plot.png"
    conda:
        "envs/py_visualization.yaml"
    script:
        "scripts/admixture_error_line_vis.py"


rule plot_admixture:
    input:
        fam = "data/{project}/ld_pruning/pruned.fam",
        qs = "data/{project}/admixture/K{K}/pruned.{K}.Q"
    output:
        plot_admixture_png = "data/{project}/figures/admixture/plots/admix_plot_K_{K}.png"
    params:
        q_prefix = "data/{project}/admixture",  
        out_dir = "data/{project}/figures/admixture/plots",
        max_k = len(K_VALUES)
    conda:
        "envs/admixture_plot_env.yaml"
    shell:
        """
        Rscript scripts/plot_admixture.R {input.fam} {params.q_prefix} {params.max_k} {params.out_dir}
        """

rule fst:
    input:
        input_vcf = get_vcf_path()
    output:
        windowed_fst = "data/{project}/fst/{project}.windowed.weir.fst"
    params:
        window_size = config["fst"]["window-size"],
        window_step = config["fst"]["window-step"],
        weir_fst_pops = " ".join([f"--weir-fst-pop {pop_file}" for pop_file in pop_files]),
        prefix = "data/{project}/fst/{project}"
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
        windowed_fst = "data/{project}/fst/{project}.windowed.weir.fst"
    output:
        fst_manhattan_image = "data/{project}/figures/fst/fst_manhattan_locus.png"
    conda:
        "envs/py_visualization.yaml"
    script:
        "scripts/fst_manhattan_windowed.py"

rule pi_dxy_calculation:
    input:
        filtered_vcf = "/home/asgurkan/Documents/population_genomics/march_data/raw/mysdav_renamed.vcf.gz",
        population1 = pop_files[0],
        population2 = pop_files[1]
    output:
        dxy_text_file = "data/{project}/dxy_results.csv",
        pi_text_file = "data/{project}/pi_results.csv"
    conda:
        "envs/pca_graph.yaml"
    log:
        dxy_log_file = "logs/{project}/dxy_log_file.txt",
        pi_log_file = "logs/{project}/pi_log_file.txt"
    script:
        "scripts/pi_dxy_calculation.R"

rule report:
    input:
        pca = directory("data/{project}/figures/pca"),
        adx = directory("data/{project}/admixture/K{K}"),     
        adx_error_table = "data/{project}/admixture/admixture_cv_errors.csv",
        fst_fig = "data/{project}/figures/fst/fst_manhattan_locus.png"
    output:
        dummy_output = "data/{project}/report_complete{K}.txt"  
    shell:
        """
        echo "Generating report for project {wildcards.project} with K={wildcards.K}" > {output}
        echo "PCA dir: {input.pca}" >> {output}
        echo "Admixture dir: {input.adx}" >> {output}
        echo "Admixture error file: {input.adx_error_table}" >> {output}
        echo "FST figure: {input.fst_fig}" >> {output}
        """

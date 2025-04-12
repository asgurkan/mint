
# MINT Pipeline

This pipeline is designed for population genomics analysis, including tasks such as SNP analysis, PCA generation, and other bioinformatics workflows.

## Requirements

Before running the pipeline, make sure you have **Conda** installed on your system. Conda is used to manage the environment and dependencies for this pipeline.

### Setup Instructions

1. **Create the Conda Environment**:
   
   The pipeline depends on specific packages and software. To set up the Conda environment with all required dependencies, use the provided `setup.yaml` file. This file contains all necessary information to create the environment.

   To create the Conda environment, run the following command:

   ```bash
   conda env create -f setup.yaml
   ```

   This will create a Conda environment called `mint` and install all the required dependencies, including `bcftools`, `vcftools`, and other tools necessary for the pipeline.

2. **Activate the Conda Environment**:
   
   Once the environment is created, activate it using:

   ```bash
   conda activate mint
   ```

   This will ensure that all the required software and dependencies are available when running the pipeline.

### Running the Sample File Generation Script

Before running the Snakemake pipeline, you need to **generate the sample file** and **population files**. These files are required as input for the pipeline.

1. **Run the Sample File Generation Script**:

   The `sample_file_generation.sh` script generates the required `sample.tsv` and population files, which are needed for the pipeline.

   Usage:

   ```bash
   ./sample_file_generation.sh <input_vcf> <project_name>
   ```

   Example:

   ```bash
   ./sample_file_generation.sh example.vcf MyProject
   ```

   This script will:
   - Validate the VCF file.
   - Generate the `sample.tsv` file containing sample, population, project, and path columns.
   - Create individual population files based on the population information in the `sample.tsv`.

### Running the Snakemake Pipeline

After generating the necessary files, you can now run the Snakemake pipeline to execute the genomic analysis.

To start the pipeline:

```bash
snakemake --use-conda --conda-frontend conda --cores 1
```

This will:
- Use the Conda environment specified in the `setup.yaml`.
- Execute the pipeline with the defined workflow and rules.

---

## Summary of Commands

1. **Create Conda environment**:
   ```bash
   conda env create -f setup.yaml
   ```

2. **Activate the environment**:
   ```bash
   conda activate mint
   ```

3. **Run the sample file generation script**:
   ```bash
   ./sample_file_generation.sh <input_vcf> <project_name>
   ```

4. **Run the Snakemake pipeline**:
   ```bash
   snakemake --use-conda --conda-frontend conda --cores 1
   ```

---

### Troubleshooting

- If you encounter issues during environment setup or package installation, make sure that your Conda installation is up to date. You can update Conda with:

  ```bash
  conda update conda
  ```

- Ensure that all required tools (`bcftools`, `vcf-validator`) are available in the Conda environment. If not, try reactivating the environment:

  ```bash
  conda activate mint
  ```
#!/usr/bin/env bash
#
# Sample TSV creation script.
# This script validates a VCF file, extracts sample names,
# and generates a "sample.tsv" file with specific columns.
# Additionally, it creates population files based on the population column.
#
# Usage:
#   ./create_sample_tsv.sh <input_vcf> <project_name>
#
# Example:
#   ./create_sample_tsv.sh example.vcf MyProject
#

# -------- Script Settings --------
set -o errexit   # Exit immediately if a command exits with a non-zero status
set -o nounset   # Treat unset variables as an error when substituting
set -o pipefail  # Pipeline returns the exit status of the last command that fails

# Default LOG_LEVEL
LOG_LEVEL="${LOG_LEVEL:-INFO}"

# -------- Logging Function -------
log() {
    local level="$1"
    shift
    local message="$*"
    local timestamp
    timestamp=$(date +"%Y-%m-%d %H:%M:%S")

    case "$level" in
        INFO)
            [[ "$LOG_LEVEL" =~ ^(INFO|WARNING)$ ]] && echo "[$timestamp] [INFO]    $message"
            ;;
        WARNING)
            [[ "$LOG_LEVEL" =~ ^(INFO|WARNING)$ ]] && echo "[$timestamp] [WARNING] $message"
            ;;
        ERROR)
            echo "[$timestamp] [ERROR]   $message" >&2
            ;;
    esac
}

# -------- Usage Function ---------
usage() {
    echo "Usage: $0 <input_vcf> <project_name>"
    exit 1
}

# ------ Argument Validation ------
if [[ $# -lt 2 ]]; then
    usage
fi

INPUT_VCF="$1"
PROJECT_NAME="$2"

# ------ Command Availability -----
if ! command -v bcftools &>/dev/null; then
    log "ERROR" "bcftools not found. Please install bcftools or add it to your PATH."
    exit 1
fi

if ! command -v vcf-validator &>/dev/null; then
    log "ERROR" "vcf-validator not found. Please install vcf-validator or add it to your PATH."
    exit 1
fi

# ------ File Existence Check -----
if [[ ! -f "$INPUT_VCF" ]]; then
    log "ERROR" "VCF file not found: $INPUT_VCF"
    exit 1
fi

# ------ vcf-validator Step -------
log "INFO" "Validating the VCF file with vcf-validator..."
VALIDATION_OUTPUT=$(vcf-validator "$INPUT_VCF" 2>&1)
VALIDATION_EXIT_CODE=$?

if [[ $VALIDATION_EXIT_CODE -eq 0 ]]; then
    log "INFO" "VCF file is valid. No errors or warnings reported."
elif [[ $VALIDATION_EXIT_CODE -eq 1 ]]; then
    log "WARNING" "vcf-validator reported warnings:"
    echo "$VALIDATION_OUTPUT"
    log "INFO" "Continuing despite warnings."
else
    log "ERROR" "vcf-validator reported errors:"
    echo "$VALIDATION_OUTPUT" >&2
    exit 1
fi

# ------ bcftools Readability -----
# Suppress bcftools' own error output to ensure only one error message if invalid
if ! bcftools query -l "$INPUT_VCF" >/dev/null 2>&1; then
    log "ERROR" "Could not read from $INPUT_VCF with bcftools. Please ensure it is a valid VCF file."
    exit 1
fi

# --------- Script Start ----------
log "INFO" "Script started."
log "INFO" "VCF file: $INPUT_VCF"
log "INFO" "Project name: $PROJECT_NAME"

# ---- Creating the sample.tsv ----
log "INFO" "Creating sample.tsv..."
echo -e "sample\tpopulation\tproject\tpath" > sample.tsv

# This query is guaranteed to succeed now, as we checked above
bcftools query -l "$INPUT_VCF" | awk -v project="$PROJECT_NAME" -v path="$INPUT_VCF" -F"_" '{
    print $0 "\t" $(NF-1)"_"$NF "\t" project "\t" path
}' >> sample.tsv

log "INFO" "sample.tsv created successfully!"

# ---- Create population files ----
log "INFO" "Creating population files..."
# Read population names from the sample.tsv and create individual population files
awk -F'\t' 'NR > 1 {print $2}' sample.tsv | sort | uniq | while read population; do
    # Create a population file for each population
    population_file="population_files/$population.txt"
    mkdir -p population_files  # Ensure the directory exists
    awk -v pop="$population" -F'\t' '$2 == pop {print $1}' sample.tsv > "$population_file"
    log "INFO" "Population file created: $population_file"
done

log "INFO" "Population files created successfully!"

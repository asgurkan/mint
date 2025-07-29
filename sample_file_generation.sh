#!/usr/bin/env bash
#
# Sample TSV and population files creation script (without vcf-validator).
#
# Usage:
#   ./create_sample_tsv.sh <input_vcf> <project_name>

# -------- Script Settings --------
set -o errexit
set -o nounset
set -o pipefail

LOG_LEVEL="${LOG_LEVEL:-INFO}"

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

usage() {
    echo "Usage: $0 <input_vcf> <project_name>"
    exit 1
}

if [[ $# -lt 2 ]]; then
    usage
fi

INPUT_VCF="$1"
PROJECT_NAME="$2"

if ! command -v bcftools &>/dev/null; then
    log "ERROR" "bcftools not found. Please install bcftools or add it to your PATH."
    exit 1
fi

if [[ ! -f "$INPUT_VCF" ]]; then
    log "ERROR" "VCF file not found: $INPUT_VCF"
    exit 1
fi

# ----- VCF Structure Check (bcftools header read) -----
if ! bcftools view -h "$INPUT_VCF" >/dev/null 2>&1; then
    log "ERROR" "Invalid VCF file: $INPUT_VCF. Header could not be read."
    exit 1
fi

# ----- Sample Existence Check -----
if ! bcftools query -l "$INPUT_VCF" >/dev/null 2>&1; then
    log "ERROR" "Could not extract samples from $INPUT_VCF. Ensure it contains sample columns."
    exit 1
fi

log "INFO" "Script started."
log "INFO" "VCF file: $INPUT_VCF"
log "INFO" "Project name: $PROJECT_NAME"

log "INFO" "Creating sample.tsv..."
echo -e "sample\tpopulation\tproject\tpath" > sample.tsv

bcftools query -l "$INPUT_VCF" | awk -v project="$PROJECT_NAME" -v path="$INPUT_VCF" -F"_" '{
    print $0 "\t" $(NF-1)"_"$NF "\t" project "\t" path
}' >> sample.tsv

log "INFO" "sample.tsv created successfully!"

log "INFO" "Creating population files..."
awk -F'\t' 'NR > 1 {print $2}' sample.tsv | sort | uniq | while read population; do
    population_file="population_files/$population.txt"
    mkdir -p population_files
    awk -v pop="$population" -F'\t' '$2 == pop {print $1}' sample.tsv > "$population_file"
    log "INFO" "Population file created: $population_file"
done

log "INFO" "Population files created successfully!"
log "INFO" "Script completed successfully!"

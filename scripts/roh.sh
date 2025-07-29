#!/usr/bin/env bash

# -----------------------------------
# ROH Analysis Script (bcftools)
# Usage: ./roh_analysis.sh input.vcf(.gz) output_prefix
# Output:
#   - {prefix}_roh.txt
#   - {prefix}_af.vcf.gz
#   - {prefix}_roh.log
# -----------------------------------

VCF_FILE="$1"
PREFIX="$2"
LOG_FILE="${PREFIX}_roh.log"
AF_FILE="${PREFIX}_af.vcf.gz"
ROH_FILE="${PREFIX}_roh.txt"

# -------- Logging Function --------
log() {
    local level="$1"
    shift
    local message="$*"
    local timestamp
    timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[$timestamp] [$level] $message" >> "$LOG_FILE"
}

# -------- Check Inputs --------
if [[ -z "$VCF_FILE" || -z "$PREFIX" ]]; then
    echo "Usage: ./roh_analysis.sh input.vcf(.gz) output_prefix"
    exit 1
fi

if [[ ! -f "$VCF_FILE" ]]; then
    log "ERROR" "VCF file not found: $VCF_FILE"
    exit 1
fi

log "INFO" "Starting ROH analysis..."
log "INFO" "Input VCF: $VCF_FILE"
log "INFO" "Output prefix: $PREFIX"

# -------- Compress and Index VCF --------
if [[ "$VCF_FILE" != *.vcf.gz ]]; then
    log "INFO" "Compressing VCF with bgzip..."
    bgzip -c "$VCF_FILE" > "${VCF_FILE%.vcf}.vcf.gz"
    VCF_FILE="${VCF_FILE%.vcf}.vcf.gz"
    log "INFO" "Compressed VCF: $VCF_FILE"
fi

if [[ ! -f "$VCF_FILE.tbi" ]]; then
    log "INFO" "Indexing VCF..."
    bcftools index "$VCF_FILE"
fi

# -------- Generate AF File (VCF Format!) --------
# Generate AF file with empty ID field
log "INFO" "Creating allele frequency file (VCF format, ID=.)..."
bcftools +fill-tags "$VCF_FILE" -Ou -- -t AF | \
bcftools annotate -x ID -Ou | \
bcftools view -Oz -o "$AF_FILE"
bcftools index -f "$AF_FILE"


# -------- Run ROH Analysis --------
log "INFO" "Running bcftools roh..."
bcftools roh --AF-file "$AF_FILE" -G30 "$VCF_FILE" > "$ROH_FILE"

log "INFO" "ROH output written to: $ROH_FILE"
log "INFO" "✅ ROH analysis completed successfully."

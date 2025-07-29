#!/usr/bin/env bash

# -------------------------------
# Heterozygosity Calculator Script
# -------------------------------
# Calculates heterozygosity metrics per sample from a multi-sample VCF file
# Logs and results are written to output files

# --------- Setup ---------------
VCF_FILE="$1"
OUTPUT_FILE="${2:-heterozygosity_summary.tsv}"
LOG_FILE="${OUTPUT_FILE%.tsv}.log"

# --------- Logging Function ----------
log() {
    local level="$1"
    shift
    local message="$*"
    local timestamp
    timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo -e "[$timestamp] [$level] $message" >> "$LOG_FILE"
}

# --------- Input Validation ---------
if [[ -z "$VCF_FILE" ]]; then
    log "ERROR" "No VCF file provided. Usage: ./heterozygosity.sh input.vcf [output.tsv]"
    exit 1
fi

if [[ ! -f "$VCF_FILE" ]]; then
    log "ERROR" "VCF file not found: $VCF_FILE"
    exit 1
fi

log "INFO" "Processing VCF file: $VCF_FILE"
log "INFO" "Output TSV file will be: $OUTPUT_FILE"
log "INFO" "Log file: $LOG_FILE"

# --------- Get Sample Names ----------
log "INFO" "Extracting sample names..."
sample_names=($(bcftools query -l "$VCF_FILE"))
sample_count=${#sample_names[@]}

if [[ "$sample_count" -eq 0 ]]; then
    log "ERROR" "No samples found in the VCF file."
    exit 1
fi

log "INFO" "Found $sample_count samples."

# --------- Export Sample Names for AWK ----------
for i in "${!sample_names[@]}"; do
    export SN$((i+1))="${sample_names[$i]}"
done

# --------- Write Header ---------------
echo -e "Sample\tHeterozygous_Count\tHomozygous_Count\tHeterozygosity_Ratio" > "$OUTPUT_FILE"
log "INFO" "Header written to $OUTPUT_FILE"

# --------- Process Genotype Data -------
log "INFO" "Starting genotype parsing and heterozygosity calculation..."
bcftools query -f '[%GT\t]\n' "$VCF_FILE" | \
awk -v count="$sample_count" -v out="$OUTPUT_FILE" '
BEGIN {
    for (i = 1; i <= count; i++) {
        het[i] = 0;
        hom[i] = 0;
    }
}
{
    for (i = 1; i <= NF; i++) {
        if ($i ~ /0[\/|]1|1[\/|]0/) het[i]++;
        else if ($i ~ /0[\/|]0|1[\/|]1/) hom[i]++;
    }
}
END {
    for (i = 1; i <= count; i++) {
        name = ENVIRON["SN"i];
        total = het[i] + hom[i];
        ratio = (total > 0) ? het[i] / total : 0;
        printf "%s\t%d\t%d\t%.4f\n", name, het[i], hom[i], ratio >> out;
    }
}'

log "INFO" "Heterozygosity summary successfully written to $OUTPUT_FILE"

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(devtools)
install_github("https://github.com/MiGurke/GenoPop")

library(GenoPop)
library(vcfR)

# input
vcf_gz_path <- snakemake@input[["filtered_vcf"]]
population1_path <- snakemake@input[["population1"]]
population2_path <- snakemake@input[["population2"]]

# output
dxy_text_file <- snakemake@output[["dxy_text_file"]]
pi_text_file <- snakemake@output[["pi_text_file"]]

vcf_data <- read.vcfR(vcf_gz_path)

contig_info <- vcf_data@meta[grepl("^##contig", vcf_data@meta)]

if(length(contig_info) > 0) {
  # Extract chromosome names and lengths
  chrom_lengths <- data.frame(
    Chromosome = sapply(contig_info, function(x) sub(".*ID=([^,]+).*", "\\1", x)),
    Length = as.numeric(sapply(contig_info, function(x) sub(".*length=([0-9]+).*", "\\1", x)))
  )
  total_sequence_length <- sum(chrom_lengths$Length)
} else {
  print("No contig information found in VCF metadata.")
}

# ### Dxy ###
# pop1_individuals <- readLines(population1_path)
# pop2_individuals <- readLines(population2_path)

# Read population individuals with trimmed whitespace
pop1_individuals <- trimws(readLines(population1_path))
pop2_individuals <- trimws(readLines(population2_path))

# Check against VCF samples
vcf_samples <- colnames(vcf_data@gt)[-1]
missing_pop1 <- setdiff(pop1_individuals, vcf_samples)
missing_pop2 <- setdiff(pop2_individuals, vcf_samples)
if (length(missing_pop1) > 0 || length(missing_pop2) > 0) {
  stop("Samples missing in VCF: Pop1 - ", paste(missing_pop1, collapse=", "), 
       "; Pop2 - ", paste(missing_pop2, collapse=", "))
}

vcf_samples <- colnames(vcf_data@gt)[-1]  # Skip FORMAT column
print("VCF Samples:")
print(vcf_samples)
print("Pop1 Samples:")
print(pop1_individuals)
print("Pop2 Samples:")
print(pop2_individuals)

dxy_windows <- Dxy(vcf_gz_path,
                   pop1_individuals,
                   pop2_individuals,
                   seq_length = total_sequence_length,
                   window_size = 1000000, skip_size = 0)
write.csv(dxy_windows, dxy_text_file, row.names = FALSE)

### Pi ###
pi_windows <- Pi(vcf_gz_path,
                 seq_length = total_sequence_length,
                 window_size = 1000000, skip_size = 0)
                 
write.csv(pi_windows, pi_text_file, row.names = FALSE)
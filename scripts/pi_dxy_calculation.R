options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(devtools)
install_github("https://github.com/MiGurke/GenoPop")

library(GenoPop)
library(vcfR)

# input
vcf_gz_path <- snakemake@input[["merged_vcf"]]
population1_path <- snakemake@input[["population1"]]
population2_path <- snakemake@input[["population2"]]

# output
dxy_text_file <- snakemake@output[["dxy_text_file"]]
pi_text_file <- snakemake@output[["pi_text_file"]]

# logs
dxy_log_file  <- snakemake@log[["dxy_log_file"]]
pi_log_file  <- snakemake@log[["pi_log_file"]]

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

### Dxy ###
pop1_individuals <- readLines(population1_path)
pop2_individuals <- readLines(population2_path)

dxy_windows <- Dxy(vcf_gz_path,
                   pop1_individuals,
                   pop2_individuals,
                   seq_length = total_sequence_length,
                   window_size = 1000000, skip_size = 0,
                   logfile = dxy_log_file)

write.csv(dxy_windows, dxy_text_file, row.names = FALSE)

### Pi ###

pi_windows <- Pi(vcf_gz_path,
                 seq_length = total_sequence_length,
                 window_size = 1000000, skip_size = 0,
                 logfile = pi_log_file)
                 
write.csv(pi_windows, pi_text_file, row.names = FALSE)
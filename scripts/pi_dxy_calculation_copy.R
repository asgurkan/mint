install.packages("vcfR")
library(GenoPop)
library(vcfR)

vcf_gz_path <- "/Users/ahmetsametgurkan/Documents/graduation_files/data/raw_data/mysdav_renamed.vcf.gz"
#vcf_gz_path <- "/Users/ahmetsametgurkan/Documents/graduation_files/filtered.vcf.gz"

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

pop1_individuals <- readLines("/Users/ahmetsametgurkan/Documents/graduation_files/data/population_files/davidii.txt")  # Path to your pop1 individual file
pop2_individuals <- readLines("/Users/ahmetsametgurkan/Documents/graduation_files/data/population_files/mystacinus.txt") 

dxy_windows <- Dxy(vcf_gz_path, pop1_individuals, pop2_individuals, seq_length = total_sequence_length,
                   window_size = 1000000, skip_size = 0, logfile="dxy_logfile.txt")

write.csv(dxy_windows, "dxy_windows.csv", row.names = FALSE)


### Pi ### 
pi_windows <- Pi(vcf_gz_path, seq_length = total_sequence_length,
                 window_size = 1000000, skip_size = 0, logfile="pi_logfile.txt")

write.csv(pi_windows, "pi_windows.csv", row.names = FALSE)

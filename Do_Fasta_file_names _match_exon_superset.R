# Set your working directory to where your files are located
# setwd("/path/to/your/directory")
setwd("../Data/GeneSet/fasta/")

gtf_file <- "./Homo_sapiens.GRCh38.111.gtf"

fasta_files <- list.files(pattern = "\\.fa$", full.names = TRUE)
length(fasta_files)

chromosomes_from_fasta <- sapply(strsplit(basename(fasta_files), "\\."), function(x) x[3])

# Extract chromosome names from FASTA files
chromosomes_from_fasta <- sapply(strsplit(basename(fasta_files), "[._]"), function(x) {
  # Look for the element containing "chromosome" and get the element following it
  idx <- which(x == "chromosome")
  if(length(idx) > 0 && idx < length(x)) {
    return(x[idx + 1])
  } else {
    return(NA) # Return NA if no matching pattern is found
  }
})

# Clean up NA values if any
chromosomes_from_fasta <- chromosomes_from_fasta[!is.na(chromosomes_from_fasta)]

# Proceed with the rest of the script as before

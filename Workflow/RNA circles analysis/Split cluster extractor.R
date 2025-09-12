library(data.table)
library(stringr)
library(Biostrings) 

fasta_file <- "/Users/awsms1/RNA lab/EC2_big_boi/orfs.cluster_size_over_10.fasta"
alignment_dir <- "/Users/awsms1/RNA lab/EC2_big_boi/ID35_clustering/rawsize_splits/uclust_10-19.uc"
output_dir <- "/Users/awsms1/RNA lab/extracted clusters/10-19"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

fasta <- readAAStringSet(fasta_file)
names(fasta) <- str_trim(names(fasta))  


alignment_files <- alignment_dir

for (afile in alignment_files) {
 
  uc <- fread(afile, sep="\t", header=FALSE, quote="", data.table=TRUE)
  
 
  clusters <- unique(uc[V1 == "C", V2])
  
 
  for (cl in clusters) {
    
    source_ids <- unique(uc[V2 == cl, V9])
    
    
    fasta_ids <- names(fasta)
    matches <- fasta[fasta_ids %in% source_ids]
    
    if (length(matches) > 0) {

      out_file <- file.path(output_dir, paste0("cluster_", cl, ".fasta"))
      
      writeXStringSet(matches, out_file)
      
      cat("Written cluster", cl, "to", out_file, "with", length(matches), "sequences\n")
    }
    
    }
  }


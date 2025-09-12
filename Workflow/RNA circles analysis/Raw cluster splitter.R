library(data.table)

input_file <- "/Users/awsms1/RNA lab/EC2 big boi/ID35 clustering/circ_clustersID35.uc"
output_dir <- "/Users/awsms1/RNA lab/EC2 big boi/ID35 clustering/rawsize_splits"
dir.create(output_dir, showWarnings = FALSE)

uc <- fread(input_file, sep="\t", header=FALSE, quote="", data.table=TRUE)

clusters <- uc[V1 == "C", .(Cluster = V2, Size = as.numeric(V3))]

clusters[, Bin := fifelse(Size >= 500, "500+",
                          fifelse(Size >= 100,
                                  paste0((Size %/% 100) * 100, "-", ((Size %/% 100) * 100 + 99)),
                                  paste0((Size %/% 10) * 10, "-", ((Size %/% 10) * 10 + 9))))]

print_cluster_counts <- function(cluster_dt) {
  cat("Unique clusters per raw size bin:\n")
  counts <- cluster_dt[, .(NumClusters = uniqueN(Cluster)), by = Bin][order(-Bin)]
  print(counts)
}

print_cluster_counts(clusters)

for (bin_name in unique(clusters$Bin)) {
  clusters_in_bin <- clusters[Bin == bin_name, Cluster]
  
  subset_uc <- uc[V2 %in% clusters_in_bin]
  subset_uc <- subset_uc[order(V2)]
  
  output_file <- file.path(output_dir, paste0("uclust_", bin_name, ".uc"))
  fwrite(subset_uc, output_file, sep="\t", quote=FALSE, col.names=FALSE)
  
  cat("Written bin", bin_name, "to", output_file, "with", nrow(subset_uc), "rows\n")
}



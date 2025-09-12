library(data.table)
library(Biostrings)

uc_file <- "/Users/awsms1/RNA lab/EC2_big_boi/ID35_clustering/circ_clustersID35.uc"
fasta_file <- "/Users/awsms1/RNA lab/EC2_big_boi/ID35_clustering/circ_centroidsID35.fasta"
fasta_out <- "/Users/awsms1/RNA lab/EC2_big_boi/ID35_clustering/cenDB.fasta"

fasta <- readAAStringSet(fasta_file)
header_dt <- data.table(
  SeqID = sub("^>", "", names(fasta)),  
  Sequence = as.character(fasta)
)

header_dt[, CoreID := sub(" .*", "", SeqID)]

uc <- fread(uc_file, sep="\t", header=FALSE, data.table=TRUE)
colnames(uc) <- c("Type", "Cluster", "Col3", "Col4", "Col5",
                  "Col6", "Col7", "Col8", "SeqID", "Col10")

centroids_uc <- uc[Type == "S", .(Cluster, SeqID)]
centroids_uc[, CoreID := sub(" .*", "", SeqID)]

cluster_sizes <- uc[, .N, by=Cluster]
setnames(cluster_sizes, "N", "ClusterSize")

centroids_uc <- merge(centroids_uc, cluster_sizes, by="Cluster", all.x=TRUE)
final_dt <- merge(header_dt, centroids_uc, by="CoreID", all.x=TRUE)

final_dt[, FastaHeader := paste0(">", Cluster, "|", CoreID, "|size=", ClusterSize)]

fasta_lines <- c(rbind(final_dt$FastaHeader, final_dt$Sequence))
writeLines(fasta_lines, fasta_out)

cat("FASTA saved to", fasta_out, "\n")

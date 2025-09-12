library(data.table)
library(Biostrings)

uc_file <- "/Users/awsms1/RNA lab/EC2 big boi/ID35 clustering/circ_clustersID35.uc"
fasta_file <- "/Users/awsms1/RNA lab/EC2 big boi/ID35 clustering/circ_centroidsID35.fasta"
output_file <- "/Users/awsms1/RNA lab/EC2 big boi/ID35 clustering/linked_centroids.tsv"

fasta <- readAAStringSet(fasta_file)
header_dt <- data.table(
  SeqID = sub("^>", "", names(fasta)),
  Sequence = as.character(fasta)
)

# Extract CoreID
header_dt[, CoreID := sub(" .*", "", SeqID)]
header_dt[, c("Accession", "CircleInfo") := tstrsplit(CoreID, "_circle_", fixed=TRUE)]

uc <- fread(uc_file, sep="\t", header=FALSE, data.table=TRUE)
colnames(uc) <- c("Type", "Cluster", "Col3", "Col4", "Col5",
                  "Col6", "Col7", "Col8", "SeqID", "Col10")

centroids_uc <- uc[Type == "S", .(Cluster, SeqID)]
centroids_uc[, CoreID := sub(" .*", "", SeqID)]

# Compute cluster sizes
cluster_sizes <- uc[, .N, by=Cluster]
setnames(cluster_sizes, "N", "Size")

# Merge cluster info
centroids_uc <- merge(centroids_uc, cluster_sizes, by="Cluster", all.x=TRUE)

final_dt <- merge(header_dt, centroids_uc, by="CoreID", all.x=TRUE)

final_dt[, SizeBin := fifelse(Size >= 500, "500+",
                              fifelse(Size >= 400, "400-499",
                                      fifelse(Size >= 300, "300-399",
                                              fifelse(Size >= 200, "200-299",
                                                      fifelse(Size >= 100, "100-199",
                                                              "<100")))))]

# Rename CoreID to linked_centroids
setnames(final_dt, "CoreID", "linked_centroids")

fwrite(final_dt[, .(linked_centroids, Cluster, Size, Accession, CircleInfo, Sequence, SizeBin)],
       output_file, sep="\t")
cat("Combined TSV saved to", output_file, "\n")

size_bins <- unique(final_dt$SizeBin)

for (bin in size_bins) {
  bin_dt <- final_dt[SizeBin == bin]
  bin_file <- sub("\\.tsv$", paste0("_", bin, ".tsv"), output_file)
  fwrite(bin_dt[, .(linked_centroids, Cluster, Size, Accession, CircleInfo, Sequence, SizeBin)],
         bin_file, sep="\t")
  cat("Saved", nrow(bin_dt), "centroids to", bin_file, "\n")
}


library(data.table)
library(stringr)

# ----- Step 1: Read DIAMOND output -----
diamond_file <- "/Users/awsms1/RNA lab/Centroid obelisk hits/Oblin1_results.m8"
diamond_dt <- fread(diamond_file, header = FALSE, sep = "\t")

# ----- Step 2: Name columns -----
colnames(diamond_dt) <- c(
  "QueryID", "ClusterID", "PercIdentity", "AlignLength", "Mismatches", 
  "GapOpens", "QStart", "QEnd", "SStart", "SEnd", "Evalue", "BitScore", "Extra"
)

# ----- Step 3: Extract cluster size from ClusterID -----
diamond_dt[, ClusterSize := as.integer(str_extract(ClusterID, "(?<=size=)\\d+"))]

# ----- Step 4: Extract SRA accession from ClusterID -----
diamond_dt[, SRA := str_extract(ClusterID, "[A-Z]{3}\\d+_\\d+")]

# ----- Step 5: Convert numeric cluster ID for sorting -----
diamond_dt[, cluster := as.numeric(str_extract(ClusterID, "^\\d+"))]

# ----- Step 6: Remove unnecessary columns -----
diamond_dt[, c("ClusterID", "Extra") := NULL]

# ----- Step 7: Sort by cluster and BitScore -----
setorder(diamond_dt, cluster, -BitScore)

# ----- Step 8: Save sorted output & top hits -----
fwrite(
  diamond_dt,
  "/Users/awsms1/RNA lab/Centroid obelisk hits/Oblin_Hits_by_cluster.tsv",
  sep = "\t"
)

top_hits <- diamond_dt[, .SD[1], by = cluster]

fwrite(
  top_hits,
  "/Users/awsms1/RNA lab/Centroid obelisk hits/Oblin_top_hits_per_cluster.tsv",
  sep = "\t"
)

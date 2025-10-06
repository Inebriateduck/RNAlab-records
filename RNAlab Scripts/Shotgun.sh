# Input files
UNCOMPRESSED=$1
COMPRESSED=$2
# Temp file to hold SRRs
TMP_SRR=$(mktemp)
# Extract SRRs from column 9 of uncompressed TSV (everything before first "_")
awk '{split($9,a,"_"); print a[1]}' "$UNCOMPRESSED" > "$TMP_SRR"
# Stream compressed TSV, match SRRs, and print yeast/C_elegans flags
zstd -d -c "$COMPRESSED" | awk -v SRR_FILE="$TMP_SRR" '
BEGIN {
    # Load SRRs into array
    while ((getline line < SRR_FILE) > 0) {
        srr[line] = 1
    }
    print "SRR\tYeast\tC_elegans"  # header
}
{
    srr_id = $1
    if (srr[srr_id]) {
        yeast_flag = ($4 == 1) ? 1 : 0
        celeg_flag = ($5 == 1) ? 1 : 0
        print srr_id "\t" yeast_flag "\t" celeg_flag
        delete srr[srr_id]  # remove SRR from future searches
    }
}'
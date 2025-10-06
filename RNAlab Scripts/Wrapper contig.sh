#!/bin/bash

# Wrapper to read 'Accession' from a CSV/TSV and extract nucleotide sequence
# Usage: ./getContig_wrapper_clean.sh -i input_file.csv -o output_file.csv

INPUT_FILE=""
OUTPUT_CSV="all_contigs.csv"

while getopts "i:o:" opt; do
  case $opt in
    i) INPUT_FILE="$OPTARG" ;;
    o) OUTPUT_CSV="$OPTARG" ;;
    *) echo "Usage: $0 -i input_file -o output_file"; exit 1 ;;
  esac
done

if [ -z "$INPUT_FILE" ] || [ ! -f "$INPUT_FILE" ]; then
    echo "Error: input file not specified or does not exist."
    exit 1
fi

# Clean line endings (Windows -> Unix)
INPUT_FILE_CLEAN=$(mktemp)
tr -d '\r' < "$INPUT_FILE" | sed '/^$/d' > "$INPUT_FILE_CLEAN"

# Detect delimiter (comma or tab)
HEADER=$(head -n 1 "$INPUT_FILE_CLEAN")
if [[ "$HEADER" == *","* ]]; then
    DELIM=","
else
    DELIM=$'\t'
fi

# Find Accession column index
COL_ACC=$(echo "$HEADER" | tr "$DELIM" '\n' | grep -nx "Accession" | cut -d: -f1)
if [ -z "$COL_ACC" ]; then
    echo "Error: 'Accession' column not found in header."
    rm "$INPUT_FILE_CLEAN"
    exit 1
fi

echo "Contig_ID,Sequence" > "$OUTPUT_CSV"

# Process each row (skip header)
tail -n +2 "$INPUT_FILE_CLEAN" | while IFS="$DELIM" read -r -a FIELDS; do
    CONTIG_ID="${FIELDS[$((COL_ACC-1))]}"
    [ -z "$CONTIG_ID" ] && continue

    echo "Processing $CONTIG_ID..."

    TMP=$(mktemp)

    # Run getContig.bash for this accession
    ./getContig.bash -i "$CONTIG_ID" -o "$TMP" 2>/dev/null

    if [ -s "$TMP" ]; then
        # Extract only nucleotide sequence (remove header, join lines)
        SEQ_OUTPUT=$(grep -v '^>' "$TMP" | tr -d '\n' | sed 's/"/""/g')
        echo "\"$CONTIG_ID\",\"$SEQ_OUTPUT\"" >> "$OUTPUT_CSV"
    else
        echo "\"$CONTIG_ID\",\"\"" >> "$OUTPUT_CSV"
    fi

    rm "$TMP"
done

rm "$INPUT_FILE_CLEAN"

echo "All contigs processed. Output saved to $OUTPUT_CSV"

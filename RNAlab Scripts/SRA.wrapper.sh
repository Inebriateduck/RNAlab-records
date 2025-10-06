#!/bin/bash

# Wrapper to read 'Accession' from one or more CSV/TSV files
# Each input file produces a FASTA file named <basename>_contigs.fa
# Usage: ./getContig_wrapper_clean.sh file1.csv file2.tsv ...

if [ $# -eq 0 ]; then
    echo "Usage: $0 input_file1 [input_file2 ...]"
    exit 1
fi

for INPUT_FILE in "$@"; do
    if [ ! -f "$INPUT_FILE" ]; then
        echo "Skipping missing file: $INPUT_FILE"
        continue
    fi

    BASENAME=$(basename "$INPUT_FILE")
    BASENAME_NOEXT="${BASENAME%.*}"
    OUTPUT_FA="${BASENAME_NOEXT}_contigs.fa"

    echo "Processing $INPUT_FILE → $OUTPUT_FA"

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
        echo "Error: 'Accession' column not found in $INPUT_FILE"
        rm "$INPUT_FILE_CLEAN"
        continue
    fi

    # Empty output file
    > "$OUTPUT_FA"

    # Process each row (skip header)
    tail -n +2 "$INPUT_FILE_CLEAN" | while IFS="$DELIM" read -r -a FIELDS; do
        CONTIG_ID="${FIELDS[$((COL_ACC-1))]}"
        [ -z "$CONTIG_ID" ] && continue

        echo "  Processing $CONTIG_ID..."

        TMP=$(mktemp)

        # Run getContig.bash for this accession
        ./getContig.bash -i "$CONTIG_ID" -o "$TMP" 2>/dev/null

        if [ -s "$TMP" ]; then
            # Extract nucleotide sequence only
            SEQ_OUTPUT=$(grep -v '^>' "$TMP" | tr -d '\n')
            {
                echo ">$CONTIG_ID"
                echo "$SEQ_OUTPUT"
            } >> "$OUTPUT_FA"
        else
            echo ">$CONTIG_ID" >> "$OUTPUT_FA"
        fi

        rm "$TMP"
    done

    rm "$INPUT_FILE_CLEAN"
    echo "Finished $INPUT_FILE → $OUTPUT_FA"
done

echo "All files processed."

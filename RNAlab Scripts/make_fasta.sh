#!/bin/bash

INPUT="$1"
OUTPUT_SMALL="$2"  # clusters_size_le5.fasta
OUTPUT_LARGE="$3"  # clusters_size_gt5.fasta

# Validate inputs
if [[ -z "$INPUT" || -z "$OUTPUT_SMALL" || -z "$OUTPUT_LARGE" ]]; then
    echo "Usage: $0 <input_file> <small_output> <large_output>"
    exit 1
fi

if [[ ! -f "$INPUT" ]]; then
    echo "Error: Input file '$INPUT' not found"
    exit 1
fi

echo "Processing $INPUT..."
echo "This may take a while for large files..."

# Clear output files
> "$OUTPUT_SMALL"
> "$OUTPUT_LARGE"

awk -v small="$OUTPUT_SMALL" -v large="$OUTPUT_LARGE" '
BEGIN { 
    FS="\t"; OFS="\t"
    processed = 0
}

NR>1 { 
    cluster=$2
    type=$1
    seq=$17
    
    # Progress indicator every 10000 lines
    if(++processed % 10000 == 0) {
        printf "Processed %d lines...\n", processed > "/dev/stderr"
    }
    
    # Count cluster size (include S and H, exclude C)
    if(type!="C") cluster_count[cluster]++
    
    # Store sequence only if type is S
    if(type=="S") {
        # Clean sequence: remove any character not A,C,G,T (case-insensitive)
        gsub(/[^ACGTacgt]/,"",seq)
        cluster_seq[cluster]=seq
    }
}

END {
    printf "Finished reading input. Writing output files...\n" > "/dev/stderr"
    
    # Process each cluster and write immediately to avoid memory issues
    cluster_num = 0
    for(cluster in cluster_seq) {
        cluster_num++
        if(cluster_num % 1000 == 0) {
            printf "Writing cluster %d...\n", cluster_num > "/dev/stderr"
        }
        
        size = cluster_count[cluster]
        header = ">cluster_num=" cluster " cluster_size=" size
        sequence = cluster_seq[cluster]
        
        # Choose output file
        output_file = (size > 5) ? large : small
        
        # Write header
        print header > output_file
        
        # Write sequence in 80-character lines (more efficient method)
        seq_len = length(sequence)
        for(i = 1; i <= seq_len; i += 80) {
            end_pos = (i + 79 <= seq_len) ? i + 79 : seq_len
            print substr(sequence, i, end_pos - i + 1) > output_file
        }
    }
    
    # Close files explicitly
    close(small)
    close(large)
    
    printf "Output writing complete.\n" > "/dev/stderr"
}
' "$INPUT"

echo "Splitting complete. Checking output files..."
echo "Small clusters (≤5): $(grep -c "^>" "$OUTPUT_SMALL" 2>/dev/null || echo 0) sequences"
echo "Large clusters (>5): $(grep -c "^>" "$OUTPUT_LARGE" 2>/dev/null || echo 0) sequences"
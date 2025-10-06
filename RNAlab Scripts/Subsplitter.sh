#!/bin/bash

# Default species if not set via environment variable
SPECIES="${FILTER_SPECIES:-saccharomyces}"

# Default cutoff (can be percentage or count)
CUTOFF="${FILTER_CUTOFF:-50}"

# Default mode: "percent" or "count"
MODE="${FILTER_MODE:-percent}"

INPUT="$1"
OUTPUT_SMALL="$2"   # clusters_size_lt3.fasta (< 3)
OUTPUT_MEDIUM="$3"  # clusters_size_3to5.fasta (3-5)
OUTPUT_LARGE="$4"   # clusters_size_gt5.fasta (> 5)

# Validate inputs
if [[ -z "$INPUT" || -z "$OUTPUT_SMALL" || -z "$OUTPUT_MEDIUM" || -z "$OUTPUT_LARGE" ]]; then
    echo "Usage: $0 <input_file> <small_output> <medium_output> <large_output>"
    echo ""
    echo "To change species filter, set FILTER_SPECIES environment variable:"
    echo "  FILTER_SPECIES=saccharomyces $0 ..."
    echo "  FILTER_SPECIES=celegans $0 ..."
    echo "  FILTER_SPECIES=drosophila $0 ..."
    echo ""
    echo "To filter by percentage (default), set FILTER_MODE=percent:"
    echo "  FILTER_CUTOFF=75 $0 ...  (include only clusters with >75% presence)"
    echo "  FILTER_CUTOFF=25 $0 ...  (include only clusters with >25% presence)"
    echo ""
    echo "To filter by number of hits, set FILTER_MODE=count:"
    echo "  FILTER_MODE=count FILTER_CUTOFF=5 $0 ...  (include only clusters with >5 hits)"
    echo "  FILTER_MODE=count FILTER_CUTOFF=10 $0 ... (include only clusters with >10 hits)"
    echo ""
    echo "Defaults: saccharomyces, cutoff=50, mode=percent"
    exit 1
fi

# Set column and label based on species
case "$SPECIES" in
    saccharomyces)
        COLUMN=11
        LABEL="Saccharomyces"
        LABEL_SHORT="saccharomyces"
        ;;
    celegans)
        COLUMN=12
        LABEL="C. elegans"
        LABEL_SHORT="c_elegans"
        ;;
    drosophila)
        COLUMN=13
        LABEL="Drosophila"
        LABEL_SHORT="drosophila"
        ;;
    *)
        echo "Error: Unknown species '$SPECIES'"
        echo "Valid options: saccharomyces, celegans, drosophila"
        exit 1
        ;;
esac

# Validate mode
if [[ "$MODE" != "percent" && "$MODE" != "count" ]]; then
    echo "Error: FILTER_MODE must be 'percent' or 'count' (got '$MODE')"
    exit 1
fi

# Validate cutoff is a number
if ! [[ "$CUTOFF" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    echo "Error: FILTER_CUTOFF must be a number (got '$CUTOFF')"
    exit 1
fi

if [[ ! -f "$INPUT" ]]; then
    echo "Error: Input file '$INPUT' not found"
    exit 1
fi

echo "Processing $INPUT..."
if [[ "$MODE" == "percent" ]]; then
    echo "Filtering for clusters with >$CUTOFF% $LABEL presence (column $COLUMN)..."
else
    echo "Filtering for clusters with >$CUTOFF $LABEL hits (column $COLUMN)..."
fi

# Clear output files
> "$OUTPUT_SMALL"
> "$OUTPUT_MEDIUM"
> "$OUTPUT_LARGE"

awk -v col="$COLUMN" -v label="$LABEL_SHORT" -v cutoff="$CUTOFF" -v mode="$MODE" -v small="$OUTPUT_SMALL" -v medium="$OUTPUT_MEDIUM" -v large="$OUTPUT_LARGE" '
BEGIN { 
    FS="\t"; OFS="\t"
    processed = 0
}

NR>1 { 
    cluster=$2
    type=$1
    seq=$17
    target_species=$col  # Column determined by input variable
    
    # Progress indicator every 10000 lines
    if(++processed % 10000 == 0) {
        printf "Processed %d lines...\n", processed > "/dev/stderr"
    }
    
    # Count cluster size and target species presence (include S and H, exclude C)
    if(type!="C") {
        cluster_count[cluster]++
        if(target_species == 1) {
            cluster_target_count[cluster]++
        }
    }
    
    # Store sequence only if type is S
    if(type=="S") {
        # Clean sequence: remove any character not A,C,G,T (case-insensitive)
        gsub(/[^ACGTacgt]/,"",seq)
        cluster_seq[cluster]=seq
    }
}

END {
    printf "Finished reading input. Calculating and writing output files...\n" > "/dev/stderr"
    
    # Process each cluster and write immediately to avoid memory issues
    cluster_num = 0
    filtered_out = 0
    
    for(cluster in cluster_seq) {
        cluster_num++
        if(cluster_num % 1000 == 0) {
            printf "Processing cluster %d...\n", cluster_num > "/dev/stderr"
        }
        
        size = cluster_count[cluster]
        target_count = cluster_target_count[cluster] + 0  # +0 converts null to 0
        
        # Apply filter based on mode
        if(mode == "percent") {
            # Calculate percentage of target species presence
            target_percent = (size > 0) ? (target_count / size * 100) : 0
            
            # Skip clusters with <=cutoff% target species presence
            if(target_percent <= cutoff) {
                filtered_out++
                continue
            }
            
            header = ">cluster_num=" cluster " cluster_size=" size " " label "_pct=" sprintf("%.1f", target_percent)
        } else {
            # Filter by count
            if(target_count <= cutoff) {
                filtered_out++
                continue
            }
            
            # Calculate percentage for header
            target_percent = (size > 0) ? (target_count / size * 100) : 0
            header = ">cluster_num=" cluster " cluster_size=" size " " label "_hits=" target_count " " label "_pct=" sprintf("%.1f", target_percent)
        }
        
        sequence = cluster_seq[cluster]
        
        # Choose output file based on size ranges
        if(size > 5) {
            output_file = large
        } else if(size >= 3 && size <= 5) {
            output_file = medium
        } else {
            output_file = small
        }
        
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
    close(medium)
    close(large)
    
    printf "Output writing complete.\n" > "/dev/stderr"
    if(mode == "percent") {
        printf "Filtered out %d clusters with <=%g%% %s presence\n", filtered_out, cutoff, label > "/dev/stderr"
    } else {
        printf "Filtered out %d clusters with <=%g %s hits\n", filtered_out, cutoff, label > "/dev/stderr"
    }
}
' "$INPUT"

echo "Splitting complete. Checking output files..."
echo "Small clusters (<3): $(grep -c "^>" "$OUTPUT_SMALL" 2>/dev/null || echo 0) sequences"
echo "Medium clusters (3-5): $(grep -c "^>" "$OUTPUT_MEDIUM" 2>/dev/null || echo 0) sequences"
echo "Large clusters (>5): $(grep -c "^>" "$OUTPUT_LARGE" 2>/dev/null || echo 0) sequences"
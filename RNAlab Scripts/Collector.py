#!/usr/bin/env python3
"""
Script to filter FASTA sequences based on Saccharomyces presence in corresponding TSV data.
Keeps clusters with >50% Saccharomyces presence.
"""

import re
import argparse
from collections import defaultdict
from pathlib import Path


def parse_fasta(fasta_file):
    """Parse FASTA file and extract cluster information."""
    clusters = {}
    current_header = None
    current_sequence = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_header:
                    clusters[current_header] = ''.join(current_sequence)
                
                # Extract cluster number from header
                cluster_match = re.search(r'cluster_num=(\d+)', line)
                if cluster_match:
                    cluster_num = int(cluster_match.group(1))
                    current_header = (cluster_num, line)
                    current_sequence = []
                else:
                    print(f"Warning: Could not extract cluster number from: {line}")
                    current_header = None
                    current_sequence = []
            elif current_header:
                current_sequence.append(line)
        
        # Don't forget the last sequence
        if current_header:
            clusters[current_header] = ''.join(current_sequence)
    
    return clusters


def parse_tsv(tsv_file, target_organism_col=10):
    """Parse TSV file and count target organism presence per cluster."""
    cluster_stats = defaultdict(lambda: {'total': 0, 'target_organism': 0})
    debug_organism_found = 0
    debug_total_lines = 0
    debug_header_skipped = False
    header_columns = []
    
    with open(tsv_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
                
            parts = line.split('\t')
            
            # Skip header line (if it contains 'record_type' in first column)
            if line_num == 1 and len(parts) > 0 and parts[0] == 'record_type':
                debug_header_skipped = True
                header_columns = parts
                print(f"DEBUG: Header columns found: {', '.join(f'{i}:{col}' for i, col in enumerate(parts))}")
                if target_organism_col < len(parts):
                    print(f"DEBUG: Target organism column {target_organism_col}: '{parts[target_organism_col]}'")
                else:
                    print(f"WARNING: Target column {target_organism_col} not found in header (max: {len(parts)-1})")
                continue
                
            debug_total_lines += 1
            
            if len(parts) <= target_organism_col:  # Need enough columns to access target column
                continue
                
            # Skip rows where column 1 (0-indexed column 0) is "C"
            if parts[0] == 'C':
                continue
                
            try:
                cluster_num = int(parts[1])  # Column 2 (0-indexed column 1)
            except ValueError:
                continue
                
            cluster_stats[cluster_num]['total'] += 1
            
            # Check target organism column
            try:
                organism_value = int(parts[target_organism_col])
                if organism_value > 0:  # Any positive value indicates organism presence
                    cluster_stats[cluster_num]['target_organism'] += 1
                    debug_organism_found += 1
                    if debug_organism_found <= 5:  # Show first 5 matches for debugging
                        organism_name = header_columns[target_organism_col] if header_columns else f"column_{target_organism_col}"
                        print(f"DEBUG: Found {organism_name} ({organism_value}) in cluster {cluster_num}")
            except (ValueError, IndexError):
                # If we can't parse the target organism column, skip this entry
                continue
    
    organism_name = header_columns[target_organism_col] if header_columns and target_organism_col < len(header_columns) else f"column_{target_organism_col}"
    print(f"DEBUG: Header skipped: {debug_header_skipped}")
    print(f"DEBUG: Total data lines processed: {debug_total_lines}")
    print(f"DEBUG: Lines with {organism_name} > 0: {debug_organism_found}")
    return cluster_stats, organism_name


def calculate_organism_percentage(cluster_stats):
    """Calculate target organism percentage for each cluster."""
    percentages = {}
    for cluster_num, stats in cluster_stats.items():
        if stats['total'] > 0:
            percentage = (stats['target_organism'] / stats['total']) * 100
            percentages[cluster_num] = {
                'percentage': percentage,
                'organism_count': stats['target_organism'],
                'total_count': stats['total']
            }
    return percentages


def filter_clusters(fasta_clusters, cluster_percentages, threshold=50.0):
    """Filter FASTA clusters based on Saccharomyces percentage threshold."""
    filtered_clusters = {}
    
    for (cluster_num, header), sequence in fasta_clusters.items():
        if cluster_num in cluster_percentages:
            percentage = cluster_percentages[cluster_num]['percentage']
            if percentage > threshold:
                filtered_clusters[(cluster_num, header)] = {
                    'sequence': sequence,
                    'percentage': percentage,
                    'stats': cluster_percentages[cluster_num]
                }
    
    return filtered_clusters


def write_filtered_fasta(filtered_clusters, output_file, organism_name):
    """Write filtered clusters to output FASTA file."""
    with open(output_file, 'w') as f:
        for (cluster_num, header), data in filtered_clusters.items():
            # Add organism percentage info to header
            stats = data['stats']
            new_header = f"{header} {organism_name}_percentage={data['percentage']:.1f}% ({stats['organism_count']}/{stats['total_count']})"
            f.write(f"{new_header}\n")
            
            # Write sequence with line wrapping (80 characters per line)
            sequence = data['sequence']
            for i in range(0, len(sequence), 80):
                f.write(f"{sequence[i:i+80]}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Filter FASTA sequences based on organism presence in TSV data'
    )
    parser.add_argument('fasta_file', help='Input FASTA file')
    parser.add_argument('tsv_file', help='Input TSV file')
    parser.add_argument('-o', '--output', default='filtered_organism.fasta',
                        help='Output FASTA file (default: filtered_organism.fasta)')
    parser.add_argument('-t', '--threshold', type=float, default=50.0,
                        help='Minimum organism percentage threshold (default: 50.0)')
    parser.add_argument('-c', '--column', type=int, default=10,
                        help='Target organism column index (0-based). Default: 10 (Saccharomyces)')
    parser.add_argument('--organism', type=str, 
                        help='Target organism name (for column lookup by name instead of index)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print detailed statistics')
    
    args = parser.parse_args()
    
    # Validate input files
    if not Path(args.fasta_file).exists():
        print(f"Error: FASTA file '{args.fasta_file}' not found")
        return 1
    
    if not Path(args.tsv_file).exists():
        print(f"Error: TSV file '{args.tsv_file}' not found")
        return 1
    
    # Determine target column
    target_column = args.column
    if args.organism:
        # Read header to find organism column by name
        with open(args.tsv_file, 'r') as f:
            header_line = f.readline().strip()
            if header_line:
                header_parts = header_line.split('\t')
                try:
                    target_column = header_parts.index(args.organism)
                    print(f"Found organism '{args.organism}' at column {target_column}")
                except ValueError:
                    print(f"Error: Organism '{args.organism}' not found in header")
                    print(f"Available columns: {', '.join(header_parts)}")
                    return 1
    
    print("Parsing FASTA file...")
    fasta_clusters = parse_fasta(args.fasta_file)
    print(f"Found {len(fasta_clusters)} clusters in FASTA file")
    
    print("Parsing TSV file...")
    cluster_stats, organism_name = parse_tsv(args.tsv_file, target_column)
    print(f"Found {len(cluster_stats)} clusters in TSV file")
    
    print(f"Calculating {organism_name} percentages...")
    cluster_percentages = calculate_organism_percentage(cluster_stats)
    
    # Debug: Show some cluster statistics
    clusters_with_organism = {k: v for k, v in cluster_percentages.items() if v['organism_count'] > 0}
    print(f"DEBUG: Clusters with any {organism_name}: {len(clusters_with_organism)}")
    if clusters_with_organism:
        print(f"DEBUG: Top 5 clusters with {organism_name}:")
        for i, (cluster_num, stats) in enumerate(list(clusters_with_organism.items())[:5]):
            print(f"  Cluster {cluster_num}: {stats['percentage']:.1f}% ({stats['organism_count']}/{stats['total_count']})")
    
    # Show clusters that exist in both FASTA and TSV
    fasta_cluster_nums = {cluster_num for cluster_num, _ in fasta_clusters.keys()}
    overlap = fasta_cluster_nums.intersection(cluster_percentages.keys())
    print(f"DEBUG: FASTA clusters also in TSV: {len(overlap)}")
    if overlap:
        print(f"DEBUG: Sample overlapping clusters: {list(overlap)[:10]}")
    
    print(f"Filtering clusters with >{args.threshold}% {organism_name} presence...")
    filtered_clusters = filter_clusters(fasta_clusters, cluster_percentages, args.threshold)
    
    print(f"Writing {len(filtered_clusters)} filtered clusters to {args.output}")
    write_filtered_fasta(filtered_clusters, args.output, organism_name)
    
    if args.verbose:
        print("\nDetailed Statistics:")
        print("=" * 50)
        for (cluster_num, _), data in filtered_clusters.items():
            stats = data['stats']
            print(f"Cluster {cluster_num}: {data['percentage']:.1f}% "
                  f"({stats['organism_count']}/{stats['total_count']} {organism_name})")
    
    print(f"\nSummary:")
    print(f"- Input clusters: {len(fasta_clusters)}")
    print(f"- Clusters with TSV data: {len([c for c in fasta_clusters.keys() if c[0] in cluster_percentages])}")
    print(f"- Clusters passing filter (>{args.threshold}% {organism_name}): {len(filtered_clusters)}")
    print(f"- Output written to: {args.output}")
    
    return 0


if __name__ == "__main__":
    exit(main())
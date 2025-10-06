#!/usr/bin/env python3
"""
Script to remove Diamond hits from original FASTA dataset based on cluster numbers
"""

import re
import argparse
from typing import Set, List, Tuple

def parse_diamond_hits(diamond_file: str) -> Set[str]:
    """
    Parse Diamond output file and extract cluster numbers from hits
    
    Args:
        diamond_file: Path to Diamond output file
        
    Returns:
        Set of cluster numbers that had hits
    """
    hit_clusters = set()
    
    with open(diamond_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and line.startswith('cluster_num='):
                # Extract cluster number from first column
                cluster_num = line.split('\t')[0].replace('cluster_num=', '')
                hit_clusters.add(cluster_num)
    
    return hit_clusters

def filter_fasta(input_fasta: str, output_fasta: str, hit_clusters: Set[str]):
    """
    Filter FASTA file to remove sequences with cluster numbers that had Diamond hits
    
    Args:
        input_fasta: Path to input FASTA file
        output_fasta: Path to output filtered FASTA file
        hit_clusters: Set of cluster numbers to remove
    """
    sequences_kept = 0
    sequences_removed = 0
    
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        current_header = None
        current_sequence = []
        
        for line in infile:
            line = line.strip()
            
            if line.startswith('>'):
                # Process previous sequence if exists
                if current_header is not None:
                    if should_keep_sequence(current_header, hit_clusters):
                        outfile.write(current_header + '\n')
                        outfile.write('\n'.join(current_sequence) + '\n')
                        sequences_kept += 1
                    else:
                        sequences_removed += 1
                
                # Start new sequence
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)
        
        # Process final sequence
        if current_header is not None:
            if should_keep_sequence(current_header, hit_clusters):
                outfile.write(current_header + '\n')
                outfile.write('\n'.join(current_sequence) + '\n')
                sequences_kept += 1
            else:
                sequences_removed += 1
    
    print(f"Filtering complete:")
    print(f"  Sequences kept: {sequences_kept}")
    print(f"  Sequences removed: {sequences_removed}")
    print(f"  Total processed: {sequences_kept + sequences_removed}")

def should_keep_sequence(header: str, hit_clusters: Set[str]) -> bool:
    """
    Determine if a sequence should be kept based on its cluster number
    
    Args:
        header: FASTA header line
        hit_clusters: Set of cluster numbers to remove
        
    Returns:
        True if sequence should be kept, False if it should be removed
    """
    # Extract cluster number from header using regex
    cluster_match = re.search(r'cluster_num=(\d+)', header)
    if cluster_match:
        cluster_num = cluster_match.group(1)
        return cluster_num not in hit_clusters
    else:
        # If no cluster number found, keep the sequence by default
        print(f"Warning: No cluster number found in header: {header}")
        return True

def main():
    parser = argparse.ArgumentParser(description='Remove Diamond hits from FASTA dataset')
    parser.add_argument('diamond_file', help='Diamond output file')
    parser.add_argument('input_fasta', help='Input FASTA file')
    parser.add_argument('output_fasta', help='Output filtered FASTA file')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    print("Parsing Diamond hits...")
    hit_clusters = parse_diamond_hits(args.diamond_file)
    print(f"Found {len(hit_clusters)} unique cluster numbers with hits")
    
    if args.verbose:
        print("Cluster numbers to remove:")
        for cluster in sorted(hit_clusters, key=int):
            print(f"  cluster_num={cluster}")
    
    print("\nFiltering FASTA file...")
    filter_fasta(args.input_fasta, args.output_fasta, hit_clusters)

if __name__ == "__main__":
    main()

# Example usage:
# python Diamond.breaker.py diamond_hits.txt input.fasta filtered_output.fasta

# You can also use the functions directly:
"""
# Parse diamond hits
hit_clusters = parse_diamond_hits('diamond_hits.txt')

# Filter FASTA
filter_fasta('original_dataset.fasta', 'filtered_dataset.fasta', hit_clusters)
"""

#!/usr/bin/env python3
"""
FASTA-TSV Sequence Merger
Merges FASTA sequences into TSV file based on header matching with SRR column.
"""

import argparse
import sys
import csv
from collections import defaultdict

def parse_fasta(fasta_file):
    """
    Parse FASTA file and return a dictionary mapping header keys to sequences.
    Key is extracted as everything up to the second underscore.
    """
    sequences = {}
    current_header = None
    current_sequence = []
    warning_count = 0
    max_warnings = 10
    
    print(f"Reading FASTA file: {fasta_file}")
    
    with open(fasta_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
            
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_header and current_sequence:
                    key = extract_key_from_header(current_header)
                    if key:
                        sequences[key] = ''.join(current_sequence)
                    else:
                        print(f"Warning: Could not extract key from header '{current_header}' at line {line_num}")
                
                # Start new sequence
                current_header = line[1:]  # Remove '>'
                current_sequence = []
                
            elif line and current_header:
                # Valid sequence line
                current_sequence.append(line)
            elif line and not current_header:
                # Sequence data without header - problematic
                if warning_count < max_warnings:
                    print(f"Warning: Sequence data found before header at line {line_num}: '{line[:50]}{'...' if len(line) > 50 else ''}'")
                elif warning_count == max_warnings:
                    print(f"Warning: Too many orphaned sequence lines. Suppressing further warnings...")
                warning_count += 1
    
    # Save last sequence
    if current_header and current_sequence:
        key = extract_key_from_header(current_header)
        if key:
            sequences[key] = ''.join(current_sequence)
        else:
            print(f"Warning: Could not extract key from final header '{current_header}'")
    
    print(f"Loaded {len(sequences)} sequences from FASTA")
    if warning_count > 0:
        print(f"Total orphaned sequence lines: {warning_count}")
    
    return sequences

def debug_fasta_lines(fasta_file, center_line, context=10):
    """
    Debug function to show lines around a specific line number in FASTA file.
    """
    print(f"Debugging FASTA file '{fasta_file}' around line {center_line}")
    print(f"Showing {context} lines before and after:")
    print("-" * 80)
    
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
    
    start_line = max(0, center_line - context - 1)
    end_line = min(len(lines), center_line + context)
    
    for i in range(start_line, end_line):
        line_num = i + 1
        line = lines[i].rstrip('\n\r')
        marker = " >>> " if line_num == center_line else "     "
        print(f"{marker}{line_num:6}: {repr(line)}")
    
    print("-" * 80)

def extract_key_from_header(header):
    """
    Extract key from FASTA header (everything up to second underscore).
    For headers like 'DRR220096_248502_circle_248502_1', extract 'DRR220096_248502'
    """
    parts = header.split('_')
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    return None

def extract_key_from_query_label(query_label):
    """
    Extract key from TSV query_label.
    For query_label like 'DRR220096_248502_circle_248502_1 [103 - 375] ...', 
    extract 'DRR220096_248502'
    """
    # Take the first part before any space
    first_part = query_label.split()[0] if query_label else ""
    parts = first_part.split('_')
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    return None

def merge_data(tsv_file, fasta_sequences, output_file, delimiter='\t', srr_column='SRR', verbose=False):
    """
    Merge TSV data with FASTA sequences based on query_label matching.
    Only processes rows where record_type = 'S'.
    """
    print(f"Reading TSV file: {tsv_file}")
    
    matches_found = 0
    total_rows = 0
    s_rows_processed = 0
    
    with open(tsv_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        # Use csv.Sniffer to detect delimiter if not specified
        if delimiter == 'auto':
            sample = infile.read(1024)
            infile.seek(0)
            sniffer = csv.Sniffer()
            delimiter = sniffer.sniff(sample).delimiter
            print(f"Auto-detected delimiter: '{delimiter}'")
        
        reader = csv.DictReader(infile, delimiter=delimiter)
        
        # Strip whitespace from fieldnames and create mapping
        original_fieldnames = reader.fieldnames
        stripped_fieldnames = [name.strip() for name in original_fieldnames]
        
        # Check if required columns exist
        required_columns = ['record_type', 'query_label']
        for req_col in required_columns:
            if req_col not in stripped_fieldnames:
                print(f"Error: Required column '{req_col}' not found in TSV file.")
                print("Available columns: " + ', '.join([f'"{name.strip()}"' for name in original_fieldnames]))
                sys.exit(1)
        
        # Find column mappings
        record_type_col = None
        query_label_col = None
        for orig_name, stripped_name in zip(original_fieldnames, stripped_fieldnames):
            if stripped_name == 'record_type':
                record_type_col = orig_name
            elif stripped_name == 'query_label':
                query_label_col = orig_name
        
        # Add 'contig' to fieldnames
        fieldnames = reader.fieldnames + ['contig']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        
        for row in reader:
            total_rows += 1
            
            # Only process rows where record_type = 'S'
            record_type = row.get(record_type_col, '').strip()
            if record_type != 'S':
                row['contig'] = ''
                writer.writerow(row)
                continue
            
            s_rows_processed += 1
            query_label = row.get(query_label_col, '').strip()
            
            # Extract key from query_label
            key = extract_key_from_query_label(query_label)
            
            # Add contig sequence if match found
            if key and key in fasta_sequences:
                row['contig'] = fasta_sequences[key]
                matches_found += 1
                if verbose and matches_found <= 3:
                    print(f"  Match {matches_found}: TSV key '{key}' -> FASTA sequence found")
            else:
                row['contig'] = ''
                if verbose and s_rows_processed <= 5:
                    print(f"  No match: TSV key '{key}' not found in FASTA")
            
            writer.writerow(row)
    
    print(f"Processing complete!")
    print(f"Total TSV rows processed: {total_rows}")
    print(f"Rows with record_type 'S': {s_rows_processed}")
    print(f"Matches found: {matches_found}")
    print(f"Output written to: {output_file}")
    
    return matches_found, total_rows

def main():
    parser = argparse.ArgumentParser(
        description="Merge FASTA sequences into TSV file based on header matching",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 merge_fasta_tsv.py input.tsv sequences.fasta -o merged_output.tsv
  python3 merge_fasta_tsv.py input.tsv sequences.fasta -o merged_output.tsv -d ","
  python3 merge_fasta_tsv.py input.tsv sequences.fasta -o merged_output.tsv --delimiter auto

The script matches FASTA header prefixes (up to 2nd underscore) with the SRR column in TSV.
Example: '>DRR000111_114_circle_114' matches SRR column value 'DRR000111_114'
        """
    )
    
    parser.add_argument('tsv_file', help='Input TSV file with SRR column')
    parser.add_argument('fasta_file', help='Input FASTA file with sequences')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file name')
    parser.add_argument('-d', '--delimiter', default='\t', 
                       help='TSV delimiter (default: tab). Use "auto" to auto-detect')
    parser.add_argument('-v', '--verbose', action='store_true', 
                       help='Verbose output showing sample matches')
    parser.add_argument('--debug-fasta', type=int, metavar='LINE_NUM',
                       help='Show lines around specified line number in FASTA for debugging')
    parser.add_argument('--srr-column', default='SRR', 
                       help='Column name containing SRR identifiers (default: SRR)')
    
    args = parser.parse_args()
    
    try:
        # Debug FASTA file if requested
        if args.debug_fasta:
            debug_fasta_lines(args.fasta_file, args.debug_fasta)
            return
        
        # Parse FASTA file
        fasta_sequences = parse_fasta(args.fasta_file)
        
        if not fasta_sequences:
            print("Error: No sequences found in FASTA file")
            sys.exit(1)
        
            if verbose:
                print("\nSample FASTA keys (first 5):")
                for i, key in enumerate(list(fasta_sequences.keys())[:5]):
                    seq_preview = fasta_sequences[key][:50] + "..." if len(fasta_sequences[key]) > 50 else fasta_sequences[key]
                    print(f"  {key} -> {seq_preview}")
                
                # Show some SRR/ERR/DRR breakdown
                srr_count = sum(1 for k in fasta_sequences.keys() if k.startswith('SRR'))
                drr_count = sum(1 for k in fasta_sequences.keys() if k.startswith('DRR'))
                err_count = sum(1 for k in fasta_sequences.keys() if k.startswith('ERR'))
                other_count = len(fasta_sequences) - srr_count - drr_count - err_count
                print(f"\nFASTA key breakdown:")
                print(f"  SRR: {srr_count}")
                print(f"  DRR: {drr_count}")
                print(f"  ERR: {err_count}")
                print(f"  Other: {other_count}")
                print()
        
        # Merge data
        matches, total = merge_data(args.tsv_file, fasta_sequences, args.output, args.delimiter, args.srr_column, args.verbose)
        
        if matches == 0:
            print("\nWarning: No matches found! Please check:")
            print("1. That the SRR column exists in your TSV file")
            print("2. That FASTA headers match the expected format")
            print("3. That there are overlapping values between SRR column and FASTA headers")
            
            if args.verbose:
                print(f"\nFirst few FASTA keys: {list(fasta_sequences.keys())[:10]}")
                
                # Show sample TSV SRR values for comparison
                print(f"\nSample TSV SRR values:")
                with open(args.tsv_file, 'r') as f:
                    reader = csv.DictReader(f, delimiter=args.delimiter if args.delimiter != 'auto' else '\t')
                    original_fieldnames = reader.fieldnames
                    srr_column_name = None
                    for orig_name in original_fieldnames:
                        if orig_name.strip() == srr_column:
                            srr_column_name = orig_name
                            break
                    
                    sample_srr_values = []
                    for i, row in enumerate(reader):
                        if i >= 5:
                            break
                        if srr_column_name:
                            sample_srr_values.append(row.get(srr_column_name, '').strip())
                    
                    for val in sample_srr_values:
                        print(f"  TSV SRR: '{val}'")
                        if val in fasta_sequences:
                            print(f"    ✓ Found in FASTA")
                        else:
                            print(f"    ✗ NOT found in FASTA")
        
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

import csv
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Convert a CSV to FASTA")
parser.add_argument("-i", "--input", required=True, help="Path to input CSV file")
parser.add_argument("-o", "--output", required=True, help="Path to output FASTA file")
args = parser.parse_args()

# Read CSV and write FASTA
with open(args.input, "r") as csvfile, open(args.output, "w") as fasta:
    reader = csv.DictReader(csvfile)
    for row in reader:
        fasta.write(f">{row['Contig_ID']}\n")
        fasta.write(f"{row['Sequence']}\n")

print(f"FASTA file saved as {args.output}")

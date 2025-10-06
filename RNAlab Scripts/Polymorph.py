#!/usr/bin/env python3
import argparse
import csv

def parse_uc(uc_file, tsv_file):
    header = [
        "RecordType", "Cluster", "Length", "PctId", "Strand",
        "Mismatch", "GapOpen", "Qlo", "Qhi", "Tlo", "Thi",
        "Evalue", "BitScore", "Query", "Target"
    ]

    with open(uc_file, "r") as infile, open(tsv_file, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(header)
        for line in infile:
            if line.startswith("#") or not line.strip():
                continue  # skip comments and blanks
            fields = line.strip().split("\t")
            # pad to header length if shorter
            while len(fields) < len(header):
                fields.append("")
            writer.writerow(fields[:len(header)])

def main():
    parser = argparse.ArgumentParser(description="Convert .uc file to .tsv")
    parser.add_argument("--input", "-i", required=True, help="Input .uc file")
    parser.add_argument("--output", "-o", required=True, help="Output .tsv file")
    args = parser.parse_args()

    parse_uc(args.input, args.output)

if __name__ == "__main__":
    main()

import argparse
import gzip
import os
import sys


# Authors: Betty Liu, Ben Parks
# Last updated: 11/21/22

# Convert a tab-delimited counts matrix (gene, cell, count) from stdin 
#  to a matrix market file format of (row, col, count). Does not include a
#  header in the output

# Usage: python mtx_from_counts.py counts.tsv.gz features.tsv.gz barcodes.tsv.gz > matrix.mtx 
# Inputs:
# - Counts: stdin counts tsv with columns gene_id, cell_id, umi_count
# - Features: TSV of gene IDs, where gene_id is the first column
# - Barcodes: TSV of cell barcodes, where cell_id is the first column
# Outputs:
# - mtx: TSV of (row, col, count) written to stdout

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("features", type=str, help="Path features.tsv file")
    parser.add_argument("barcodes", type=str, help="Path features.tsv file")
    args = parser.parse_args()

    # Load index corresponding to each feature name + cell name
    feature_dict = {}
    for idx, line in enumerate(gzip.open(args.features, "rb")):
        id = line.strip().split(b"\t")[0]
        feature_dict[id] = str(idx+1).encode()

    barcode_dict = {}
    for idx, line in enumerate(gzip.open(args.barcodes, "rb")):
        id = line.strip().split(b"\t")[0]
        barcode_dict[id] = str(idx+1).encode()

    for line in sys.stdin.buffer:
        gene, cell, umi_count = line.split(b"\t")
        sys.stdout.buffer.write(feature_dict[gene])
        sys.stdout.buffer.write(b"\t")
        sys.stdout.buffer.write(barcode_dict[cell])
        sys.stdout.buffer.write(b"\t")
        sys.stdout.buffer.write(umi_count) # umi_count has a newline on the end of it

if __name__ == "__main__":
    main()
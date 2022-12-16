import argparse
import gzip
import os
import sys

# Convert a tab-delimited counts matrix (gene, cell, count) from stdin 
#  to a matrix market file format. Must provide pre-saved features file. 

# Authors: Betty Liu, Ben Parks
# Last updated: 12/5/22

# Usage: python counts_to_mtx.py counts.tsv.gz features.tsv.gz --mtx matrix.mtx.gz --barcodes barcodes.tsv.gz
# Inputs:
# - Counts: Gzipped TSV of counts with columns gene_id, cell_id, umi_count
# - Features: Gzipped TSV of gene IDs, where gene_id is the first column

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("counts", type=str, help="Path to counts.tsv.gz file")
    parser.add_argument("features", type=str, help="Path to features.tsv.gz file")
    parser.add_argument("--mtx", type=str, help="Path of gzipped mtx output", default="matrix.mtx.gz")
    parser.add_argument("--barcodes", type=str, help="Path of gzipped barcodes output", default="barcodes.tsv.gz")
    args = parser.parse_args()

    # Load index corresponding to each feature name
    feature_dict = {}
    for idx, line in enumerate(gzip.open(args.features, "rb")):
        gene_id = line[:line.find(b"\t")]
        feature_dict[gene_id] = str(idx+1).encode()
    

    # Get the set of cell barcodes and non-zero count in first pass through the input
    next_barcode = 0
    barcode_dict = {}
    barcodes = gzip.open(args.barcodes, "wb")

    nonzero_count = 0
    for line in gzip.open(args.counts):
        nonzero_count += 1
        gene, cell, umi_count = line.split(b"\t")
        if cell not in barcode_dict:
            barcode_dict[cell] = str(next_barcode+1).encode()
            barcodes.write(cell)
            barcodes.write(b"\n")
            next_barcode += 1

    # Write header
    mtx = gzip.open(args.mtx, "wb")
    mtx.write(b"%%MatrixMarket matrix coordinate real general\n")
    mtx.write(b"%\n")
    mtx.write(f"{len(feature_dict)}\t{len(barcode_dict)}\t{nonzero_count}\n".encode())

    for line in gzip.open(args.counts, "rb"):
        gene, cell, umi_count = line.split(b"\t")
        mtx.write(feature_dict[gene])
        mtx.write(b"\t")
        mtx.write(barcode_dict[cell])
        mtx.write(b"\t")
        mtx.write(umi_count) # umi_count has a newline on the end of it


if __name__ == "__main__":
    main()

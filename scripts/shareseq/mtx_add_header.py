import argparse
import gzip
import shutil
import sys

# Authors: Ben Parks
# Last updated: 11/21/22

# Add mtx header to a tab-delimited file of (row, col, value) entries

# Usage: python mtx_merge.py input.mtx features.tsv.gz barcodes.tsv.gz > output.mtx
#  python mtx_merge.py counts.tsv.gz features.tsv --mtx matrix.mtx --barcodes barcodes.tv
# Inputs:
# - mtx: (row, col, value) tab-separated mtx entries (no header)
# - features: TSV of gene IDs, where gene_id is the first column
# - barcodes: TSV of cell barcodes, where cell_id is the first column
# Outputs:
# - mtx: TSV of (row, col, count) written to stdout with appropriate Matrix-Market header

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mtx", type=str, help="Path to mtx.gz file")
    parser.add_argument("features", type=str, help="Path features.tsv file")
    parser.add_argument("barcodes", type=str, help="Path barcodes.tsv file")
    args = parser.parse_args()

    feature_count = sum(1 for line in gzip.open(args.features, "rb"))
    barcode_count = sum(1 for line in gzip.open(args.barcodes, "rb"))

    # Slow step
    nonzero_count = sum(1 for line in gzip.open(args.mtx, "rb"))

    # Write header
    sys.stdout.buffer.write(b"%%MatrixMarket matrix coordinate real general\n")
    sys.stdout.buffer.write(b"%\n")
    sys.stdout.buffer.write(f"{feature_count}\t{barcode_count}\t{nonzero_count}\n".encode())

    # Write contents
    shutil.copyfileobj(gzip.open(args.mtx, "rb"), sys.stdout.buffer)

if __name__ == "__main__":
    main()

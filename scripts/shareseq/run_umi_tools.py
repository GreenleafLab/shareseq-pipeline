# Authors: Ben Parks
# Last updated: 11/21/22

# Run UMI tools with custom I/O formats
# Input: tsv of gene, cell, umi, count. Sorted by gene, cell, umi
# Output: tsv of gene, cell, count. Deduplicated by UMI with 1bp mismatch
# Read from stdin, write to stdout. Writes a one-line log at the end to stderr

# Usage: python run_umi_tools.py < input > output 2> log

import argparse
import sys

import umi_tools

def print_line(gene, cell, count):
    sys.stdout.buffer.write(gene)
    sys.stdout.buffer.write(b"\t")
    sys.stdout.buffer.write(cell)
    sys.stdout.buffer.write(b"\t")
    sys.stdout.buffer.write(str(count).encode())
    sys.stdout.buffer.write(b"\n")

total_reads = 0
clusterer = umi_tools.UMIClusterer(cluster_method="directional")

prev_gene = None
prev_cell = None
umi_dict = {}

for line in sys.stdin.buffer:
    gene, cell, umi, count = line.split(b"\t")
    if gene != prev_gene or cell != prev_cell:
        if prev_gene is not None:
            dedup_count = len(clusterer(umi_dict, threshold=1))
            print_line(prev_gene, prev_cell, dedup_count)
            total_reads += dedup_count
        umi_dict.clear()
        prev_gene, prev_cell = gene, cell
    umi_dict[umi] = int(count)
print_line(prev_gene, prev_cell, dedup_count)
total_reads += dedup_count

print(f"INFO Number of (post deduplication) reads counted: {total_reads}", file=sys.stderr)
# Authors: Ben Parks
# Last updated: 11/21/22

# Convert featureCounts CORE-format output to a tsv with gene, cell, umi
# Read from stdin, write to stdout
import sys

for line in sys.stdin:
    read, status, count, feature = line.strip().split("\t")
    if status != "Assigned":
        continue
    _, cell, umi = read.rsplit("_", maxsplit=2)
    print(feature, cell, umi, sep="\t")
    
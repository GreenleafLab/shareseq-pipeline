
# Merge duplicate fragments sorted by (chr, start, end, cell) 
# Fragment should be format chr, start, end, cell, count (where count can be 1)
# In output, the counts of adjacent of (chr, start, end, cell) will be added together
# so only one line is output per unique (chr, start, end, cell) combination

# Authors: Ben Parks
# Last updated: 9/27/22

# Usage: python reverse_comp_barcodes.py input.tsv output.tsv
import argparse

from typing import *

def main():
    parser = argparse.ArgumentParser(description="Reverse complement the barcode file")
    parser.add_argument("input", type=str, help="input path")
    parser.add_argument("output", type=str, help="output path")
    args = parser.parse_args()
    input_lines = open(args.input, "rb").readlines()
    
    output = open(args.output, "wb")

    assert input_lines[0] == b"Name\tSequence\n"
    output.write(input_lines[0])
    for l in input_lines[1:]:
        name, seq = l.strip().split(b"\t")
        seq = reverse_complement(seq.upper())
        output.write(name + b"\t" + seq + b"\n")

complement = bytes.maketrans(b"ATGC", b"TACG")
def reverse_complement(seq):
    return seq.translate(complement)[::-1]


if __name__ == "__main__":
    main()

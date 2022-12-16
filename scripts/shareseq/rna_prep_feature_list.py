import argparse
import re
import sys

# Generates a list of features from gtf (stdin)

# Authors: Betty Liu, Ben Parks
# Last updated: 11/17/22

# Usage: python rna_prep_feature_list.py gencode.v41.annotations.gt > features.tsv

#gtf_colnames = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"] # for reference

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gtf", type=str, help="Path of input gtf file")
    parser.add_argument("--feature", type=str, help="feature type to output, default to gene", default="gene")
    parser.add_argument("--attribute", type=str, nargs="+", help="attribute columns to output, default to gene_id, gene_name", default=["gene_id", "gene_name"])
    parser.add_argument("--nomito", action="store_true", help="filter out mitochondrial genes")
    args = parser.parse_args()

    attribute_patterns = [re.compile(f'{attribute} "([^"]*)"') for attribute in args.attribute]

    for line in open(args.gtf, "r"):
        if line.startswith("#"): # skip comment lines
            continue
        
        line_ls = line.strip().split('\t')
        if args.nomito and line_ls[0] == "chrM":
            continue
        if line_ls[2] != args.feature:
            continue
        
        attributes = [attr.search(line_ls[-1]).group(1) for attr in attribute_patterns]
        print("\t".join(attributes))       

if __name__ == "__main__":
    main()

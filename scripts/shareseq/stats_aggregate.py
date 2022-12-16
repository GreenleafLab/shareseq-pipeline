import argparse
import glob
import gzip
import json
import re
import sys

def main():
    parser = argparse.ArgumentParser(description="Collect statistics from ATAC pipeline")
    parser.add_argument("--input", type=str, nargs="+", help="Paths of input json summary files")
    parser.add_argument("--output", type=str, help="Path of output json summary file")
    args = parser.parse_args()

    stats = {}
    for f in args.input:
        stats = add_stat_counts(stats, json.load(open(f)))
    
    ### Output summary
    json.dump(stats, open(args.output, "w"), indent = 2)

def add_stat_counts(d1, d2):
    """
    Merge to stat dictionaries together, where all corresponding integer
    values are added
    """
    out = {}
    for k, v in d1.items():
        if k not in d2:
            out[k] = v
        elif isinstance(v, int) or isinstance(v, float):
            out[k] = v + d2[k]
        elif isinstance(v, dict):
            out[k] = add_stat_counts(v, d2[k])
        elif isinstance(v, list):
            assert len(d2[k]) == len(v)
            out[k] = v
            for i in range(len(v)):
                v[i] += d2[k][i]
    for k, v in d2.items():
        if k not in d1:
            out[k] = v
    return out


if __name__ == "__main__":
    main()

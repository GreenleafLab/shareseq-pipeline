import argparse
import glob
import gzip
import json
import re
import sys

def main():
    parser = argparse.ArgumentParser(description="Collect statistics from ATAC pipeline")
    parser.add_argument("--barcode_stats", type=str, nargs="+", help="Paths of input barcode matching stats")
    parser.add_argument("--bowtie2_log", type=str, nargs="+", help="Paths of bowtie2 logs")
    parser.add_argument("--fragment_file", type=str, help="Path of fragment file (uncompressed)")
    parser.add_argument("--summary_output", type=str, help="Path of output json summary file")
    parser.add_argument("--barcodes_output", type=str, help="Path of detailed barcode matching summary file")
    args = parser.parse_args()

    ### Barcode-matching stats
    barcode_stats = {}
    for f in args.barcode_stats:
        stats = json.load(open(f))
        barcode_stats = add_stat_counts(barcode_stats, stats)

    total_reads = sum(barcode_stats["total_mismatch_histogram"])
    matching_reads = sum(barcode_stats["total_mismatch_histogram"][:-1])

    ### Bowtie2 stats
    bowtie2_stats = {}
    for f in args.bowtie2_log:
        bowtie2_stats = add_stat_counts(bowtie2_stats, parse_bowtie2_log(f))

    aligned_reads = (bowtie2_stats["concordant_1"] + bowtie2_stats["concordant_multi"] + bowtie2_stats["discordant_1"] + \
                     0.5*bowtie2_stats["aligned_mates_1"] + 0.5*bowtie2_stats["aligned_mates_multi"])
    
    ### Fragment duplication stats (slow part -- scans through fragments.tsv.gz)
    total_fragments = 0
    unique_fragments = 0
    chrM_fragments = 0
    for l in open(args.fragment_file, "rb"):
        total_fragments += int(l[l.rfind(b"\t"):])
        unique_fragments += 1
        if l[:l.find(b"\t")] == b"chrM":
            chrM_fragments += 1

    ### Output summary
    json.dump({
        "total_reads": total_reads,
        "valid_barcode_reads": matching_reads,
        "total_aligned_reads": aligned_reads,
        "multimapped_reads": bowtie2_stats["concordant_multi"],
        "total_fragments": total_fragments,
        "unique_fragments": unique_fragments,
        "mitochondrial_fragments": chrM_fragments
    }, open(args.summary_output, "w"), indent = 2)

    json.dump(barcode_stats, open(args.barcodes_output, "w"), indent=2)

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

#############################
### Bowtie2 parsing
#############################

percentage_regex = "\\(\d+\.\d+%\\)"
bowtie_log_template = """(?P<total_reads>\d+) reads; of these:
  (?P<paired_reads>\d+) {percentage_regex} were paired; of these:
    (?P<concordant_0>\d+) {percentage_regex} aligned concordantly 0 times
    (?P<concordant_1>\d+) {percentage_regex} aligned concordantly exactly 1 time
    (?P<concordant_multi>\d+) {percentage_regex} aligned concordantly >1 times
    ----
    (?P<concordant_0_b>\d+) pairs aligned concordantly 0 times; of these:
      (?P<discordant_1>\d+) {percentage_regex} aligned discordantly 1 time
    ----
    (?P<unaligned_pairs>\d+) pairs aligned 0 times concordantly or discordantly; of these:
      (?P<unaligned_mates>\d+) mates make up the pairs; of these:
        (?P<aligned_mates_0>\d+) {percentage_regex} aligned 0 times
        (?P<aligned_mates_1>\d+) {percentage_regex} aligned exactly 1 time
        (?P<aligned_mates_multi>\d+) {percentage_regex} aligned >1 times
\d+\\.\d+% overall alignment rate""".format(percentage_regex=percentage_regex)

def parse_bowtie2_log(path):
    output_groups = ["total_reads", "concordant_0", "concordant_1", "concordant_multi", "discordant_1", "aligned_mates_1", "aligned_mates_multi"]
    log = open(path).read()
    match = re.search(bowtie_log_template, log)
    if match is None:
        raise Exception(f"Path {path} doesn't match template regex for bowtie2 log format")
    return {group: int(match.group(group)) for group in output_groups}


if __name__ == "__main__":
    main()

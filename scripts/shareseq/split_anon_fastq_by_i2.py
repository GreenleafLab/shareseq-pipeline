#!/usr/bin/env python3
"""Split a per-sample anonymized SHARE-seq FASTQ pair into per-sublibrary
.fastq.zst files, demultiplexing by the I2 index in the read header.

Input headers (produced by anonymize.smk -> merge_samples -> SRA upload, then
fastq-dump/fasterq-dump with `--defline-seq '@$sn 1:N:0:$sg'` /
`--seq-defline '@$sn 1:N:0:$sg'`) look like:

    @A00509:707:H5HMLDSX7:1:1113:1000:11710 1:N:0:<99bp_I1>+GTTATCGT

The last 8 bp after the final `+` is the I2 sublibrary index. Each read pair
is routed to the sublibrary whose I2 matches up to 1bp (mirroring the
1-mismatch tolerance bcl2fastq applies for I2 in the non-anon pipeline, and
the per-barcode 1bp correction match_barcodes.py applies for BC1+BC2+BC3).
1bp variants that are equidistant from two canonical I2s are dropped to
avoid mis-assignment.

Output naming: {out_dir}/{sublibrary}_{R1,R2}.fastq.zst -- one pair of
streams per sublibrary plus optional unmatched bucket. Reads land in the same
TAB-as-space-normalized form ingest_anonymized.smk's old concatenate rule
produced, so downstream shareseq.smk match_barcodes is unaffected.
"""

import argparse
import itertools
import json
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path


def open_reader(path, threads):
    """Open a .fastq.gz file via pigz/gzip subprocess for parallel decompression."""
    decompressor = "pigz" if shutil.which("pigz") else "gzip"
    cmd = [decompressor, "-dc", str(path)] if decompressor == "gzip" else [decompressor, "-dcp", str(threads), str(path)]
    return subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=4 * 1024 * 1024)


def open_writer(path, threads_per_stream):
    """Open a zstd subprocess that writes to `path`. Returns Popen."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    cmd = ["zstd", "--fast=1", "-q", f"-T{threads_per_stream}", "-o", str(path)]
    return subprocess.Popen(cmd, stdin=subprocess.PIPE, bufsize=4 * 1024 * 1024)


def single_mismatches(seq):
    """Yield all 1bp variants of seq, including N substitutions.

    Mirrors single_mismatches() in scripts/shareseq/match_barcodes.py so the
    I2 demultiplexing uses the same convention as BC1+BC2+BC3 matching.
    """
    for idx, base in itertools.product(range(len(seq)), (b"A", b"T", b"G", b"C", b"N")):
        if base == seq[idx:idx + 1]:
            continue
        yield seq[:idx] + base + seq[idx + 1:]


def build_i2_lookup(i2_map):
    """Pre-compute a lookup {variant_bytes -> (sublibrary_str, mismatches_int)}.

    Includes each canonical I2 (mismatches=0) and all 1bp variants (mismatches=1),
    with ambiguous variants (1bp from two canonicals) dropped. Mirrors
    add_mismatches() in scripts/shareseq/match_barcodes.py.
    """
    canonical = {v.encode(): k for k, v in i2_map.items()}
    if len(canonical) != len(i2_map):
        sys.exit(f"error: --i2-map has duplicate I2 sequences: {i2_map}")
    out = {}
    for seq, name in canonical.items():
        out[seq] = (name, 0)
        for mutant in single_mismatches(seq):
            if mutant in canonical:
                # Two canonical I2s are pairwise within Hamming 1 - shouldn't
                # happen for SHARE-seq batch indices but guard anyway.
                sys.exit(f"error: canonical I2 sequences are within Hamming 1: {canonical}")
            if mutant in out:
                out[mutant] = None
            else:
                out[mutant] = (name, 1)
    return {k: v for k, v in out.items() if v is not None}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r1-in", required=True, help="merged anon R1 fastq.gz")
    ap.add_argument("--r2-in", required=True, help="merged anon R2 fastq.gz")
    ap.add_argument("--i2-map", required=True,
                    help="JSON object mapping {sublibrary_name: I2_sequence}, e.g. '{\"CL65\":\"CTTAATGC\",...}'")
    ap.add_argument("--out-dir", required=True, help="directory for per-sublibrary fastq.zst output")
    ap.add_argument("--stats", required=True, help="output JSON path with per-sublibrary read counts")
    ap.add_argument("--threads", type=int, default=4, help="threads for input pigz decompression")
    args = ap.parse_args()

    i2_map = json.loads(args.i2_map)
    if not i2_map:
        sys.exit("error: --i2-map is empty")
    # Sanity-check: I2 sequences must be unique and same length
    i2_values = list(i2_map.values())
    if len(set(i2_values)) != len(i2_values):
        sys.exit(f"error: --i2-map has duplicate I2 sequences: {i2_map}")
    i2_lens = {len(v) for v in i2_values}
    if len(i2_lens) != 1:
        sys.exit(f"error: --i2-map I2 sequences have inconsistent lengths: {i2_lens}")
    i2_len = i2_lens.pop()
    # Build lookup including 1bp variants (with collision-drop for ambiguous variants)
    rev = build_i2_lookup(i2_map)

    out_dir = Path(args.out_dir)
    # 1 thread per zstd is enough; pigz on input gets the rest
    writers = {
        sublib: (open_writer(out_dir / f"{sublib}_R1.fastq.zst", 1),
                 open_writer(out_dir / f"{sublib}_R2.fastq.zst", 1))
        for sublib in i2_map
    }

    r1_proc = open_reader(args.r1_in, args.threads)
    r2_proc = open_reader(args.r2_in, args.threads)
    r1 = r1_proc.stdout
    r2 = r2_proc.stdout

    exact_counts = Counter()
    mismatch_counts = Counter()
    unmatched = 0
    line_no = 0

    try:
        while True:
            # Read one 4-line record from each
            h1 = r1.readline()
            if not h1:
                # confirm R2 also exhausted
                if r2.readline():
                    sys.exit("error: R2 has more reads than R1 (length mismatch)")
                break
            s1 = r1.readline()
            p1 = r1.readline()
            q1 = r1.readline()
            h2 = r2.readline()
            s2 = r2.readline()
            p2 = r2.readline()
            q2 = r2.readline()
            if not (s1 and p1 and q1 and h2 and s2 and p2 and q2):
                sys.exit(f"error: truncated record near R1 line {line_no + 1}")
            line_no += 4

            # Parse I2 from the trailing bytes of the R1 header (after the last '+')
            # h1 ends with b"\n"; strip and find last '+'
            head = h1.rstrip(b"\r\n")
            plus_idx = head.rfind(b"+")
            if plus_idx == -1 or len(head) - plus_idx - 1 != i2_len:
                # malformed header; count as unmatched
                unmatched += 1
                continue
            i2 = head[plus_idx + 1:]
            hit = rev.get(i2)
            if hit is None:
                unmatched += 1
                continue
            sublib, mm = hit

            # Normalize the TAB-before-`1:N:0:` (samtools fastq artifact) to a space,
            # mirroring the awk pass in the old concatenate rule. We only need to do
            # this for the R1 header since match_barcodes parses indices from R1.
            # But to keep R2 byte-identical to the old pipeline's output, do it
            # there too. Performance hit is negligible (one tab.replace per header).
            h1_norm = h1.replace(b"\t", b" ", 1)
            h2_norm = h2.replace(b"\t", b" ", 1)

            w_r1, w_r2 = writers[sublib]
            w_r1.stdin.write(h1_norm)
            w_r1.stdin.write(s1)
            w_r1.stdin.write(p1)
            w_r1.stdin.write(q1)
            w_r2.stdin.write(h2_norm)
            w_r2.stdin.write(s2)
            w_r2.stdin.write(p2)
            w_r2.stdin.write(q2)
            if mm == 0:
                exact_counts[sublib] += 1
            else:
                mismatch_counts[sublib] += 1
    finally:
        for w_r1, w_r2 in writers.values():
            w_r1.stdin.close()
            w_r2.stdin.close()
        for w_r1, w_r2 in writers.values():
            w_r1.wait()
            w_r2.wait()
        r1.close()
        r2.close()
        r1_proc.wait()
        r2_proc.wait()

    per_sublibrary = {
        sub: {
            "exact_match": exact_counts.get(sub, 0),
            "1bp_mismatch": mismatch_counts.get(sub, 0),
            "total": exact_counts.get(sub, 0) + mismatch_counts.get(sub, 0),
        }
        for sub in i2_map
    }
    stats = {
        "per_sublibrary_read_pairs": per_sublibrary,
        "unmatched_read_pairs": unmatched,
        "total_input_read_pairs": sum(exact_counts.values()) + sum(mismatch_counts.values()) + unmatched,
        "i2_map": i2_map,
    }
    Path(args.stats).parent.mkdir(parents=True, exist_ok=True)
    Path(args.stats).write_text(json.dumps(stats, indent=2))
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

# Authors: Ben Parks
# Last updated: 9/26/22

# Run bcl2fastq for a SHARE-seq library
# Inputs:
#   - bcl run directory
#   - Simplified SampleSheet.csv file with only index 2 sequences
#     (see details for example formatting)
# Outputs:
#   - I1, I2, R1, R2 fastq.gz files
#     These are named like: SampleName1_R1.fastq.gz, SampleName2_I2.fastq.gz, etc.

# Details:
# - Because bcl2fastq doesn't like to demultiplex on only some of the index reads,
#   we use a bases mask to pretend that I1 is not an index read.
#   After that, we rename files so that we keep the normall I1/I2/R1/R2 notation,
#   and simplify the output naming scheme.
# - We also support exporting only a fraction of the total flow cell tiles, which
#   can be used for extracting a small test dataset, or for parallelizing bcl2fastq
#   across multiple nodes
#
# SampleSheet.csv example:
# Sample, I2
# ATAC_1, ATGCATGC
# ATAC_2, CATGCATG

# [Data]
# Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
# ATAC_1,,,,,,,ATGCATGC,,
# ATAC_2,,,,,,,CATGCATG,,

import argparse
from pathlib import Path
import subprocess
import sys
import xml.parsers.expat


def main():
    parser = argparse.ArgumentParser(description="Demultiplex SHARE-seq using bcl2fastq")
    # Required arguments:
    parser.add_argument("--samples", required=True, type=str, help="Path of samples csv (2-column csv with header of Sample, I2)")
    parser.add_argument("--input", required=True, type=str, help="Illumina run directory")
    parser.add_argument("--output", required=True, type=str, help="Path of output directory")
    parser.add_argument("--log", required=True, type=str, help="Path of bcl2fastq log")
    
    # Optional arguments:
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use")
    parser.add_argument("--num_tile_chunks", type=int, default=1, help="Number of chunks to split the tiles into")
    parser.add_argument("--tile_chunk", type=int, default=1, help="Which tile chunk to extract (1-based index)")
    parser.add_argument("--short_mask", type=int, default=0, help="See bcl2fastq documentation option --mask-short-adapter-reads")

    args = parser.parse_args()

    log = open(args.log, "wb")
    
    # Parse samples list
    samples = {}
    for i, line in enumerate(open(args.samples)):
        sample, I2 = line.strip().split(",")
        sample = sample.strip()
        I2 = I2.strip().upper()
        if i == 0:
            if sample.lower() != "sample" or I2 != "I2":
                print("Error: First line of --samples csv must be Sample, I2")
                sys.exit(1)
        else:
            samples[sample] = I2
    
    if len(set(samples.keys())) != len(samples.keys()):
        print("Error: --samples must have unique sample names")
        sys.exit(1)
    
    # Write samplesheet
    samplesheet_path = Path(args.output).parent / "SampleSheet.csv"
    with open(samplesheet_path, "w") as samplesheet:
        samplesheet.write("[Data]\n")
        samplesheet.write("Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description\n")
        for sample, I2 in samples.items():
            samplesheet.write(f"{sample},,,,bc,,,{I2},,\n")

    tiles = get_tiles(args.input)
    tiles = tiles[args.tile_chunk - 1::args.num_tile_chunks]

    # Run bcl2fastq
    command = [
        "bcl2fastq", 
        "--runfolder-dir", args.input,
        "--sample-sheet", str(samplesheet_path), 
        "--interop-dir", args.output + "/InterOp",
        "--no-lane-splitting",
        "--create-fastq-for-index-reads",
        "--use-bases-mask", get_bases_mask(args.input),
        "--mask-short-adapter-reads", str(args.short_mask),
        "--output", args.output,
        "--processing-threads", str(args.threads), # Empirically, vast majority of the work is spent in processing threads, so we don't add threads for other purposes
        "--tiles", ",".join(tiles),
    ]
    bcl2fastq = subprocess.run(command, stderr=log)
    log.close()

    if bcl2fastq.returncode != 0:
        sys.exit(bcl2fastq.returncode)

    # Rename output files to match the right file naming convention
    i = 0
    read_rename = {
        "R1": "R1", 
        "R2": "I1",
        "R3": "R2",
        "I1": "I2",
    }
    
    for sample in samples:
        for old_read, new_read in read_rename.items():
            bcl2fastq_path = Path(f"{args.output}/{sample}_S{i+1}_{old_read}_001.fastq.gz")
            bcl2fastq_path.rename(f"{args.output}/{sample}_{new_read}.fastq.gz")
        i += 1

def get_bases_mask(run_dir):
    """Parse out a bases mask from the run_dir, switching I1 to be treated like a read"""
    reads = []
    p = xml.parsers.expat.ParserCreate()
    
    def start_element(name, attrs):
        if name == "Read":
            reads.append(attrs)
    p.StartElementHandler = start_element
    p.ParseFile(open(Path(run_dir) / "RunInfo.xml", "rb"))
    
    assert len(reads) == 4
    reads.sort(key=lambda x: x["Number"])
    # Carefully check our parsing
    for i, read in enumerate(reads):
        assert int(read["Number"]) == i+1
        if i in [0,3]:
            assert read["IsIndexedRead"] == "N"
        else:
            assert read["IsIndexedRead"] == "Y"
    # Make the bases mask
    return (
        f"Y{reads[0]['NumCycles']},"
        f"Y{reads[1]['NumCycles']},"
        f"I{reads[2]['NumCycles']},"
        f"Y{reads[3]['NumCycles']}"
    )

def get_tiles(run_dir):
    """Parse out the list of tiles from the run_dir"""
    p = xml.parsers.expat.ParserCreate()
    
    tiles = []
    is_tile = False
    cur_tile = ""
    
    def start_element(name, attrs):
        nonlocal is_tile
        is_tile = name == "Tile"
    def end_element(name):
        nonlocal is_tile, cur_tile
        if is_tile:
            tiles.append(cur_tile)
            cur_tile = ""
        is_tile = False
    def char_data(data):
        nonlocal cur_tile
        if is_tile:
            cur_tile += data
    
    p.StartElementHandler = start_element
    p.EndElementHandler = end_element
    p.CharacterDataHandler = char_data
    p.ParseFile(open(Path(run_dir) / "RunInfo.xml", "rb"))
    
    return tiles

if __name__ == "__main__":
    main()

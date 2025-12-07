# Script to anonymize SHARE-seq FASTQ reads by removing identifying genetic information
# Author: Betty Liu
# Last Modified: 11/12/2025

# Primary outputs:
# - {assay}/sublibraries/{assay}_{sublibrary}_anon_R1/R2.fastq.gz: Anonymized fastq files per sublibrary
# - {assay}/samples/{assay}_{sample}_anon_{read}.fastq.gz: Anonymized fastq files per sample
# - {assay}/samples/raw/{assay}_{sample}_raw_{read}.fastq.gz: Raw fastq files per sample

import collections
import os
import re

import utils

workdir: config["output_dir"]

# global singularity container to use
# only set if container given in config and is not none
if "singularity" in config.keys() and config["singularity"]: 
    singularity: config["singularity"] 
    
#############################
### Config parsing and metadata helpers
#############################

chunk_size = config["chunk_size"]

# Fix up any config keys that are not strings
utils.string_only_keys(config)

def get_chunks(sequencing_path):
    """Generate chunk IDs for a sequencing path based on read count. Adds padding 0s as needed"""
    reads = int(open(f"{sequencing_path}/read_count.txt").read())
    chunk_count = (reads + chunk_size - 1) // chunk_size
    if "test_chunks" in config:
        chunk_count = min(chunk_count, config["test_chunks"])
    str_len = max(2, len(str(chunk_count)))
    return [f"{i:0{str_len}d}" for i in range(1, chunk_count+1)]

def expand_sublibrary_chunks(pattern, assay, w):
    """Generate all sequencing_path/chunks for a given sublibrary """
    results = []
    for seqpath in utils.get_sequencing_paths(assay, config, sublib=w.sublibrary):
        results += expand(pattern, sequencing_path=seqpath, chunk=get_chunks(seqpath))
    return results

def expand_assay_chunks(pattern, assay):
    """Generate all sequencing_path/chunks for a given assay """
    results = []
    for seqpath in utils.get_sequencing_paths(assay, config):
        results += expand(pattern, sequencing_path=seqpath, chunk=get_chunks(seqpath), allow_missing=True)
    return results

# Confirm that we have read counts for all the input sequences
for sequencing_path in utils.get_sequencing_paths("ATAC", config) + utils.get_sequencing_paths("RNA", config):
    if not os.path.exists(f"{sequencing_path}/read_count.txt"):
        raise RuntimeError(f"Must run prep_fastq.smk; missing read counts for: {sequencing_path}")
del sequencing_path

wildcard_constraints:
    chunk = "\d+", # Chunk is a number
    sequencing_path = "(ATAC|RNA)/([^/]+/)?[^/]+", # Sequencing path is 2-3 folders
    sample = "|".join(re.escape(s) for s in config["samples"].keys())

barcodes = utils.bc_names(srcdir("config/barcodes/Round1.tsv"))

sample_barcodes = {
    sample: [b for b in barcodes if utils.grep_regex_match(b, regex)] for sample, regex in config["samples"].items()
}

# Check that each barcode is used exactly once

used_barcodes = [b for l in sample_barcodes.values() for b in l]

if len(used_barcodes) != len(set(used_barcodes)) or set(used_barcodes) != set(barcodes):
    duplicates = [b for b in barcodes if used_barcodes.count(b) > 1]
    missing = [b for b in barcodes if used_barcodes.count(b) == 0]
    raise RuntimeError(f"Not all barcodes used exactly once! Duplicates: {duplicates}, Missing: {missing}")
del used_barcodes
del barcodes

outputs = []
if len(utils.get_sequencing_paths("ATAC", config)) > 0:
    outputs += (
        expand('ATAC/sublibraries/ATAC_{sublibrary}_anon_{read}.fastq.gz', sublibrary=utils.get_sublibraries("ATAC", config), read=["R1","R2"]) +
        expand('{assay}/samples/{assay}_{sample}_anon_{read}.fastq.gz', assay=["ATAC"], sample=config["samples"].keys(), read=["R1", "R2"]) +
        expand('{assay}/samples/raw/{assay}_{sample}_raw_{read}.fastq.gz', assay=["ATAC"], sample=config["samples"].keys(), read=["R1", "R2"])
    )

if len(utils.get_sequencing_paths("RNA", config)) > 0:
    outputs += (
        expand('RNA/sublibraries/RNA_{sublibrary}_anon_{read}.fastq.gz', sublibrary=utils.get_sublibraries("RNA", config), read=["R1", "R2"]) +
        expand('{assay}/samples/{assay}_{sample}_anon_{read}.fastq.gz', assay=["RNA"], sample=config["samples"].keys(), read=["R1", "R2"]) +
        expand('{assay}/samples/raw/{assay}_{sample}_raw_{read}.fastq.gz', assay=["RNA"], sample=config["samples"].keys(), read=["R1", "R2"])
    )

if "filter_dag" in config.keys() and config["filter_dag"]=="false":
    filtered_outputs = outputs
else:
    filtered_outputs = []
    for o in outputs:
        if os.path.exists(o):
            print(f"Skipping existing output: {o}", file=sys.stderr)
        else:
            filtered_outputs.append(o)

rule all:
    input: filtered_outputs 

#############################
### ATAC + RNA fastq processing 
#############################
# Split fastqs
rule split_fastqs:
    input:
        fastq = lambda w: utils.fastq_path(w.sequencing_path, w.read, config),
        read_count = "{sequencing_path}/read_count.txt",
    output:
        chunks = temp(directory("{sequencing_path}/split_fastqs/{read}"))
    params:
        decompress = lambda w: utils.fastq_decompress(w.sequencing_path, config),
        lines = chunk_size * 4,
        suffix_length = lambda w: len(get_chunks(w.sequencing_path)[0]),
        truncate_test_chunks = lambda w: f" | head -n {chunk_size*config['test_chunks']*4} " if "test_chunks" in config else ""
    resources:
        runtime = 60 * 5, # Be generous on time in case of large fastqs
    threads: 3
    log: '{sequencing_path}/split_fastqs/{read}.log'
    shell: "mkdir {output.chunks} && "
          " split <({params.decompress} {input.fastq} {params.truncate_test_chunks}) "
          " --numeric-suffixes=1 --lines {params.lines} "
          " --suffix-length={params.suffix_length} "
          " --additional-suffix=.fastq.zst "
          " --filter='zstd --fast=1 -q -o $FILE' "
          " {output.chunks}/ 2> {log}"


#############################
### ATAC-specific workflow 
#############################

# Remove adapter ends from the raw fastq reads.
# Discards reads that result in <15bp (--length_required=15 by default)
rule atac_trim_adapters:
    input: 
        R1 = expand(rules.split_fastqs.output.chunks, read="R1", allow_missing=True),
        R2 = expand(rules.split_fastqs.output.chunks, read="R2", allow_missing=True),
    output:
        de_R1 = temp("{sequencing_path}/{chunk}/R1.fastq.gz"), 
        de_R2 = temp("{sequencing_path}/{chunk}/R2.fastq.gz"), 
        #interleaved = temp("{sequencing_path}/{chunk}/anon_01_trim_adapters.interleaved.fastq.zst"),
        interleaved = temp("{sequencing_path}/{chunk}/anon_01_trim_adapters.interleaved.fastq"),
        report_json = "{sequencing_path}/{chunk}/qc_stats/anon_01_trim_adapters.json",
        report_html = "{sequencing_path}/{chunk}/qc_stats/anon_01_trim_adapters.html",
    params:
        R1_in = "{sequencing_path}/split_fastqs/R1/{chunk}.fastq.zst",
        R2_in = "{sequencing_path}/split_fastqs/R2/{chunk}.fastq.zst",
    threads: 4
    log: '{sequencing_path}/{chunk}/anon_01_trim_adapters.log'
    shell: 
        "zstd -dc {params.R1_in} | gzip -c > {output.de_R1}; zstd -dc {params.R2_in} | gzip -c > {output.de_R2}; "
        "fastp --in1 {output.de_R1} --in2 {output.de_R2} "
        " --adapter_sequence    CTGTCTCTTATACACATCTCCGAGCCCACGAGAC "
        " --adapter_sequence_r2 CTGTCTCTTATACACATCTGACGCTGCCGACGA "
        " -j {output.report_json} -h {output.report_html} "
        " -G -Q -w {threads} 2> {log} "
        #" --stdout | zstd --fast=1 -q -o {output.interleaved}"
        " --stdout > {output.interleaved}"


# Align ATAC reads with bowtie2, and no filtering
rule atac_bowtie2:
    input: 
        fastq = rules.atac_trim_adapters.output.interleaved
    output: 
        bam = temp('{sequencing_path}/{chunk}/anon_02_atac_bowtie2.bam'),
        index = temp('{sequencing_path}/{chunk}/anon_02_atac_bowtie2.bam.bai')
    params:
        index = config["genome"]["bowtie2"]
    resources:
        runtime = min(60, 5 * config["chunk_size"] // 1_000_000), # 5 minutes-per 1M read time estimate
        mem_mb = 64000
    threads: 8
    log: '{sequencing_path}/{chunk}/anon_02_atac_bowtie2.log',
    shell: 
           #"bowtie2 --interleaved <(zstd -dc {input.fastq}) -x {params.index} " # this sometimes leads to truncated files due to piping glitches
           "bowtie2 --interleaved {input.fastq} -x {params.index} "
           " --sam-append-comment --maxins 2000 --threads {threads} 2> {log} | "
           " samtools sort -@ {threads} > {output.bam} ; "
           " samtools index {output.bam}"

# Anonymize bams
rule atac_anon_bam:
    input:
        bam = rules.atac_bowtie2.output.bam,
        index = rules.atac_bowtie2.output.index
    output:
        bam = temp('{sequencing_path}/{chunk}/anon_03_atac_bowtie2_anon.bam'),
        index = temp('{sequencing_path}/{chunk}/anon_03_atac_bowtie2_anon.bam.bai')
    params:
        fasta = config["genome"]["fasta"]
    resources:
        mem_mb = 64000
    threads: 16
    shell: "BAMboozle --bam {input.bam} --out {output.bam} --fa {params.fasta} "
           " --p {threads}"

# Convert anonymized bams to fastqs
rule atac_convert_bamtofastq:
    input:
        bam = rules.atac_anon_bam.output.bam
    output:
        R1 = temp('{sequencing_path}/{chunk}/anon_04_atac_anon_R1.fastq.gz'),
        R2 = temp('{sequencing_path}/{chunk}/anon_04_atac_anon_R2.fastq.gz'),
    shell:
        "samtools sort -n {input.bam} | "
        " samtools fastq -T BC "
        " -1 >(sed 's/BC:Z://g' | gzip -c > {output.R1}) "
        " -2 >(sed 's/BC:Z://g' | gzip -c > {output.R2})"

# re-pair read 1 and read 2
rule atac_match_r2:
    input:
        R1 = rules.atac_convert_bamtofastq.output.R1,
        R2 = rules.atac_convert_bamtofastq.output.R2,
    output:
        R1 = temp('{sequencing_path}/{chunk}/anon_04_atac_anon_R1.paired.fastq.gz'),
        R2 = temp('{sequencing_path}/{chunk}/anon_04_atac_anon_R2.paired.fastq.gz')
    resources:
        mem_mb = 64000
    shell:
        "seqkit pair -1 {input.R1} -2 {input.R2}"

rule atac_merge_chunks_fastq:
    input:
        R1s = lambda w: expand_sublibrary_chunks(rules.atac_match_r2.output.R1, "ATAC", w),
        R2s = lambda w: expand_sublibrary_chunks(rules.atac_match_r2.output.R2, "ATAC", w)
    output:
        R1 = 'ATAC/sublibraries/ATAC_{sublibrary}_anon_R1.fastq.gz',
        R2 = 'ATAC/sublibraries/ATAC_{sublibrary}_anon_R2.fastq.gz'
    shell: "cat {input.R1s} > {output.R1} && "
           "cat {input.R2s} > {output.R2}  "

# Perform barcode matching (round 1 only)
rule atac_match_barcodes:
    input: 
        R1 = rules.atac_match_r2.output.R1,
        R2 = rules.atac_match_r2.output.R2
    output:
        R1 = temp("{sequencing_path}/{chunk}/anon_ATAC_match_barcodes_R1.fastq.zst"),
        R2 = temp("{sequencing_path}/{chunk}/anon_ATAC_match_barcodes_R2.fastq.zst"),
        stats = "{sequencing_path}/{chunk}/qc_stats/anon_match_barcodes.json",
    params:
        script = srcdir("scripts/shareseq/match_barcodes_r1_only.py"),
        BC1 = srcdir("config/barcodes/Round1.tsv"),
    threads: 2
    log: "{sequencing_path}/{chunk}/anon_match_barcodes.log"
    shell: "python3 {params.script} "
        " --R1_in <(gzip -dc {input.R1}) --R2_in <(gzip -dc {input.R2}) "
        " --R1_out {output.R1} "
        " --R2_out {output.R2} "
        " --output-cmd 'zstd --fast=1 -q -o $FILE' "
        " --BC1 {params.BC1} "
        " --json_stats {output.stats} "
        " 2> {log} "

rule atac_match_barcodes_raw:
    input: 
        R1 = expand(rules.split_fastqs.output.chunks, read="R1", allow_missing=True),
        R2 = expand(rules.split_fastqs.output.chunks, read="R2", allow_missing=True),
    output:
        R1 = temp("{sequencing_path}/{chunk}/raw_ATAC_match_barcodes_R1.fastq.zst"),
        R2 = temp("{sequencing_path}/{chunk}/raw_ATAC_match_barcodes_R2.fastq.zst"),
        stats = "{sequencing_path}/{chunk}/qc_stats/raw_match_barcodes.json",
    params:
        script = srcdir("scripts/shareseq/match_barcodes_r1_only.py"),
        BC1 = srcdir("config/barcodes/Round1.tsv"),
        R1_in = "{sequencing_path}/split_fastqs/R1/{chunk}.fastq.zst",
        R2_in = "{sequencing_path}/split_fastqs/R2/{chunk}.fastq.zst",
    threads: 2
    log: "{sequencing_path}/{chunk}/raw_match_barcodes.log"
    shell: "python3 {params.script} "
        " --R1_in <(zstd -dc {params.R1_in}) --R2_in <(zstd -dc {params.R2_in}) "
        " --R1_out {output.R1} "
        " --R2_out {output.R2} "
        " --output-cmd 'zstd --fast=1 -q -o $FILE' "
        " --BC1 {params.BC1} "
        " --json_stats {output.stats} "
        " 2> {log} "

#############################
### RNA-specific workflow 
#############################

# Remove adapter ends from the raw R1 fastq reads
rule rna_trim_adapters:
    input: 
        R1 = expand(rules.split_fastqs.output.chunks, read="R1", allow_missing=True),
    output:
        R1 = temp("{sequencing_path}/{chunk}/anon_01_trim_adapters.R1.fastq.zst"),
        report_json = "{sequencing_path}/{chunk}/qc_stats/anon_01_trim_adapters.json",
        report_html = "{sequencing_path}/{chunk}/qc_stats/anon_01_trim_adapters.html",
    params:
        R1_in = "{sequencing_path}/split_fastqs/R1/{chunk}.fastq.zst"
    threads: 4
    log: '{sequencing_path}/{chunk}/anon_01_trim_adapters.log'
    shell: "fastp --in1 <(zstd -dc {params.R1_in}) "
        " --adapter_sequence    CTGTCTCTTATACACATCTCCGAGCCCACGAGAC "
        " -j {output.report_json} -h {output.report_html} "
        " -G -Q -L -w {threads} 2> {log} "
        " --stdout | zstd --fast=1 -q -o {output.R1}"

# Change readname format so the indices are connected with coordinates, otherwise STAR will omit anything after white space
rule rna_format_BC:
    input: 
        R1 = rules.rna_trim_adapters.output.R1
    output:
        R1 = temp("{sequencing_path}/{chunk}/anon_02_trim_adapters_formatted.R1.fastq.zst"),
    shell:
        "zstd -dc {input.R1} | "
        " awk 'NR % 4 == 1 && substr($0,1,1) == \"@\" {{sub(/ /, \"__\", $0)}}; {{print}}' | "
        " zstd -q -o {output.R1}"

# Align RNA reads with star, no filtering
rule rna_star:
    input: 
        fastq = rules.rna_format_BC.output.R1
    output: 
        bam = temp('{sequencing_path}/{chunk}/anon_03_rna_star_Aligned.sortedByCoord.out.bam'),
        index = temp('{sequencing_path}/{chunk}/anon_03_rna_star_Aligned.sortedByCoord.out.bam.bai'),
        sj = temp('{sequencing_path}/{chunk}/anon_03_rna_star_SJ.out.tab'),
        log_prog = temp('{sequencing_path}/{chunk}/anon_03_rna_star_Log.progress.out'),
    params:
        index = config["genome"]["star"],
        prefix = "{sequencing_path}/{chunk}/anon_03_rna_star_"
    resources:
        runtime = min(60, 5 * config["chunk_size"] // 1_000_000), # 5 minutes-per 1M read time estimate
        mem_mb = 64000,
    threads: 16
    log: 
        setup = '{sequencing_path}/{chunk}/anon_03_rna_star_Log.out',
        summary = '{sequencing_path}/{chunk}/anon_03_rna_star_Log.final.out'
    shell: " STAR --chimOutType WithinBAM "
           " --runThreadN {threads} "
           " --genomeDir {params.index} "
           " --readFilesIn {input.fastq} "
           " --readFilesCommand zstd -dc "   
           " --outSAMattributes NH HI AS NM MD"    
           " --outSAMtype BAM SortedByCoordinate " 
           " --outSAMunmapped Within "
           " --outSAMstrandField intronMotif "
           " --outReadsUnmapped None "
           " --outFileNamePrefix {params.prefix} "     
           " --outFilterType Normal "
           " --outFilterMultimapNmax 999999 "
           " --outFilterMismatchNmax 999 "
           " --outFilterMismatchNoverReadLmax 1.0 "
           " --outFilterScoreMinOverLread 0 "
           " --outFilterMatchNminOverLread 0 "
           " --alignIntronMin 1 "
           " --alignIntronMax 1000000 "
           " --alignMatesGapMax 1000000 "
           " --alignSJoverhangMin 1 "
           " --alignSJDBoverhangMin 1 "
           " --sjdbScore 1 "
           " --limitOutSJcollapsed 5000000 ;"
           " samtools index {output.bam}"

# Anonymize bams
rule rna_anon_bam:
    input:
        bam = rules.rna_star.output.bam
    output:
        bam = temp('{sequencing_path}/{chunk}/anon_04_rna_star_Aligned_anon.out.bam'),
    params:
        fasta = config["genome"]["fasta"]
    resources:
        mem_mb = 16000
    threads: 8
    shell: "BAMboozle --bam {input.bam} --out {output.bam} --fa {params.fasta} "
           " --p {threads}"

# Convert anonymized bams to fastqs
rule rna_convert_bamtofastq:
    input:
        bam = rules.rna_anon_bam.output.bam
    output:
        R1 = temp('{sequencing_path}/{chunk}/anon_05_rna_anon_R1.fastq.gz'),
    shell:
        "samtools sort -n {input.bam} | "
        " samtools fastq -0 >(sed 's/__/\\t/g' | gzip -c > {output.R1})"

# Anonymize R2 for RNA, only keep the first 10bp of read 2 (UMI)
#   the rest of read 2 is difficult to align and anonymize due to polyA so we truncate to first 10bp only
rule rna_match_r2:
    input:
        R1 = rules.rna_convert_bamtofastq.output.R1,
        R2 = expand(rules.split_fastqs.output.chunks, read="R2", allow_missing=True)
    output:
        de_R2 = temp("{sequencing_path}/{chunk}/anon_06_rna_anon_R2_untrimmed.fastq.gz"),
        R1 = temp('{sequencing_path}/{chunk}/anon_05_rna_anon_R1.paired.fastq.gz'),
        R2_tmp = temp('{sequencing_path}/{chunk}/anon_06_rna_anon_R2_untrimmed.paired.fastq.gz'),
        R2 = temp('{sequencing_path}/{chunk}/anon_06_rna_anon_R2.paired.fastq.gz')
    params:
        R2_in = "{sequencing_path}/split_fastqs/R2/{chunk}.fastq.zst",
    resources:
        mem_mb = 64000
    shell:
        "zstd -dc {params.R2_in} | seqkit sort -n | gzip -c > {output.de_R2};"
        "seqkit pair -1 {input.R1} -2 {output.de_R2};"
        "cat {output.R2_tmp} | seqkit subseq -r 1:10 | gzip -c > {output.R2}"
    
rule rna_merge_chunks_fastq:
    input:
        R1s = lambda w: expand_sublibrary_chunks(rules.rna_match_r2.output.R1, "RNA", w),
        R2s = lambda w: expand_sublibrary_chunks(rules.rna_match_r2.output.R2, "RNA", w),
    output:
        R1 = 'RNA/sublibraries/RNA_{sublibrary}_anon_R1.fastq.gz',
        R2 = 'RNA/sublibraries/RNA_{sublibrary}_anon_R2.fastq.gz',
    shell: "cat {input.R1s} > {output.R1} && "
           "cat {input.R2s} > {output.R2} "

# Perform barcode matching (round 1 only)
rule rna_match_barcodes:
    input: 
        R1 = rules.rna_match_r2.output.R1,
        R2 = rules.rna_match_r2.output.R2
    output:
        R1 = temp("{sequencing_path}/{chunk}/anon_RNA_match_barcodes_R1.fastq.zst"),
        R2 = temp("{sequencing_path}/{chunk}/anon_RNA_match_barcodes_R2.fastq.zst"),
        stats = "{sequencing_path}/{chunk}/qc_stats/anon_match_barcodes.json"
    params:
        script = srcdir("scripts/shareseq/match_barcodes_r1_only.py"),
        BC1 = srcdir("config/barcodes/Round1.tsv"),
    threads: 2
    log: "{sequencing_path}/{chunk}/anon_match_barcodes.log"
    shell: "python3 {params.script} "
        " --R1_in <(gzip -dc {input.R1}) --R2_in <(gzip -dc {input.R2}) "
        " --R1_out {output.R1} "
        " --R2_out {output.R2} "
        " --output-cmd 'zstd --fast=1 -q -o $FILE' "
        " --BC1 {params.BC1} "
        " --json_stats {output.stats} "
        " 2> {log} "

rule rna_match_barcodes_raw:
    input: 
        R1 = expand(rules.split_fastqs.output.chunks, read="R1", allow_missing=True),
        R2 = expand(rules.split_fastqs.output.chunks, read="R2", allow_missing=True),
    output:
        R1 = temp("{sequencing_path}/{chunk}/raw_RNA_match_barcodes_R1.fastq.zst"),
        R2 = temp("{sequencing_path}/{chunk}/raw_RNA_match_barcodes_R2.fastq.zst"),
        stats = "{sequencing_path}/{chunk}/qc_stats/raw_match_barcodes.json"
    params:
        script = srcdir("scripts/shareseq/match_barcodes_r1_only.py"),
        BC1 = srcdir("config/barcodes/Round1.tsv"),
        R1_in = "{sequencing_path}/split_fastqs/R1/{chunk}.fastq.zst",
        R2_in = "{sequencing_path}/split_fastqs/R2/{chunk}.fastq.zst",
    threads: 2
    log: "{sequencing_path}/{chunk}/raw_match_barcodes.log"
    shell: "python3 {params.script} "
        " --R1_in <(zstd -dc {params.R1_in}) --R2_in <(zstd -dc {params.R2_in}) "
        " --R1_out {output.R1} "
        " --R2_out {output.R2} "
        " --output-cmd 'zstd --fast=1 -q -o $FILE' "
        " --BC1 {params.BC1} "
        " --json_stats {output.stats} "
        " 2> {log} "

#############################
### Joint processing
#############################
# localrules: split_samples # if hitting slurm job submission limits, use localrules
rule split_samples:
    input: 
        R1 = expand("{sequencing_path}/{chunk}/anon_{assay}_match_barcodes_R1.fastq.zst", allow_missing=True),
        R2 = expand("{sequencing_path}/{chunk}/anon_{assay}_match_barcodes_R2.fastq.zst", allow_missing=True)
    output:
        R1 = temp('{sequencing_path}/{chunk}/split_samples/anon_{assay}_{sample}_R1.fastq.gz'),
        R2 = temp('{sequencing_path}/{chunk}/split_samples/anon_{assay}_{sample}_R2.fastq.gz')
    params:
        barcode_pattern = lambda w: f"_CB:Z:({config['samples'][w.sample]})",
        keep_ids = lambda w: f"{w.sequencing_path}/{w.chunk}/split_samples/{w.sample}_keep_ids.txt"
    resources:
        runtime= 60 * 2
    shell: "zstd -dc {input.R1} | "
        "seqkit grep -r -n -p '{params.barcode_pattern}' | "
        "seqkit seq -n | awk '{{print $1}}' > {params.keep_ids} && "
        "zstd -dc {input.R1} | seqkit grep -f {params.keep_ids} | seqkit replace -p '_CB:Z:.*$' | gzip > {output.R1} && "
        "zstd -dc {input.R2} | seqkit grep -f {params.keep_ids} | seqkit replace -p '_CB:Z:.*$' | gzip > {output.R2} && "
        "rm -f {params.keep_ids}"

rule merge_samples:
    input: 
        R1s = lambda w: expand_assay_chunks(rules.split_samples.output.R1, w.assay),
        R2s = lambda w: expand_assay_chunks(rules.split_samples.output.R2, w.assay)
    output:
        R1 = '{assay}/samples/{assay}_{sample}_anon_R1.fastq.gz',
        R2 = '{assay}/samples/{assay}_{sample}_anon_R2.fastq.gz',
    threads: 8
    resources:
        runtime= 60 * 2
    shell:"cat {input.R1s} > {output.R1} && "
          "cat {input.R2s} > {output.R2} "

# localrules: split_samples_raw # if hitting slurm job submission limits, use localrules
rule split_samples_raw:
    input: 
        R1 = expand("{sequencing_path}/{chunk}/raw_{assay}_match_barcodes_R1.fastq.zst", allow_missing=True),
        R2 = expand("{sequencing_path}/{chunk}/raw_{assay}_match_barcodes_R2.fastq.zst", allow_missing=True)
    output:
        R1 = temp('{sequencing_path}/{chunk}/split_samples/raw_{assay}_{sample}_R1.fastq.gz'),
        R2 = temp('{sequencing_path}/{chunk}/split_samples/raw_{assay}_{sample}_R2.fastq.gz')
    params:
        barcode_pattern = lambda w: f"_CB:Z:({config['samples'][w.sample]})",
        keep_ids = lambda w: f"{w.sequencing_path}/{w.chunk}/split_samples/{w.sample}_keep_ids.txt"
    resources:
        runtime= 60 * 2
    shell: "zstd -dc {input.R1} | "
        "seqkit grep -r -n -p '{params.barcode_pattern}' | "
        "seqkit seq -n | awk '{{print $1}}' > {params.keep_ids} && "
        "zstd -dc {input.R1} | seqkit grep -f {params.keep_ids} | seqkit replace -p '_CB:Z:.*$' | gzip > {output.R1} && "
        "zstd -dc {input.R2} | seqkit grep -f {params.keep_ids} | seqkit replace -p '_CB:Z:.*$' | gzip > {output.R2} && "
        "rm -f {params.keep_ids}"

rule merge_samples_raw:
    input: 
        R1s = lambda w: expand_assay_chunks(rules.split_samples_raw.output.R1, w.assay),
        R2s = lambda w: expand_assay_chunks(rules.split_samples_raw.output.R2, w.assay)
    output:
        R1 = '{assay}/samples/raw/{assay}_{sample}_raw_R1.fastq.gz',
        R2 = '{assay}/samples/raw/{assay}_{sample}_raw_R2.fastq.gz',
    threads: 8
    resources:
        runtime= 60 * 2
    shell:"cat {input.R1s} > {output.R1} && "
          "cat {input.R2s} > {output.R2} "
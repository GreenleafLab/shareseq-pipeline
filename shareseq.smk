# Main shareseq pipeliine
# Author: Ben Parks, Betty Liu
# Last Modified: 03/29/2023

# Primary outputs:
# - ATAC/samples/{sample}.fragments.tsv.gz: 10x-compatible fragment file (though end coordinates are +1bp relative to 10x)
# - RNA/samples/{sample}.matrix.mtx.gz: 10x-compatible MatrixMarket files

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

# Confirm that we have read counts for all the input sequences
for sequencing_path in utils.get_sequencing_paths("ATAC", config) + utils.get_sequencing_paths("RNA", config):
    if not os.path.exists(f"{sequencing_path}/read_count.txt"):
        raise RuntimeError(f"Must run prep_fastq.smk; missing read counts for: {sequencing_path}")
del sequencing_path

wildcard_constraints:
    chunk = "\d+", # Chunk is a number
    barcode_chunk = "barcodes_\d+", # Barcode chunk is barcode_ folowed by a number
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

# Make groupings of last-round barcodes for use in RNA-seq umi deduplication
end_barcodes = utils.bc_names(srcdir("config/barcodes/Round3.tsv"))
barcode_group_size = 8
end_barcode_groups = [
    end_barcodes[i:i+barcode_group_size] for i in range(0, len(end_barcodes), barcode_group_size)
]
barcode_chunks = [
    f"barcodes_{i:02d}" for i in range(1, len(end_barcode_groups) + 1)
]

outputs = []
if len(utils.get_sequencing_paths("ATAC", config)) > 0:
    outputs += (
        ["ATAC/samples/alignment_stats.json", "ATAC/samples/barcode_stats.json"] +
        expand('ATAC/samples/{sample}.fragments.tsv.gz', sample=config["samples"].keys()) +
        expand('ATAC/sublibraries/{sublibrary}/fragments.tsv.gz', sublibrary=utils.get_sublibraries("ATAC", config))
    )

if len(utils.get_sequencing_paths("RNA", config)) > 0:
    outputs += (
        ["RNA/samples/alignment_stats.json", "RNA/samples/barcode_stats.json"] +
        expand('RNA/samples/{sample}.{file}.gz', sample=config["samples"].keys(), file=["matrix.mtx", "barcodes.tsv", "features.tsv"]) +
        expand('RNA/sublibraries/{sublibrary}/{file}.gz', sublibrary=utils.get_sublibraries("RNA", config), file=["matrix.mtx", "barcodes.tsv", "features.tsv"])
    )

rule all:
    input: outputs 

#############################
### Build C-based scripts
#############################
localrules: build_count_unique
rule build_count_unique:
    input: srcdir("scripts/shareseq/count_unique.c")
    output: "bin/count_unique"
    shell: "gcc -O3 -o {output} {input}"

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

# Perform barcode matching
rule match_barcodes:
    input: 
        R1 = expand(rules.split_fastqs.output.chunks, read="R1", allow_missing=True),
        R2 = expand(rules.split_fastqs.output.chunks, read="R2", allow_missing=True)
    output:
        R1 = temp("{sequencing_path}/{chunk}/01_match_barcodes_R1.fastq.zst"),
        R2 = temp("{sequencing_path}/{chunk}/01_match_barcodes_R2.fastq.zst"),
        stats = "{sequencing_path}/{chunk}/qc_stats/01_match_barcodes.json",
    params:
        script = srcdir("scripts/shareseq/match_barcodes.py"),
        R1_in = "{sequencing_path}/split_fastqs/R1/{chunk}.fastq.zst",
        R2_in = "{sequencing_path}/split_fastqs/R2/{chunk}.fastq.zst",
        BC1 = srcdir("config/barcodes/Round1.tsv"),
        BC2 = srcdir("config/barcodes/Round2.tsv"),
        BC3 = srcdir("config/barcodes/Round3.tsv"),
        assay = lambda w: (w.sequencing_path).split("/")[0]
    threads: 2
    log: '{sequencing_path}/{chunk}/01_match_barcodes.log'
    shell: "python3 {params.script} "
        " --R1_in <(zstd -dc {params.R1_in}) --R2_in <(zstd -dc {params.R2_in}) "
        " --R1_out {output.R1} "
        " --R2_out {output.R2} "
        " --output-cmd 'zstd --fast=1 -q -o $FILE' "
        " --BC1 {params.BC1} --BC2 {params.BC2} --BC3 {params.BC3} "
        " --json_stats {output.stats} "
        " --assay {params.assay} "
        " 2> {log} "

#############################
### ATAC-specific workflow 
#############################

# Remove adapter ends from the raw fastq reads.
# Discards reads that result in <15bp (--length_required=15 by default)
rule atac_trim_adapters:
    input: 
        R1 = rules.match_barcodes.output.R1,
        R2 = rules.match_barcodes.output.R2,
    output:
        interleaved = temp("{sequencing_path}/{chunk}/02_trim_adapters.interleaved.fastq.zst"),
        report_json = "{sequencing_path}/{chunk}/qc_stats/02_trim_adapters.json",
        report_html = "{sequencing_path}/{chunk}/qc_stats/02_trim_adapters.html",
    threads: 4
    log: '{sequencing_path}/{chunk}/02_trim_adapters.log'
    shell: "fastp --in1 <(zstd -dc {input.R1}) --in2 <(zstd -dc {input.R2}) "
        " --adapter_sequence    CTGTCTCTTATACACATCTCCGAGCCCACGAGAC "
        " --adapter_sequence_r2 CTGTCTCTTATACACATCTGACGCTGCCGACGA "
        " -j {output.report_json} -h {output.report_html} "
        " -G -Q -w {threads} 2> {log} "
        " --stdout | zstd --fast=1 -q -o {output.interleaved}"

# Alternative trim command
# "SeqPurge -a1 CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a2 CTGTCTCTTATACACATCTGACGCTGCCGACGA "
#     " -qcut 0 -ncut 0 "
#     " -threads {threads} -out1 {output.R1} -out2 {output.R2} "
#     " -in1 {input.R1} -in2 {input.R2} > {log}"

# fastp args from ENCODE sc_atac pipeline
# "fastp -i {input.fastq1_bc} -I {input.fastq2_bc} -o {output.fastq1_trim} -O {output.fastq2_trim}"
#        " -h {log.html} -j {log.json} -G -Q -L -w {threads} 2> {output.stats}"

# Align ATAC reads with bowtie2, and filter to good quality reads only
rule atac_bowtie2:
    input: 
        fastq = rules.atac_trim_adapters.output.interleaved
    output: 
        bam = temp('{sequencing_path}/{chunk}/03_atac_bowtie2.bam'),
    params:
        index = config["genome"]["bowtie2"]
    resources:
        runtime = min(60, 2 * config["chunk_size"] // 1_000_000) # 2 minutes-per 1M read time estimate
    threads: 16
    log: '{sequencing_path}/{chunk}/03_atac_bowtie2.log',
    shell: "bowtie2 --interleaved <(zstd -dc {input.fastq}) -x {params.index} "
           " --sam-append-comment --maxins 2000 --threads {threads} 2> {log} | "
           # -F 1804: exclude flag, exludes unmapped, next segment unmapped, secondary alignments, not passing platform q, PCR or optical duplicates
           # -f 2: flags to require, properly aligned
           # -q 30: exlude low MAPQ, set as adjustable configuration parameter
           "samtools view -F 1804 -f 2 -q 30 -1 - > {output} "
        
# Convert bam to fragments format and sort for first pass
rule atac_convert_fragments:
    input:
        bam = rules.atac_bowtie2.output.bam
    output: 
        fragments = temp('{sequencing_path}/{chunk}/04_atac_convert_fragments.fragments.tsv.zst'),
    params:
        script = srcdir("scripts/shareseq/bam_to_fragments.py"),
        memory = "4G",
    threads: 4
    log: '{sequencing_path}/{chunk}/04_atac_convert_fragments.log',
    shell: "python {params.script} {input} 2> {log} | "
           "LC_ALL=C sort -k1,1V -k2,2n -k3,3n -k4,4 -t$'\\t' "
           "-S {params.memory} --parallel={threads} | "
           "zstd --fast=1 -q -o {output.fragments} "

rule atac_merge_chunks:
    input: 
        fragments = lambda w: expand_sublibrary_chunks(rules.atac_convert_fragments.output, "ATAC", w), 
        script = rules.build_count_unique.output
    output: 
        fragments = temp('ATAC/sublibraries/{sublibrary}/fragments.tsv.zst'),
    params:
        memory = "4G",
        fragments_decompress = lambda w: expand_sublibrary_chunks(f"<(zstd  -dc {rules.atac_convert_fragments.output})", "ATAC", w)
    threads: 5
    resources:
        runtime = 60 * 2
    shell: "LC_ALL=C sort -k1,1V -k2,2n -k3,3n -k4,4 -t$'\\t' --parallel=4 "
        " --merge --batch-size=200 -S {params.memory} {params.fragments_decompress} | "
        " {input.script} | "
        " zstd --fast=1 -q -o {output.fragments} "

rule atac_export_sublibrary:
    input:
        fragments = rules.atac_merge_chunks.output.fragments
    output:
        compressed = 'ATAC/sublibraries/{sublibrary}/fragments.tsv.gz',
        indexed = 'ATAC/sublibraries/{sublibrary}/fragments.tsv.gz.tbi',
    threads: 4
    resources: 
        runtime = 60 * 2
    shell: " zstd -dc {input.fragments} | "
           " bgzip -@ 4 -c > {output.compressed} && "
           " tabix --preset bed {output.compressed} "

rule atac_split_samples:
    input: 
        fragments = rules.atac_merge_chunks.output.fragments
    output:
        fragments = temp('ATAC/sublibraries/{sublibrary}/{sample}.tsv.zst'),
    params:
        barcode_pattern = lambda w: f"\t{config['samples'][w.sample]}\\+",
        sublibrary_id = lambda w: w.sublibrary
    threads: 4
    shell: "zstd -dc {input.fragments} | "
           "grep -E '{params.barcode_pattern}' | "
           "awk -c 'BEGIN {{OFS=\"\t\"}} {{$4=(\"{params.sublibrary_id}_\" $4); print $0}}' | " # Prefix sublibrary ID
           "zstd --fast=1 -q -o {output.fragments} "


rule atac_merge_samples:
    input: 
        fragments = lambda w: expand(rules.atac_split_samples.output.fragments, sublibrary=utils.get_sublibraries("ATAC", config), allow_missing=True),
    output:
        compressed = 'ATAC/samples/{sample}.fragments.tsv.gz',
        indexed = 'ATAC/samples/{sample}.fragments.tsv.gz.tbi',
    params:
        memory = "4G",
        fragments_decompress = lambda w: expand(f"<(zstd  -dc {rules.atac_split_samples.output.fragments})", sublibrary=utils.get_sublibraries("ATAC", config), sample=w.sample),
    threads: 8
    resources:
        runtime= 60 * 2
    shell:"LC_ALL=C sort -k1,1V -k2,2n -k3,3n -k4,4 -t$'\\t' "
        " --merge --batch-size=100 -S {params.memory} {params.fragments_decompress} --parallel=4 | "
        " bgzip -@ 4 -c > {output.compressed} && "
        " tabix --preset bed {output.compressed} "

rule atac_stats_libraries:
    input: 
        fragments = rules.atac_merge_chunks.output.fragments,
        bowtie2_log = lambda w: expand_sublibrary_chunks(rules.atac_bowtie2.log, "ATAC", w),
        barcode_stats = lambda w: expand_sublibrary_chunks(rules.match_barcodes.output.stats, "ATAC", w)        
    output: 
        summary = "ATAC/sublibraries/{sublibrary}/alignment_stats.json",
        barcodes = "ATAC/sublibraries/{sublibrary}/barcode_stats.json",
    params:
        script = srcdir("scripts/shareseq/stats_collect_atac.py")
    wildcard_constraints:
        sequencing_path = "ATAC/.*"
    shell: "python {params.script} "
          " --barcode_stats {input.barcode_stats} "
          " --bowtie2_log {input.bowtie2_log} "
          " --fragment_file <(zstd -dc {input.fragments}) "
          " --summary_output {output.summary} "
          " --barcodes_output {output.barcodes} "

localrules: atac_stats_merge
rule atac_stats_merge:
    input: 
        sequencing = expand(rules.atac_stats_libraries.output.summary, sublibrary=utils.get_sublibraries("ATAC", config)),
        barcodes = expand(rules.atac_stats_libraries.output.barcodes, sublibrary=utils.get_sublibraries("ATAC", config)),
    output: 
        sequencing = "ATAC/samples/alignment_stats.json",
        barcodes = "ATAC/samples/barcode_stats.json"
    params: 
        script = srcdir("scripts/shareseq/stats_aggregate.py")
    shell: "python {params.script} --input {input.sequencing} --output {output.sequencing};"
           "python {params.script} --input {input.barcodes} --output {output.barcodes};"

#############################
### RNA-specific workflow 
#############################

# Remove adapter ends from the raw R1 fastq reads
rule rna_trim_adapters:
    input: 
        R1 = rules.match_barcodes.output.R1
    output:
        R1 = temp("{sequencing_path}/{chunk}/02_trim_adapters.R1.fastq.zst"),
        report_json = "{sequencing_path}/{chunk}/qc_stats/02_trim_adapters.json",
        report_html = "{sequencing_path}/{chunk}/qc_stats/02_trim_adapters.html",
    threads: 4
    log: '{sequencing_path}/{chunk}/02_trim_adapters.log'
    shell: "fastp --in1 <(zstd -dc {input.R1}) "
        " --adapter_sequence    CTGTCTCTTATACACATCTCCGAGCCCACGAGAC "
        " -j {output.report_json} -h {output.report_html} "
        " -G -Q -L -w {threads} 2> {log} "
        " --stdout | zstd --fast=1 -q -o {output.R1}"

# Align RNA reads with star
rule rna_star:
    input: 
        fastq = rules.rna_trim_adapters.output.R1
    output: 
        bam = temp('{sequencing_path}/{chunk}/03_rna_star_Aligned.out.bam'),
        sj = temp('{sequencing_path}/{chunk}/03_rna_star_SJ.out.tab'),
        log_prog = temp('{sequencing_path}/{chunk}/03_rna_star_Log.progress.out'),
    params:
        index = config["genome"]["star"],
        prefix = "{sequencing_path}/{chunk}/03_rna_star_"
    resources:
        runtime = min(60, 2 * config["chunk_size"] // 1_000_000), # 2 minutes-per 1M read time estimate
        mem_mb = 64000,
    threads: 4
    log: 
        setup = '{sequencing_path}/{chunk}/03_rna_star_Log.out',
        summary = '{sequencing_path}/{chunk}/03_rna_star_Log.final.out'
    shell: " STAR --chimOutType WithinBAM "
           " --runThreadN {threads} "
           " --genomeDir {params.index} "
           " --readFilesIn {input.fastq} "
           " --readFilesCommand zstd -dc "   
           " --outFilterMultimapNmax 50 "
           " --outFilterScoreMinOverLread 0.3 "
           " --outFilterMatchNminOverLread 0.3 "
           " --outSAMattributes NH HI AS NM MD "    
           " --outSAMtype BAM Unsorted " 
           " --outSAMunmapped Within "
           " --outSAMstrandField intronMotif "
           " --outReadsUnmapped None "
           " --outFileNamePrefix {params.prefix}"     
           " --outFilterType BySJout "
           " --outFilterMismatchNmax 999 "
           " --outFilterMismatchNoverReadLmax 0.04 "
           " --alignIntronMin 10 "
           " --alignIntronMax 1000000 "
           " --alignMatesGapMax 1000000 " 
           " --alignSJoverhangMin 8 "
           " --alignSJDBoverhangMin 1 " 
           " --sjdbScore 1 "          
           " --limitOutSJcollapsed 5000000 "
        
# Add genome annotation
rule rna_feature_counts:
    input:
        bam = rules.rna_star.output.bam
    output:
        counts = temp('{sequencing_path}/{chunk}/04_rna_featureCounts.tsv.zst'),
        genes = temp('{sequencing_path}/{chunk}/04_rna_featureCounts.genes'), 
        summary = temp('{sequencing_path}/{chunk}/04_rna_featureCounts.genes.summary'),
    params:
        annot = config["genome"]["gene_annotation"],
        features = "gene",
        gene_name = "gene_id",
        tmp_bam = lambda w, input: input.bam + ".featureCounts.bam",
        script = srcdir("scripts/shareseq/featurecounts_to_tsv.py")
    threads: 16
    log: '{sequencing_path}/{chunk}/04_rna_featureCounts.log'
    shell: "featureCounts -T {threads} -Q 30 -a {params.annot} -t {params.features} -g {params.gene_name} -s 1"
           " -o {output.genes} -R CORE {input.bam} 2>{log}; " 
           " python {params.script} < {input.bam}.featureCounts | " # Convert counts output
           " zstd --fast=1 -q -o {output.counts}; "
           " rm {input.bam}.featureCounts "

rule rna_collapse_umis:
    input:
        counts = rules.rna_feature_counts.output.counts,
        script = rules.build_count_unique.output
    output:
        counts = temp('{sequencing_path}/{chunk}/05_umi_counts.tsv.zst')
    params:
        memory = "4G"
    threads: 4
    shell: "LC_ALL=C sort --parallel={threads} -S {params.memory} "
           " <(zstd -dc {input.counts}) | "
           " {input.script} | "
           " zstd --fast=1 -q -o {output.counts} "

def rna_group_barcodes_input(w):
    """Get subshell inputs to decompress and filter inputs by barcode group"""
    if "_" not in w.barcode_chunk:
        print("ERROR:", w)
    barcode_index = int(w.barcode_chunk.split("_")[1]) - 1
    if barcode_index >= len(end_barcode_groups):
        print("ERROR:", barcode_index, len(end_barcode_groups))
    grep_pattern = b"|".join(end_barcode_groups[barcode_index]).decode()    
    return expand_sublibrary_chunks(f"<(zstd -dc {rules.rna_collapse_umis.output.counts} | " +
                                 f" grep -E '\+({grep_pattern})\t')", "RNA", w)    


rule rna_group_barcodes:
    input:
        counts = lambda w: expand_sublibrary_chunks(rules.rna_collapse_umis.output.counts, "RNA", w), 
        script = rules.build_count_unique.output
    output:
        counts = temp('RNA/sublibraries/{sublibrary}/{barcode_chunk}/counts.tsv.zst')
    params:
        counts_decompress = rna_group_barcodes_input,
        memory = "4G"
    threads: 6
    shell: "LC_ALL=C sort --parallel={threads} --merge --batch-size=100 -S {params.memory} {params.counts_decompress} | "
           " {input.script} merge | "
           " zstd --fast=1 -q -o {output.counts} "

rule rna_dedup_count:
    input:
        counts = rules.rna_group_barcodes.output.counts
    output:
        counts = temp('RNA/sublibraries/{sublibrary}/{barcode_chunk}/counts_dedup.tsv.zst')
    params:
        script = srcdir("scripts/shareseq/run_umi_tools.py")
    log: temp('RNA/sublibraries/{sublibrary}/{barcode_chunk}/counts_dedup.log')
    shell: "zstd -dc {input.counts} | "
           " python {params.script} 2> {log} | " 
           " zstd --fast=1 -q -o {output.counts}"

rule rna_merge_sublibrary:
    input:
        counts = expand(rules.rna_dedup_count.output.counts, barcode_chunk=barcode_chunks, allow_missing=True)
    output:
        counts = temp('RNA/sublibraries/{sublibrary}/counts.tsv.gz')
    shell: "zstd -dc {input.counts} | gzip --fast > {output.counts}"
    
# generate 10x-compatible matrix per sublibrary   

# Feature list per sample and per sublibrary
localrules: rna_prep_features
rule rna_prep_features:
    input:
        annot = config["genome"]["gene_annotation"],
    output:
        features = temp('RNA/samples/features.tsv.gz'),
    params:
        script = srcdir("scripts/shareseq/rna_prep_feature_list.py")
    shell: "python {params.script} {input.annot} | gzip > {output.features}"

localrules: rna_features_sample
rule rna_features_sample:
    input: rules.rna_prep_features.output
    output: "RNA/samples/{sample}.features.tsv.gz"
    shell: "cp {input} {output}"

localrules: rna_features_sublibrary
rule rna_features_sublibrary:
    input: rules.rna_prep_features.output
    output: "RNA/sublibraries/{sublibrary}/features.tsv.gz"
    shell: "cp {input} {output}"

# Cell barcode list per sample and per sublibrary
rule rna_unique_cells_chunk:
    input:
        counts = rules.rna_feature_counts.output.counts
    output:
        cells = temp('{sequencing_path}/{chunk}/06_barcodes.txt.gz')
    shell: "zstd -dc {input.counts} | cut -f 2 | sort --unique | gzip > {output.cells}"

rule rna_unique_cells_sublibrary:
    input:
        counts = lambda w: expand_sublibrary_chunks(rules.rna_unique_cells_chunk.output.cells, "RNA", w) 
    output:
        cells = 'RNA/sublibraries/{sublibrary}/barcodes.tsv.gz'
    shell: "gzip -dc {input.counts} | sort --unique | gzip > {output.cells}"

    
rule rna_unique_cells_sublibrary_prefix:
    input: 
        cells = rules.rna_unique_cells_sublibrary.output.cells
    output:
        cells = temp('RNA/sublibraries/{sublibrary}/prefixed_barcodes.tsv.gz')
    params:
        sublibrary_id = lambda w: w.sublibrary
    shell: "gzip -dc {input.cells} | "
           "awk -c '{{print \"{params.sublibrary_id}_\" $0;}}' | "
           "gzip --fast > {output.cells}"

rule rna_unique_cells_sample:
    input:
        cells = lambda w: expand(rules.rna_unique_cells_sublibrary_prefix.output.cells, sublibrary = utils.get_sublibraries("RNA", config))
    output:
        cells = 'RNA/samples/{sample}.barcodes.tsv.gz'
    params:
        barcode_pattern = lambda w: f"_{config['samples'][w.sample]}\\+"
    shell: "gzip -dc {input.cells} | "
           "grep -E '{params.barcode_pattern}' | "
           "sort --unique | gzip > {output.cells}"
        
# mtx values -- collate per-sublibrary and per-sample
rule rna_mtx_chunk_sublibrary:
    input: 
        counts = rules.rna_dedup_count.output.counts,
        features = rules.rna_prep_features.output.features,
        barcodes = rules.rna_unique_cells_sublibrary.output.cells
    output:
        mtx = temp('RNA/sublibraries/{sublibrary}/{barcode_chunk}/mtx_entries.zst')
    params:
        script = srcdir("scripts/shareseq/mtx_from_counts.py"),
        memory = "4G"
    threads: 4
    shell: "zstd -dc {input.counts} | "
           "python {params.script} {input.features} {input.barcodes}  | "
           "LC_ALL=C sort -k2,2n -k1,1n -t$'\\t' -S {params.memory} --parallel=2 | "
           "zstd --fast=1 -q -o {output.mtx}"

rule rna_mtx_chunk_sample:
    input:
        counts = rules.rna_dedup_count.output.counts,
        features = rules.rna_prep_features.output.features,
        barcodes = rules.rna_unique_cells_sample.output.cells,
    output:
        mtx = temp('RNA/sublibraries/{sublibrary}/{barcode_chunk}/{sample}_mtx_entries.zst')
    params:
        script = srcdir("scripts/shareseq/mtx_from_counts.py"),
        memory = "4G",
        sublibrary_id = lambda w: w.sublibrary,
        barcode_pattern = lambda w: f"\t{config['samples'][w.sample]}\\+"
    threads: 4
    shell: "zstd -dc {input.counts} | "
           "grep -E '{params.barcode_pattern}' | " # Filter to sample
           "awk -c 'BEGIN {{OFS=\"\t\"}} {{$2=(\"{params.sublibrary_id}_\" $2); print $0;}}' | " # Prepend sublibrary ID
           "python {params.script} {input.features} {input.barcodes} | "
           "LC_ALL=C sort -k2,2n -k1,1n -t$'\\t' -S {params.memory} --parallel=2 | "
           "zstd --fast=1 -q -o {output.mtx}"

rule rna_mtx_merge_sublibrary:
    input:
        counts = expand(rules.rna_mtx_chunk_sublibrary.output.mtx, barcode_chunk=barcode_chunks, allow_missing=True),
    output:
        mtx = temp('RNA/sublibraries/{sublibrary}/mtx_entries.gz')
    params:
        mtx_decompress = expand(f"<(zstd -dc {rules.rna_mtx_chunk_sublibrary.output.mtx})", barcode_chunk=barcode_chunks, allow_missing=True),
        memory = "4G",
    threads: 4
    shell: "LC_ALL=C sort -k2,2n -k1,1n -t$'\\t' "
           " -S {params.memory} --batch-size=100 --parallel={threads} {params.mtx_decompress} | "
           " gzip -c > {output.mtx} "

rule rna_mtx_merge_sample:           
    input:
        counts = expand(rules.rna_mtx_chunk_sample.output.mtx, barcode_chunk=barcode_chunks, sublibrary=utils.get_sublibraries("RNA", config), allow_missing=True),
    output:
        mtx = temp('RNA/samples/{sample}.mtx_entries.gz')
    params:
        mtx_decompress = expand(f"<(zstd -dc {rules.rna_mtx_chunk_sample.output.mtx})", barcode_chunk=barcode_chunks, sublibrary=utils.get_sublibraries("RNA", config), allow_missing=True),
        memory = "4G",
    threads: 4
    shell: "LC_ALL=C sort -k2,2n -k1,1n -t$'\\t' "
           " -S {params.memory} --batch-size=100 --parallel={threads} {params.mtx_decompress} | "
           " gzip -c > {output.mtx} "

rule rna_mtx_sublibrary:
    input:
        mtx = rules.rna_mtx_merge_sublibrary.output.mtx,
        features = rules.rna_mtx_chunk_sublibrary.input.features,
        barcodes = rules.rna_mtx_chunk_sublibrary.input.barcodes,
    output:
        mtx = 'RNA/sublibraries/{sublibrary}/matrix.mtx.gz'
    params:
        script = srcdir("scripts/shareseq/mtx_add_header.py")
    threads: 2
    shell: "python {params.script} {input.mtx} {input.features} {input.barcodes} | "
           "gzip -c --fast > {output.mtx}"

rule rna_mtx_sample:
    input:
        mtx = rules.rna_mtx_merge_sample.output.mtx,
        features = rules.rna_mtx_chunk_sample.input.features,
        barcodes = rules.rna_mtx_chunk_sample.input.barcodes,
    output:
        mtx = 'RNA/samples/{sample}.matrix.mtx.gz'
    params:
        script = srcdir("scripts/shareseq/mtx_add_header.py")
    threads: 2
    shell: "python {params.script} {input.mtx} {input.features} {input.barcodes} | "
           "gzip -c --fast > {output.mtx}"

rule rna_stats_libraries:
    input: 
        barcode_stats = lambda w: expand_sublibrary_chunks(rules.match_barcodes.output.stats, "RNA", w),    
        star_log = lambda w: expand_sublibrary_chunks(rules.rna_star.log.summary, "RNA", w),
        feature_counts_log = lambda w: expand_sublibrary_chunks(rules.rna_feature_counts.output.summary, "RNA", w),
        dedup_log = expand(rules.rna_dedup_count.log, barcode_chunk=barcode_chunks, allow_missing=True)
    output: 
        barcodes = "RNA/sublibraries/{sublibrary}/barcode_stats.json",
        summary = "RNA/sublibraries/{sublibrary}/alignment_stats.json"
    params:
        script = srcdir("scripts/shareseq/stats_collect_rna.py")
    wildcard_constraints:
        sequencing_path = "RNA/.*"
    shell: "python {params.script} "
          " --barcode_stats {input.barcode_stats} "
          " --star_log {input.star_log} "
          " --feature_counts_log {input.feature_counts_log} "
          " --dedup_log {input.dedup_log} "
          " --barcodes_output {output.barcodes} "
          " --summary_output {output.summary} "
          
localrules: rna_stats_merge
rule rna_stats_merge:
    input: 
        sequencing = expand(rules.rna_stats_libraries.output.summary, sublibrary=utils.get_sublibraries("RNA", config)),
        barcodes = expand(rules.rna_stats_libraries.output.barcodes, sublibrary=utils.get_sublibraries("RNA", config)),
    output: 
        sequencing = "RNA/samples/alignment_stats.json",
        barcodes = "RNA/samples/barcode_stats.json"
    params: 
        script = srcdir("scripts/shareseq/stats_aggregate.py")
    shell: "python {params.script} --input {input.sequencing} --output {output.sequencing};"
           "python {params.script} --input {input.barcodes} --output {output.barcodes};"

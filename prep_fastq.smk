# Fastq preparation to run before main shareseq pipeline
# Author: Ben Parks
# Last Modified: 12/7/22

# 1. Demultiplex a run using bcl2fastq
# 2. Count the number of reads per fastq file
import utils

workdir: config["output_dir"]

wildcard_constraints:
    sequencing_path = "(ATAC|RNA)/([^/]+/)?[^/]+", # Sequencing path is 2-3 folders
    tile_chunk = "[0-9]+"

utils.string_only_keys(config)

outputs = (
    expand("{sequencing_path}/read_count.txt", sequencing_path=utils.get_sequencing_paths("ATAC", config) + utils.get_sequencing_paths("RNA", config)) + 
    expand("bcl2fastq/{sequencing_path}_{read}.fastq.zst", sequencing_path=utils.get_sequencing_paths("ATAC", config, run_types=["bcl"]) + utils.get_sequencing_paths("RNA", config, run_types="bcl"), read=["R1", "R2"])
)
filtered_outputs = []
for o in outputs:
    if os.path.exists(o):
        print(f"Skipping existing output: {o}", file=sys.stderr)
    else:
        filtered_outputs.append(o)



localrules: all
rule all:
    input: filtered_outputs
        

#############################
### Direct BCL extraction
#############################

localrules: bcl2fastq_samples
rule bcl2fastq_samples:
    output: 
        samples = "bcl2fastq/samples_{sequencing_run}.csv"
    script: "scripts/prep_fastq/prep_samplesheet.py"

rule bcl2fastq:
    input:
        run_dir = lambda w: config["sequencing"][w.sequencing_run]["run_dir"],
        samples = rules.bcl2fastq_samples.output,
    output:
        results_dir = temp(directory("bcl2fastq/raw/{sequencing_run}/{tile_chunk}")),
    params:
        script = srcdir("scripts/prep_fastq/run_bcl2fastq.py"),
        num_tile_chunks = lambda w: config["sequencing"][w.sequencing_run]["tile_chunks"]
    resources:
        runtime = 4 * 60,
        mem_mb = 32_000,
    threads: 16
    log: "bcl2fastq/logs/{sequencing_run}_{tile_chunk}.log"
    benchmark: "bcl2fastq/logs/{sequencing_run}_{tile_chunk}.runtime"
    shell: "python {params.script}"
        " --samples {input.samples} "
        " --input {input.run_dir} "
        " --output {output.results_dir} "
        " --threads {threads} "
        " --num_tile_chunks {params.num_tile_chunks} "
        " --tile_chunk {wildcards.tile_chunk} "
        " --log {log}"

def bcl2fastq_dependency(sequencing_path, tile_chunk):
    """Take a sublibrary path and return a list of input dependencies"""
    run_id = sequencing_path.split("/")[1]
    assert config["sequencing"][run_id]["type"] == "bcl"
    return f"bcl2fastq/raw/{run_id}/{tile_chunk}"


def bcl2fastq_output(sequencing_path, tile_chunk, read):
    """Take a bcl sublibrary path and return the path to its bcl2fastq path"""
    assay_type, run_id = sequencing_path.split("/")[:2]
    assert config["sequencing"][run_id]["type"] == "bcl"
    sublib_id = sequencing_path.split("/")[2]
    return f"bcl2fastq/raw/{run_id}/{tile_chunk}/{assay_type}_{sublib_id}_{read}.fastq.gz"


localrules: build_fastq_index_to_readname
rule build_fastq_index_to_readname:
    input: srcdir("scripts/prep_fastq/fastq_index_to_readname.c")
    output: "bin/fastq_index_to_readname"
    shell: "gcc -O3 -o {output} {input}"

rule bcl2fastq_index_to_read_names:
    input: 
        fastqs = lambda w: bcl2fastq_dependency(w.sequencing_path, w.tile_chunk),
        script = rules.build_fastq_index_to_readname.output,
    output:
        R1 = temp("bcl2fastq/{sequencing_path}/{tile_chunk}/R1.fastq.zst"),
        R2 = temp("bcl2fastq/{sequencing_path}/{tile_chunk}/R2.fastq.zst"),
    params: 
        R1 = lambda w: bcl2fastq_output(w.sequencing_path, w.tile_chunk, "R1"),
        R2 = lambda w: bcl2fastq_output(w.sequencing_path, w.tile_chunk, "R2"),
        I1 = lambda w: bcl2fastq_output(w.sequencing_path, w.tile_chunk, "I1"),
        I2 = lambda w: bcl2fastq_output(w.sequencing_path, w.tile_chunk, "I2"),
    resources:
        runtime = 2 * 60, # Be generous on time in case of large fastqs
    threads: 3
    log: "bcl2fastq/logs/{sequencing_path}/index_to_read_names_{tile_chunk}.log"
    shell: "{input.script} <(gzip -dc {params.R1}) <(gzip -dc {params.R2}) "
           "  <(gzip -dc {params.I1}) <(gzip -dc {params.I2}) "
           "  {output.R1} {output.R2} "
           " --output-cmd 'zstd --fast -qo $FILE' 2> {log}"

def get_tile_chunks(sequencing_path):
    """Generate tile chunk IDs for a sequencing path. Adds padding 0s as needed"""
    run_id = sequencing_path.split("/")[1]
    assert config["sequencing"][run_id]["type"] == "bcl"
    tile_chunks = config["sequencing"][run_id]["tile_chunks"]
    str_len = max(2, len(str(tile_chunks)))
    return [f"{i:0{str_len}d}" for i in range(1, tile_chunks+1)]

rule bcl2fastq_merge_tile_chunks:
    input:
        fastq = lambda w: expand("bcl2fastq/{sequencing_path}/{tile_chunk}/{read}.fastq.zst", tile_chunk=get_tile_chunks(w.sequencing_path), allow_missing=True),
    output:
        fastq = "bcl2fastq/{sequencing_path}_{read}.fastq.zst",
    threads: 4
    log: "bcl2fastq/logs/{sequencing_path}/merge_tile_chunks_{read}.log"
    shell: "zstd -dc {input.fastq} | pzstd -1 -c --processes 3 -o {output.fastq} 2> {log}"

#############################
### Count reads
#############################

rule count_reads:
    input: lambda w: utils.fastq_path(w.sequencing_path, "R1", config)
    output: "{sequencing_path}/read_count.txt"
    params: 
        decompress = lambda w: utils.fastq_decompress(w.sequencing_path, config),
        truncate_test_chunks = lambda w: f" | head -n {config['chunk_size']*config['test_chunks']*4} " if "test_chunks" in config else ""
    resources:
        runtime = 5 * 60, # Be generous on time in case of large fastqs
    shell: "{params.decompress} {input} {params.truncate_test_chunks} | awk -c 'END{{print int(NR/4)}}' > {output}"
           "|| if [[ $? -eq 141 ]]; then true; else exit $?; fi" # Ignore spurious exit code 141 produced with test_chunks set


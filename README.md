# shareseq-pipeline
Snakemake-based pipeline for processing SHARE-seq data 

> [!NOTE]
> This branch is dedicated to processing anonymized raw sequence data as produced by the `anonymize` branch. It **is not synchronized with the main branch**.


## Processing anonymized SHARE-seq data from SRA

This branch processes already-anonymized SHARE-seq FASTQs into per-sample
fragments and matrices, with optional SRA download as the first stage. The
pipeline runs `prep_sra.smk` → `ingest_anonymized.smk` → `shareseq.smk`,
bypassing `bcl2fastq` and the demultiplexing/anonymization stages on `main`.
We use data from HDMA (Liu et al, Nature, 2026 as an exampl).

### Inputs
- SRR accessions for the anonymized FASTQs (or the FASTQs already on disk)
- Genome references: bowtie2 + STAR indexes, GTF annotation, FASTA. Use
  `bash scripts/references/prep_genome.sh hg38` to build these if needed.
- A per-batch config YAML — see [runs/share_sra_demo.yaml](runs/share_sra_demo.yaml).

### Outputs (per sample, under `output_dir`)
- ATAC fragments: `ATAC/samples/{sample}.fragments.tsv.gz`
- RNA matrix:     `RNA/samples/{sample}.{matrix.mtx,barcodes.tsv,features.tsv}.gz`
- QC stats:       `{ATAC,RNA}/samples/{alignment,barcode}_stats.json`
- Per-sublibrary fragments/matrices under `{ATAC,RNA}/sublibraries/`.

### Running on Sherlock

This section explains how to run the process on Stanford's Sherlock HPC; it 
can be adapted for other settings.

> [!NOTE]
> For general Sherlock setup (dependencies, profile config, containerized
> runs, building genome references), see the `main` branch
> [README](https://github.com/GreenleafLab/shareseq-pipeline/blob/main/README.md).
> The steps below cover only what's specific to the anonymized + SRA workflow.

1. Load conda environment and modules (e.g. sra-toolkit).

2. Copy [runs/share_sra_demo.yaml](runs/share_sra_demo.yaml) to a new
   `runs/MY_CONFIG.yaml` and edit:
   - `output_dir` — pipeline working directory
   - `samples` — Round1 BC1 regex per sample (must cover all 96 barcodes
     exactly once across entries; see [shareseq.smk:69-72](shareseq.smk#L69-L72))
   - `sequencing.<run>.data_dir` — where SRA fetches will land
   - `sequencing.<run>.{ATAC,RNA}_sra` — `{sample: SRR}` for samples to
     fetch from SRA (omit for samples already in `data_dir`)
   - `sequencing.<run>.{ATAC,RNA}_samples` — every sample to process

  `runs/share_sra_demo.yaml` is an example for downloading and processing ATAC and RNA runs
  for one sample from one batch in the Human Development Multiomic Atlas (HDMA).

3. Do a small test run with truncated inputs:
   ```yaml
   chunk_size: 2_000_000
   test_chunks: 2
   ```
   ```bash
   sbatch -p wjg,sfgf,biochem run_process_anonymize.sh runs/MY_CONFIG.yaml
   ```

4. Full pass: remove `test_chunks` and use `chunk_size: 20_000_000`, then
   re-run the same `sbatch` command.

`run_process_anonymize.sh` will skip the SRA stage when no `*_sra`
keys are declared, so it also works for FASTQs that are already present locally
(named `{ATAC|RNA}_{SampleID}_anon_{R1|R2}.fastq.gz` in `data_dir`).

**Cleanup of intermediates** (optional — only needed if `run_process_anonymize.sh`
was run with `--notemp`):
```bash
snakemake --profile=$(pwd)/profile -s ingest_anonymized.smk --configfile runs/MY_CONFIG.yaml --delete-temp-output --config filter_dag=false
snakemake --profile=$(pwd)/profile -s shareseq.smk          --configfile runs/MY_CONFIG.yaml --delete-temp-output --config filter_dag=false
```

**Containerized runs:** uncomment the `singularity:` line in the config
(pointing to a `shareseq_latest.sif` built from
`docker://bettybliu/shareseq:latest`); the wrapper picks it up automatically.

---


## Required dependencies
- bcl2fastq
- bgzip
- bowtie2 version >=2.4.2 (flag --sam-append-comment)
- fastp
- featureCounts (subread)
- pysam
- python3
- STAR
- samtools
- snakemake
- sra-tools (`prefetch`, `fastq-dump`) — only needed if running `prep_sra.smk`
- tabix
- umi_tools
- zstd
- BAMboozle
- seqkit

Builtin unix tools:
- awk
- gcc
- grep
- gzip
- sort
- split

### Additional dependencies needed for plotting
- poppler
- R
    - Seurat
    - ggplot2
    - patchwork
    - ggrastr
    - gridExtra
    - dplyr
    - rjson
    - [BPCells](https://bnprks.github.io/BPCells/index.html)


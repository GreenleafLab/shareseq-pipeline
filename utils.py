import subprocess

# Helper functions for snakemake pipeline logic
# Authors: Ben Parks, Betty Liu
# Last Modified: 12/7/22

#####################################
# Normalize input config object
#####################################

def string_only_keys(data):
    """Recursively edit a dictionary to convert all keys to string types"""
    assert isinstance(data, dict)
    for k, v in list(data.items()):
        if not isinstance(k, str):
            data[str(k)] = v
            del data[k]
        if isinstance(v, dict):
            string_only_keys(v)

#####################################
# Gather Pipeline Inputs
#####################################
        
def get_sequencing_paths(assay, config, run_types=["bcl", "anonymized"], sublib=None):
    """Get a list of all paths for sublibraries for the current assay.

    For type "bcl" runs, returns paths of the form "{assay}/{run_id}/{sublib_id}".
    For type "anonymized" runs (already-demuxed-per-sample SRA inputs), returns
    "{assay}/samples/{sample_id}" — each SampleID is treated as its own sublibrary.
    """
    assert assay in ["ATAC", "RNA"]
    sequencing_paths = []
    for run_id, run in config["sequencing"].items():
        assert run["type"] in ["bcl", "anonymized"]
        if run["type"] == "bcl" and "bcl" in run_types:
            if (f"{assay}_I2" in run.keys()) and run[f"{assay}_I2"]:
                whitelist = run[f"{assay}_I2"] if not sublib else [sublib]
                sequencing_paths += [
                    f"{assay}/{run_id}/{sublib_id}" for sublib_id in run[f"{assay}_I2"] if sublib_id in whitelist
                ]
        elif run["type"] == "anonymized" and "anonymized" in run_types:
            if (f"{assay}_samples" in run.keys()) and run[f"{assay}_samples"]:
                whitelist = run[f"{assay}_samples"] if not sublib else [sublib]
                sequencing_paths += [
                    f"{assay}/samples/{sample_id}" for sample_id in run[f"{assay}_samples"] if sample_id in whitelist
                ]
    return sequencing_paths

def get_sublibraries(assay, config, run_types=["bcl", "anonymized"]):
    """Get a list of all unique sublibraries for the current assay.

    For type "anonymized" runs, each SampleID listed in {assay}_samples is
    treated as a sublibrary so downstream sublibrary-level rules in
    shareseq.smk operate per sample.
    """
    assert assay in ["ATAC", "RNA"]
    sublibraries = []
    for run_id, run in config["sequencing"].items():
        assert run["type"] in ["bcl", "anonymized"]
        if run["type"] == "bcl" and "bcl" in run_types:
            if (f"{assay}_I2" in run.keys()) and run[f"{assay}_I2"]:
                sublibraries += [
                    sublib_id for sublib_id in run[f"{assay}_I2"]
                ]
        elif run["type"] == "anonymized" and "anonymized" in run_types:
            if (f"{assay}_samples" in run.keys()) and run[f"{assay}_samples"]:
                sublibraries += [
                    sample_id for sample_id in run[f"{assay}_samples"]
                ]
    return list(set(sublibraries))

def fastq_path(sequencing_path, read, config):
    """Take a sublibrary path and return the path to its R1 or R2 fastq.

    BCL inputs land in bcl2fastq/ (produced by prep_fastq.smk); anonymized
    inputs land in staged_fastq/ (produced by ingest_anonymized.smk).
    """
    parts = sequencing_path.split("/")
    if parts[1] == "samples":
        return f"staged_fastq/{sequencing_path}_{read}.fastq.zst"
    run_id = parts[1]
    if config["sequencing"][run_id]["type"] == "bcl":
        return f"bcl2fastq/{sequencing_path}_{read}.fastq.zst"
    assert False

def fastq_decompress(sequencing_path, config):
    """Take a sublibrary path and return the command to decompress it"""
    parts = sequencing_path.split("/")
    if parts[1] == "samples":
        return "zstd -dc"
    run_id = parts[1]
    if config["sequencing"][run_id]["type"] == "bcl":
        return "zstd -dc"
    assert False

#####################################
# Logic for barcode-based sample demultiplexing
#####################################

def grep_regex_match(text, regex):
    """Check if a regex matches text, using grep -E as the regex engine"""
    if isinstance(text, str):
        text = text.encode()
    res = subprocess.run(["grep", "-E", f"^{regex}$"], input=text, stdout=subprocess.PIPE)
    return len(res.stdout) != 0

def bc_names(tsv_path):
    """Return list of barcode names from a barcode file"""
    lines = open(tsv_path, "rb").readlines()
    assert lines[0] == b"Name\tSequence\n"
    names = []
    for l in lines[1:]:
        names.append(l.strip().split(b"\t")[0])
    return names

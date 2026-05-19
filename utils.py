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

    For both "bcl" and "anonymized" run types, returns paths of the form
    "{assay}/{run_id}/{sublib_id}". For "anonymized" runs the sublibraries
    are taken from run["sublibraries"][assay], a {sublib_id: I2_sequence} map
    that ingest_anonymized.smk uses to split-by-I2 the merged per-sample SRA
    upload back into per-sublibrary FASTQs. The shared path scheme means
    shareseq.smk treats anon and bcl sublibraries identically and downstream
    cell barcodes inherit the original CL{N}_ prefix.
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
            sublibs = (run.get("sublibraries") or {}).get(assay) or {}
            if sublibs:
                whitelist = sublibs if not sublib else [sublib]
                sequencing_paths += [
                    f"{assay}/{run_id}/{sublib_id}" for sublib_id in sublibs if sublib_id in whitelist
                ]
    return sequencing_paths

def get_sublibraries(assay, config, run_types=["bcl", "anonymized"]):
    """Get a list of all unique sublibraries for the current assay.

    For "anonymized" runs the sublibraries come from run["sublibraries"][assay].
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
            sublibs = (run.get("sublibraries") or {}).get(assay) or {}
            sublibraries += list(sublibs.keys())
    return list(set(sublibraries))

def fastq_path(sequencing_path, read, config):
    """Take a sublibrary path and return the path to its R1 or R2 fastq.

    BCL inputs land in bcl2fastq/ (produced by prep_fastq.smk); anonymized
    inputs land in staged_fastq/ (produced by ingest_anonymized.smk).
    """
    parts = sequencing_path.split("/")
    run_id = parts[1]
    rtype = config["sequencing"][run_id]["type"]
    if rtype == "bcl":
        return f"bcl2fastq/{sequencing_path}_{read}.fastq.zst"
    if rtype == "anonymized":
        return f"staged_fastq/{sequencing_path}_{read}.fastq.zst"
    assert False

def fastq_decompress(sequencing_path, config):
    """Take a sublibrary path and return the command to decompress it"""
    parts = sequencing_path.split("/")
    run_id = parts[1]
    rtype = config["sequencing"][run_id]["type"]
    if rtype in ("bcl", "anonymized"):
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

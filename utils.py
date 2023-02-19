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
        
def get_sequencing_paths(assay, config, run_types=["bcl"]):
    """Get a list of all paths for sublibraries for the current sample"""
    assert assay in ["ATAC", "RNA"]
    sequencing_paths = []
    for run_id, run in config["sequencing"].items():
        assert run["type"] in ["bcl"]
        if run["type"] == "bcl" and "bcl" in run_types:
            if (f"{assay}_I2" in run.keys()) and run[f"{assay}_I2"]:
		sequencing_paths += [
                	f"{assay}/{run_id}/{sublib_id}" for sublib_id in run[f"{assay}_I2"]
            	]
    return sequencing_paths

def fastq_path(sequencing_path, read, config):
    """Take a sublibrary path and return the path to its R1 or R2 fastq"""
    assay_type, run_id = sequencing_path.split("/")[:2] 
    if config["sequencing"][run_id]["type"] == "bcl":
        return f"bcl2fastq/{sequencing_path}_{read}.fastq.zst"
    else:
        assert False

def fastq_decompress(sequencing_path, config):
    """Take a sublibrary path and return the command to decompress it"""
    assay_type, run_id = sequencing_path.split("/")[:2] 
    if config["sequencing"][run_id]["type"] == "bcl":
        return "zstd -dc" 
    else:
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

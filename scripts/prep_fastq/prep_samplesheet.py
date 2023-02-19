# Write samples for run_bcl2fastq.py
sequencing_run = snakemake.config["sequencing"][snakemake.wildcards.sequencing_run]
assert sequencing_run["type"] == "bcl"

def get_items(assay):
    if f"{assay}_I2" in sequencing_run.keys() and sequencing_run[f"{assay}_I2"]:
        return sequencing_run[f"{assay}_I2"].items()
    else:
        return []

# Make samplesheet
with open(snakemake.output.samples, "w") as out:
    out.write("Sample, I2\n")
    for sample_id, I2 in get_items("ATAC"):
        out.write(f"ATAC_{sample_id},{I2}\n")
    for sample_id, I2 in get_items("RNA"):
        out.write(f"RNA_{sample_id},{I2}\n")

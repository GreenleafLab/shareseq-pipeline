# Write samples for run_bcl2fastq.py
sequencing_run = snakemake.config["sequencing"][snakemake.wildcards.sequencing_run]
assert sequencing_run["type"] == "bcl"

# Make samplesheet
with open(snakemake.output.samples, "w") as out:
    out.write("Sample, I2\n")
    for sample_id, I2 in sequencing_run["ATAC_I2"].items():
        out.write(f"ATAC_{sample_id},{I2}\n")
    for sample_id, I2 in sequencing_run["RNA_I2"].items():
        out.write(f"RNA_{sample_id},{I2}\n")

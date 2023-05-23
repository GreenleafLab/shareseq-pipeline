import sys
import pysam

# Convert bam file to a fragment file format, while adding +4/-4 coordinate adjustment

# Authors: Ben Parks
# Last updated: 9/26/22

# Usage: python bam_to_fragments.py file.bam > fragments.tsv

input = pysam.AlignmentFile(sys.argv[1], "rb")

prev_read = next(input, None)
while True:
    read = next(input, None)
    if read is None or prev_read is None:
        break
    if read.query_name != prev_read.query_name:
        prev_read = read
        continue

    if prev_read.is_reverse:
        prev_read, read = read, prev_read

    chromosome = read.reference_name
    start = prev_read.reference_start + 4
    end = prev_read.reference_start + prev_read.template_length - 4
    cell_barcode = read.get_tag("CB")
    print(chromosome, start, end, cell_barcode, sep="\t")
    prev_read = next(input, None)


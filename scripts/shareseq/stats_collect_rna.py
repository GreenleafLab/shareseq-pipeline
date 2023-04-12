import argparse
import glob
import gzip
import json
import re
import sys

def main():
    parser = argparse.ArgumentParser(description="Collect statistics from RNA pipeline")
    parser.add_argument("--barcode_stats", type=str, nargs="+", help="Paths of input barcode matching stats")
    parser.add_argument("--star_log", type=str, nargs="+", help="Paths of STAR logs *final.out")
    parser.add_argument("--feature_counts_log", type=str, nargs="+", help="Paths of featureCounts logs")
    parser.add_argument("--dedup_log", type=str, nargs="+", help="Path of UMItools count log")
    parser.add_argument("--summary_output", type=str, help="Path of output json summary file")
    parser.add_argument("--barcodes_output", type=str, help="Path of detailed barcode matching summary file")
    args = parser.parse_args()

    ### Barcode-matching stats
    barcode_stats = {}
    for f in args.barcode_stats:
        stats = json.load(open(f))
        barcode_stats = add_stat_counts(barcode_stats, stats)

    total_reads = sum(barcode_stats["total_mismatch_histogram"])
    matching_reads = sum(barcode_stats["total_mismatch_histogram"][:-1])

    ### STAR stats
    star_stats = {}
    for f in args.star_log:
        star_stats = add_stat_counts(star_stats, parse_star_log(f))

    ### featureCounts stats
    feature_counts_stats = {}
    for f in args.feature_counts_log:
        feature_counts_stats = add_stat_counts(feature_counts_stats, parse_feature_counts_log(f))
    
    ### UMItools dedup and count stats
    dedup_stats = {}
    for f in args.dedup_log:
        stats = parse_dedup_log(f)
        dedup_stats = add_stat_counts(dedup_stats, stats)
 
    ### Alignment summary
    json.dump({
        "total_reads": total_reads,
        "valid_barcode_reads": matching_reads,
        "total_aligned_reads": star_stats["unique_map_reads"] + star_stats["total_multimap"] + star_stats["total_multimap_toomany"],
        "multimapped_reads": star_stats["total_multimap"] + star_stats["total_multimap_toomany"],
        "total_annotated_reads": feature_counts_stats["Assigned"],
        "total_deduped_reads": dedup_stats["deduped_reads"],
    }, open(args.summary_output, "w"), indent = 2)

    json.dump(barcode_stats, open(args.barcodes_output, "w"), indent=2)

def add_stat_counts(d1, d2):
    """
    Merge to stat dictionaries together, where all corresponding integer
    values are added
    """
    out = {}
    for k, v in d1.items():
        if k not in d2:
            out[k] = v
        elif isinstance(v, int) or isinstance(v, float):
            out[k] = v + d2[k]
        elif isinstance(v, dict):
            out[k] = add_stat_counts(v, d2[k])
        elif isinstance(v, list):
            assert len(d2[k]) == len(v)
            out[k] = v
            for i in range(len(v)):
                v[i] += d2[k][i]
    for k, v in d2.items():
        if k not in d1:
            out[k] = v
    return out

#############################
### STAR log parsing
#############################

percentage_regex = "\d+\.\d+%"

star_log_template = r"""Number of input reads \|	(?P<total_reads>\d+)
                      Average input read length \|	\d+
                                    UNIQUE READS:
                   Uniquely mapped reads number \|	(?P<unique_map_reads>\d+)
                        Uniquely mapped reads % \|	{percentage_regex}
                          Average mapped length \|	\d+\.\d+
                       Number of splices: Total \|	(?P<total_splices>\d+)
            Number of splices: Annotated \(sjdb\) \|	(?P<total_splices_annot>\d+)
                       Number of splices: GT\/AG \|	(?P<splices_GT_AG>\d+)
                       Number of splices: GC\/AG \|	(?P<splices_GC_AG>\d+)
                       Number of splices: AT\/AC \|	(?P<splices_AT_AC>\d+)
               Number of splices: Non-canonical \|	(?P<total_splices_noncanonical>\d+)
                      Mismatch rate per base, % \|	{percentage_regex}
                         Deletion rate per base \|	{percentage_regex}
                        Deletion average length \|	\d+\.\d+
                        Insertion rate per base \|	{percentage_regex}
                       Insertion average length \|	\d+\.\d+
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci \|	(?P<total_multimap>\d+)
             % of reads mapped to multiple loci \|	{percentage_regex}
        Number of reads mapped to too many loci \|	(?P<total_multimap_toomany>\d+)
             % of reads mapped to too many loci \|	{percentage_regex}""".format(percentage_regex=percentage_regex)

def parse_star_log(path):
    output_groups = ["total_reads", "unique_map_reads", "total_multimap", "total_multimap_toomany",]
    log = open(path).read()
    match = re.search(star_log_template, log)
    if match is None:
        raise Exception(f"Path {path} doesn't match template regex for STAR log format")
    return {group: int(match.group(group)) for group in output_groups}


#############################
### featureCounts log parsing
#############################

feature_counts_log_template = """Assigned	(?P<Assigned>\d+)
Unassigned_Unmapped	(?P<Unassigned_Unmapped>\d+)
Unassigned_Read_Type	(?P<Unassigned_Read_Type>\d+)
Unassigned_Singleton	(?P<Unassigned_Singleton>\d+)
Unassigned_MappingQuality	(?P<Unassigned_MappingQuality>\d+)
Unassigned_Chimera	(?P<Unassigned_Chimera>\d+)
Unassigned_FragmentLength	(?P<Unassigned_FragmentLength>\d+)
Unassigned_Duplicate	(?P<Unassigned_Duplicate>\d+)
Unassigned_MultiMapping	(?P<Unassigned_MultiMapping>\d+)
Unassigned_Secondary	(?P<Unassigned_Secondary>\d+)
Unassigned_NonSplit	(?P<Unassigned_NonSplit>\d+)
Unassigned_NoFeatures	(?P<Unassigned_NoFeatures>\d+)
Unassigned_Overlapping_Length	(?P<Unassigned_Overlapping_Length>\d+)
Unassigned_Ambiguity	(?P<Unassigned_Ambiguity>\d+)"""

def parse_feature_counts_log(path):
    output_groups = ["Assigned", "Unassigned_Unmapped", "Unassigned_Read_Type", "Unassigned_Singleton",
                    "Unassigned_MappingQuality", "Unassigned_Chimera", "Unassigned_FragmentLength",
                    "Unassigned_Duplicate", "Unassigned_MultiMapping", "Unassigned_Secondary",
                    "Unassigned_NonSplit", "Unassigned_NoFeatures", "Unassigned_Overlapping_Length",
                    "Unassigned_Ambiguity"]
    log = open(path).read()
    match = re.search(feature_counts_log_template, log)
    if match is None:
        raise Exception(f"Path {path} doesn't match template regex for featureCounts log format")
    return {group: int(match.group(group)) for group in output_groups}

#############################
### UMItools count log parsing
#############################

dedup_log_template = """INFO Number of \(post deduplication\) reads counted: (?P<deduped_reads>\d+)"""

def parse_dedup_log(path):
    output_groups = ["deduped_reads"]
    log = open(path).read()
    match = re.search(dedup_log_template, log)
    if match is None:
        raise Exception(f"Path {path} doesn't match template regex for UMItools count log format")
    return {"deduped_reads": int(match.group("deduped_reads"))}


if __name__ == "__main__":
    main()

import json
from pathlib import Path

import numpy as np
import pandas as pd

import tqdm

p = Path(".snakemake/metadata")

raw_data = []
for f in tqdm.tqdm(list(p.iterdir())):
    raw_data.append(json.load(f.open()))

d = pd.DataFrame(raw_data)
d["duration"] = d.endtime - d.starttime

stats = d.groupby("rule").agg(
    count = ("starttime", len),
    mean_time = ("duration", "mean"),
    total_time = ("duration", "sum"),
    longest = ("duration", "max")
).reset_index().sort_values("mean_time", ascending=False, ignore_index=True)
stats["count"] = stats["count"].astype(int)

def format_duration(x):
    if not type(x) is float:
        return x
    x = int(x)
    hours = x // (60*60)
    minutes = (x % (60*60)) // 60
    seconds = x % 60
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"

print(str(stats.applymap(format_duration)))

#                              rule  count mean_time total_time   longest
# 0               rna_counts_to_mtx     12  00:40:22   08:04:27  00:45:01
# 1              atac_split_samples      4  00:33:25   02:13:41  00:46:50
# 2               atac_merge_chunks      4  00:30:03   02:00:15  00:55:02
# 3                    split_fastqs     16  00:28:49   07:41:09  00:49:35
# 4                     count_reads      8  00:20:22   02:42:58  00:31:44
# 5                    atac_bowtie2    288  00:20:04   96:21:54  00:38:40
# 6     bcl2fastq_merge_tile_chunks     16  00:17:40   04:42:43  00:35:51
# 7                       bcl2fastq     20  00:17:17   05:45:46  00:21:02
# 8              rna_feature_counts    729  00:11:18  137:23:51  00:13:09
# 9            atac_stats_libraries      8  00:11:02   01:28:22  00:15:53
# 10                       rna_star    729  00:10:41  129:48:43  00:13:25
# 11                 match_barcodes   1593  00:10:37  282:01:50  00:20:17
# 12              rna_collapse_umis    243  00:10:26   42:16:49  00:13:12
# 13             atac_trim_adapters    864  00:10:08  145:59:49  00:18:58
# 14              rna_trim_adapters    729  00:08:46  106:36:02  00:11:31
# 15  bcl2fastq_index_to_read_names    320  00:08:17   44:10:44  00:14:16
# 16         atac_convert_fragments    288  00:08:06   38:55:28  00:13:09
# 17              rna_split_samples      4  00:06:03   00:24:15  00:06:39
# 18             rna_group_barcodes     48  00:06:03   04:50:36  00:09:25
# 19              rna_merge_samples      1  00:06:01   00:06:01  00:06:01
# 20                rna_dedup_count     48  00:04:41   03:45:17  00:05:55
# 21           rna_merge_sublibrary      4  00:03:19   00:13:18  00:04:18
# 22             build_count_unique      1  00:00:11   00:00:11  00:00:11
# 23              rna_prep_features      1  00:00:10   00:00:10  00:00:10
# 24              bcl2fastq_samples      1  00:00:02   00:00:02  00:00:02
# 25               atac_stats_merge      1  00:00:01   00:00:01  00:00:01
# 26  build_fastq_index_to_readname      1  00:00:00   00:00:00  00:00:00


#                                   rule  count mean_time total_time   longest
# 0                    atac_merge_chunks      9  00:35:37   05:20:33  00:38:21
# 1                         split_fastqs     16  00:31:39   08:26:31  00:46:38
# 2                   atac_split_samples      3  00:28:09   01:24:28  00:29:58
# 3                         atac_bowtie2    288  00:20:44   99:33:34  00:37:25
# 4                          count_reads      8  00:20:42   02:45:40  00:33:07
# 5                            bcl2fastq     20  00:18:03   06:01:11  00:20:58
# 6          bcl2fastq_merge_tile_chunks     16  00:16:50   04:29:33  00:33:11
# 7                             rna_star    729  00:15:22  186:48:29  00:28:08
# 8                 rna_mtx_merge_sample      1  00:14:49   00:14:49  00:14:49
# 9                   rna_feature_counts    729  00:12:58  157:39:00  00:17:21
# 10                   rna_collapse_umis    243  00:12:12   49:25:12  00:16:37
# 11              rna_unique_cells_chunk    243  00:12:10   49:19:34  00:16:38
# 12                      rna_mtx_sample      1  00:12:01   00:12:01  00:12:01
# 13                  atac_trim_adapters    864  00:11:44  168:58:46  00:17:21
# 14              atac_convert_fragments    288  00:10:21   49:41:44  00:17:02
# 15                      match_barcodes   1593  00:09:22  248:49:54  00:13:46
# 16                atac_stats_libraries      6  00:09:15   00:55:35  00:10:31
# 17       bcl2fastq_index_to_read_names    320  00:08:49   47:04:08  00:32:30
# 18                   rna_trim_adapters    729  00:08:39  105:14:56  00:14:31
# 19         rna_unique_cells_sublibrary      4  00:04:24   00:17:37  00:05:24
# 20            rna_mtx_merge_sublibrary      4  00:04:20   00:17:21  00:04:32
# 21                     rna_dedup_count     48  00:04:19   03:27:55  00:06:20
# 22                  rna_group_barcodes     48  00:04:18   03:26:42  00:05:06
# 23                  rna_mtx_sublibrary      4  00:03:49   00:15:17  00:04:47
# 24  rna_unique_cells_sublibrary_prefix      4  00:03:47   00:15:08  00:03:50
# 25                 rna_stats_libraries      8  00:02:09   00:17:16  00:02:14
# 26            rna_mtx_chunk_sublibrary     48  00:01:38   01:19:02  00:02:07
# 27                rna_mtx_chunk_sample     48  00:01:28   01:10:54  00:01:45
# 28             rna_unique_cells_sample      1  00:01:25   00:01:25  00:01:25
# 29                   rna_prep_features      1  00:00:17   00:00:17  00:00:17
# 30                  build_count_unique      1  00:00:16   00:00:16  00:00:16
# 31                   bcl2fastq_samples      1  00:00:02   00:00:02  00:00:02
# 32       build_fastq_index_to_readname      1  00:00:01   00:00:01  00:00:01
# 33                     rna_stats_merge      2  00:00:00   00:00:00  00:00:00
# 34                 rna_features_sample      1  00:00:00   00:00:00  00:00:00
# 35             rna_features_sublibrary      4  00:00:00   00:00:00  00:00:00
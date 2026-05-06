#!/usr/bin/env python3
"""
Collapse UMIs per genomic position using umi-tools, including MAPQ aggregation.

For each unique combination of the first 5 columns (genomic position),
UMIs are clustered using umi-tools (directional method).
One output row is produced per UMI cluster, using the most abundant UMI
as representative.

MAPQ is aggregated as a read-count-weighted mean.

Input columns (tab-delimited):
  1 chrom
  2 start
  3 end
  4 orientation
  5 strand
  6 UMI
  7 read_count
  8 MAPQ

Output columns:
  chrom, start, end, orientation, strand
  UMI                  : representative (most abundant UMI)
  reads_sum            : total reads in cluster
  collapsed_UMIs       : number of UMIs collapsed
  mean_MAPQ            : weighted mean MAPQ of the cluster
  umi_list             : comma-separated list of collapsed UMIs
  umi_counts           : read counts per UMI
  umi_mapqs            : mean MAPQ per UMI
"""

import argparse
from pathlib import Path
from collections import defaultdict
import pandas as pd
from umi_tools import UMIClusterer


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Collapse UMIs per position with MAPQ using umi-tools"
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input TSV/BED file with UMI, read count, and MAPQ"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output TSV file with collapsed UMIs"
    )

    parser.add_argument(
        "-t", "--threshold",
        type=int,
        default=1,
        help="Edit-distance threshold for UMI clustering (default: 1)"
    )

    parser.add_argument(
        "-r", "--rescueR2",
        default="FALSE",
        help="whether to rescued R2 reads from pairs that failed to align"
    )

    parser.add_argument(
        "--method",
        default="directional",
        help="whether to rescued R2 reads from pairs that failed to align"
    )

    return parser.parse_args()

### Define function to read from R2 rescued reads align as single ends
# this will update the number of read per UMI per position detected from paired-end alignments

def load_file_into_result(filename, result):
    """
    Read a UMI/count/MAPQ file and populate or update the result dict.
    New positions and new UMIs are created if they do not already exist.
    """

    with open(filename) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")

            group_key = tuple(cols[0:5])
            umi = cols[5].encode("ascii")

            read_count = int(cols[6])
            mapq = float(cols[7])

            # Create new position + UMI if needed
            if umi not in result[group_key]:
                result[group_key][umi] = {
                    "count": 0,
                    "mapq_sum": 0,
                    "mapq_n": 0
                }

            result[group_key][umi]["count"] += read_count
            result[group_key][umi]["mapq_sum"] += mapq * read_count
            result[group_key][umi]["mapq_n"] += read_count


def main():
    args = parse_args()

    # Initialize UMI clusterer
    clusterer = UMIClusterer(cluster_method=args.method)

    # Data structure:
    # result[(chrom, start, end, orientation, strand)][umi] = {
    #     "count": int,
    #     "mapq_sum": int,
    #     "mapq_n": int
    # }
    result = defaultdict(dict)

    # Read input file
    with open(args.input) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")

            group_key = tuple(cols[0:5])
            umi = cols[5].encode("ascii")

            read_count = int(cols[6])
            mapq = float(cols[7])

            if umi not in result[group_key]:
                result[group_key][umi] = {
                    "count": 0,
                    "mapq_sum": 0,
                    "mapq_n": 0
                }

            result[group_key][umi]["count"] += read_count
            result[group_key][umi]["mapq_sum"] += mapq * read_count
            result[group_key][umi]["mapq_n"] += read_count

    ## open the file with R2 rescued reads if specified

    
    # Optionally load rescue R2
    if args.rescueR2=="TRUE":
        library=args.input.replace(".reads_per_UMI_per_IS.bed", "")
        rescuedR2=library+"_R2rescued.reads_per_UMI_per_IS.bed"
        print("RescueR2 was set to TRUE, reading "+rescuedR2)
        
        if not Path(rescuedR2).is_file():
            raise FileNotFoundError(f"Input file not found: {rescuedR2}")

        load_file_into_result(rescuedR2, result)


    collapsed_rows = []

    # Process each genomic position
    for group_key, umi_dict in result.items():
        chrom, start, end, orientation, strand = group_key

        # umi-tools requires {UMI: count}
        umi_counts = {
            umi: d["count"]
            for umi, d in umi_dict.items()
        }

        clusters = clusterer(
            umi_counts,
            threshold=args.threshold
        )

        for cluster in clusters:
            # Sort UMIs by abundance
            cluster_sorted = sorted(
                cluster,
                key=lambda u: umi_counts[u],
                reverse=True
            )

            rep_umi = cluster_sorted[0].decode()

            counts = [umi_dict[u]["count"] for u in cluster_sorted]

            umi_mapqs = [
                umi_dict[u]["mapq_sum"] / umi_dict[u]["mapq_n"]
                for u in cluster_sorted
            ]

            cluster_mapq = (
                sum(umi_dict[u]["mapq_sum"] for u in cluster_sorted) /
                sum(umi_dict[u]["mapq_n"] for u in cluster_sorted)
            )

            collapsed_rows.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "orientation": orientation,
                "strand": strand,
                "UMI": rep_umi,
                "count": sum(counts),
                "avg_qual": round(cluster_mapq, 2),
                "UMIs": ",".join(u.decode() for u in cluster_sorted),
                "UMIs_count": ",".join(map(str, counts)),
                "nUMIs": len(cluster_sorted)
                
            })

    # Write output TSV
    df = pd.DataFrame(collapsed_rows)
    df.to_csv(
        args.output,
        sep="\t",
        index=False
    )


if __name__ == "__main__":
    main()

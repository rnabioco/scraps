import os
import gzip
import argparse
from collections import defaultdict
from multiprocessing import Pool

import pysam
from umi_tools.network import UMIClusterer

""" Count UMI-deduplicated reads per (cell, polyA site) from a featureCounts
assigned BAM (with CB/UB/XT/XS tags), replacing `umi_tools count`.

Grouping is per (cell-barcode, featureCounts XT site tag). Within each group the
distinct UMIs are collapsed with umi_tools' own directional network algorithm
(UMIClusterer, threshold=1) so the output matches
`umi_tools count --per-gene --per-cell` (default method=directional) bit-for-bit.
STAR already 1MM-directional error-corrects UMIs upstream; the directional pass
here operates per-site (not per-gene), which is the resolution scraps requires.

Output is a gzipped TSV with columns `gene\tcell\tcount`, sorted by (gene, cell),
matching the umi_tools count table consumed by scraps_to_matrix (readr::read_tsv).
"""


def _iter_groups(bam_path, contig):
    """Yield ((cell, xt), {umi_bytes: read_count}) for one contig.

    Reads are filtered to featureCounts-assigned records (XS == 'Assigned') with
    valid cell and UMI tags (STAR writes 'CB:Z:-' / 'UB:Z:-' for undetermined
    barcodes/UMIs; these are dropped, mirroring the `grep -v 'CB:Z:-|UB:Z:-'`
    prefilter in the original rules).
    """
    groups = defaultdict(lambda: defaultdict(int))
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(contig=contig):
            if read.is_unmapped:
                continue
            # count each fragment once: featureCounts -p tags both mates, but
            # umi_tools count collapses the pair (keyed on read1). Skip the
            # second mate and any secondary/supplementary alignments. Single-end
            # reads are not flagged read2, so they are retained.
            if read.is_read2 or read.is_secondary or read.is_supplementary:
                continue
            try:
                if read.get_tag("XS") != "Assigned":
                    continue
                cb = read.get_tag("CB")
                ub = read.get_tag("UB")
                xt = read.get_tag("XT")
            except KeyError:
                continue
            if cb == "-" or ub == "-":
                continue
            groups[(cb, xt)][ub.encode()] += 1
    return groups


def _count_contig(args):
    """Worker: collapse UMIs per group for one contig -> list of (xt, cb, count)."""
    bam_path, contig = args
    clusterer = UMIClusterer(cluster_method="directional")
    rows = []
    for (cb, xt), umi_counts in _iter_groups(bam_path, contig).items():
        clusters = clusterer(umi_counts, threshold=1)
        rows.append((xt, cb, len(clusters)))
    return rows


def count_sites(bam_path, out_path, threads):
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        contigs = list(bam.references)

    work = [(bam_path, c) for c in contigs]
    if threads > 1 and len(work) > 1:
        with Pool(min(threads, len(work))) as pool:
            results = pool.map(_count_contig, work)
    else:
        results = [_count_contig(w) for w in work]

    rows = [r for sub in results for r in sub]
    # umi_tools count sorts the output table by gene, then cell
    rows.sort(key=lambda x: (x[0], x[1]))

    with gzip.open(out_path, "wt") as out:
        out.write("gene\tcell\tcount\n")
        for gene, cell, count in rows:
            out.write(f"{gene}\t{cell}\t{count}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Directional UMI count per (cell, site) from an assigned BAM"
    )
    parser.add_argument("-i", "--inbam",
                        help="featureCounts assigned BAM (CB/UB/XT/XS tags)",
                        required=True)
    parser.add_argument("-o", "--out",
                        help="output gzipped counts TSV (gene/cell/count)",
                        required=True)
    parser.add_argument("-t", "--threads",
                        help="worker processes (parallel over contigs)",
                        type=int, default=1)
    args = parser.parse_args()

    if not os.path.exists(args.inbam + ".bai") and not os.path.exists(
            args.inbam.rsplit(".", 1)[0] + ".bai"):
        pysam.index(args.inbam)

    count_sites(args.inbam, args.out, args.threads)


if __name__ == "__main__":
    main()

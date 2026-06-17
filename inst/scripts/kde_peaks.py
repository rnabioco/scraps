""" Call RT priming site summits via kernel density estimation.

Input is a strand-aware, single-base, UMI-deduplicated priming pileup bed
(chrom, start, end, count, strand). For each (chrom, strand) the per-base UMI
counts are smoothed with a Gaussian kernel (fixed bandwidth, in bp) to estimate
the RT priming density. Local maxima of the density above a threshold are
called as priming sites; maxima closer than --merge-dist are collapsed, and the
reported summit is the highest-count base within the merged region.

Output is a gzipped TSV: chrom, summit, strand, kde_score, umi_support
where summit is a 0-based single-base position.
"""

import argparse
import gzip
import sys

import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input',
                        help="stranded single-base priming bed (gzipped): "
                             "chrom start end count strand",
                        required=True)
    parser.add_argument('-o', '--output',
                        help="output peaks TSV (gzipped)",
                        required=True)
    parser.add_argument('--bandwidth', type=float, default=10.0,
                        help="Gaussian KDE bandwidth in bp (default 10)")
    parser.add_argument('--min-density', type=float, default=0.0,
                        help="minimum smoothed density to call a peak "
                             "(default 0, i.e. any local maximum)")
    parser.add_argument('--merge-dist', type=int, default=24,
                        help="collapse maxima within this distance in bp "
                             "(default 24)")
    return parser.parse_args()


def read_bed(path):
    """Read stranded single-base bed into {(chrom, strand): {pos: count}}."""
    opener = gzip.open if path.endswith('.gz') else open
    data = {}
    with opener(path, 'rt') as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            count = float(fields[3])
            strand = fields[4] if len(fields) > 4 else '+'
            data.setdefault((chrom, strand), {})[pos] = count
    return data


def gaussian_kernel(bandwidth):
    """Discrete Gaussian kernel truncated at +/- 4 sigma."""
    radius = max(1, int(np.ceil(4 * bandwidth)))
    x = np.arange(-radius, radius + 1)
    kernel = np.exp(-0.5 * (x / bandwidth) ** 2)
    kernel /= kernel.sum()
    return kernel, radius


def call_peaks_one(positions, counts, kernel, radius,
                   min_density, merge_dist):
    """Call peaks on one (chrom, strand) track.

    Returns list of (summit_pos, kde_score, umi_support).
    """
    pmin = positions.min() - radius
    span = positions.max() + radius - pmin + 1
    # raw per-base UMI signal on a contiguous grid
    raw = np.zeros(span, dtype=np.float64)
    raw[positions - pmin] = counts
    # smooth (density estimate) via convolution with Gaussian kernel
    density = np.convolve(raw, kernel, mode='same')

    # candidate local maxima: strictly greater than neighbors, above threshold
    n = density.size
    peaks = []
    for i in range(1, n - 1):
        d = density[i]
        if d <= min_density:
            continue
        if d >= density[i - 1] and d > density[i + 1]:
            peaks.append(i)
    if not peaks:
        return []

    # merge maxima within merge_dist into clusters
    merged = []
    cluster = [peaks[0]]
    for idx in peaks[1:]:
        if idx - cluster[-1] <= merge_dist:
            cluster.append(idx)
        else:
            merged.append(cluster)
            cluster = [idx]
    merged.append(cluster)

    results = []
    for cluster in merged:
        lo = cluster[0]
        hi = cluster[-1]
        # summit = highest raw UMI base within the merged window; tie -> max density
        window = range(lo, hi + 1)
        summit_idx = max(window, key=lambda j: (raw[j], density[j]))
        kde_score = float(density[summit_idx])
        umi_support = float(raw[lo:hi + 1].sum())
        results.append((summit_idx + pmin, kde_score, umi_support))
    return results


def main():
    args = parse_args()
    if args.bandwidth <= 0:
        sys.exit("--bandwidth must be > 0")

    kernel, radius = gaussian_kernel(args.bandwidth)
    data = read_bed(args.input)

    n_sites = 0
    with gzip.open(args.output, 'wt') as out:
        out.write("chrom\tsummit\tstrand\tkde_score\tumi_support\n")
        for (chrom, strand) in sorted(data.keys()):
            posmap = data[(chrom, strand)]
            positions = np.fromiter(posmap.keys(), dtype=np.int64)
            counts = np.fromiter(posmap.values(), dtype=np.float64)
            order = np.argsort(positions)
            positions = positions[order]
            counts = counts[order]
            for summit, score, support in call_peaks_one(
                    positions, counts, kernel, radius,
                    args.min_density, args.merge_dist):
                out.write("{}\t{}\t{}\t{:.6g}\t{:.0f}\n".format(
                    chrom, summit, strand, score, support))
                n_sites += 1

    sys.stderr.write(
        "kde_peaks: called {} RT priming sites across {} (chrom, strand) "
        "tracks\n".format(n_sites, len(data)))


if __name__ == '__main__':
    main()

""" Call RT priming site summits via kernel density estimation.

Input is a strand-aware, single-base, UMI-deduplicated priming pileup bed
(chrom, start, end, count, strand). For each (chrom, strand) the per-base UMI
counts are smoothed with a Gaussian kernel (fixed bandwidth, in bp) to estimate
the RT priming density. Local maxima of the density above a threshold are
called as priming sites; maxima closer than --merge-dist are collapsed, and the
reported summit is the highest-count base within the merged region.

Because discovery operates on pseudobulk (UMI-deduplicated per cell, then summed
across all cells and, optionally, across pooled samples), a called peak is only
retained if its summed UMI support clears a *depth-relative* floor:

    umi_support >= max(--min-umi, --min-umi-frac * T)

where T is the total UMI count across the whole input. This scales the bar with
sequencing depth and with sample pooling, so a fixed absolute count is not
required. The smoothed-density floor (--min-density) is a secondary shape gate;
by default it is derived from the bandwidth (~2 UMIs concentrated within one
bandwidth) so that isolated single-UMI bumps cannot form a called maximum.

Density is evaluated only at occupied bases (sparse), which is equivalent to the
dense convolution at those positions but avoids allocating full-chromosome
arrays.

Output is a gzipped TSV: chrom, summit, strand, kde_score, umi_support
where summit is a 0-based single-base position.
"""

import argparse
import gzip
import sys

import numpy as np

# sentinel: --min-density not set by the user -> derive from bandwidth
_DENSITY_AUTO = -1.0
# density floor expressed as "this many UMIs concentrated within one bandwidth"
_DENSITY_UMIS_PER_BW = 2.0


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
    parser.add_argument('--min-density', type=float, default=_DENSITY_AUTO,
                        help="minimum smoothed density to call a peak. If unset "
                             "(or <0), derived from --bandwidth as "
                             "%g / (sqrt(2*pi) * bandwidth), i.e. ~%g UMIs "
                             "concentrated within one bandwidth."
                             % (_DENSITY_UMIS_PER_BW, _DENSITY_UMIS_PER_BW))
    parser.add_argument('--merge-dist', type=int, default=24,
                        help="collapse maxima within this distance in bp "
                             "(default 24)")
    parser.add_argument('--min-umi', type=int, default=10,
                        help="absolute minimum summed UMI support per peak "
                             "(default 10); safety floor for shallow datasets")
    parser.add_argument('--min-umi-frac', type=float, default=1e-6,
                        help="minimum UMI support as a fraction of total UMIs T; "
                             "the effective floor is max(--min-umi, "
                             "--min-umi-frac * T) (default 1e-6). Depth-relative, "
                             "so the bar scales with pooling.")
    return parser.parse_args()


def gaussian_kernel(bandwidth):
    """Discrete Gaussian kernel truncated at +/- 4 sigma."""
    radius = max(1, int(np.ceil(4 * bandwidth)))
    x = np.arange(-radius, radius + 1)
    kernel = np.exp(-0.5 * (x / bandwidth) ** 2)
    kernel /= kernel.sum()
    return kernel, radius


def read_bed(path):
    """Read stranded single-base bed into {(chrom, strand): {pos: count}}.

    Also returns T, the total UMI count across all bases/strands.
    """
    opener = gzip.open if path.endswith('.gz') else open
    data = {}
    total = 0.0
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
            total += count
    return data, total


def sparse_density(positions, counts, kernel, radius):
    """Smoothed density evaluated only at occupied bases.

    Equivalent to convolving the dense per-base signal with `kernel` and reading
    off the values at `positions`, but without allocating a full-span array.
    `positions` must be sorted ascending.
    """
    n = positions.size
    density = counts * kernel[radius]
    for delta in range(1, radius + 1):
        w = kernel[radius + delta]
        # find, for each base i, the base exactly `delta` bp to its left
        idx = np.searchsorted(positions, positions - delta)
        valid = (idx < n)
        # guard against out-of-range before comparing positions
        idx_clipped = np.clip(idx, 0, n - 1)
        valid &= (positions[idx_clipped] == positions - delta)
        li = np.nonzero(valid)[0]
        if li.size:
            j = idx[li]
            # symmetric: left base contributes to i, and i contributes to left
            density[li] += counts[j] * w
            density[j] += counts[li] * w
    return density


def call_peaks_one(positions, counts, kernel, radius,
                   min_density, merge_dist, support_floor):
    """Call peaks on one (chrom, strand) track.

    Returns list of (summit_pos, kde_score, umi_support). Peaks below
    `support_floor` (summed UMI support over the merged cluster) are dropped.
    `positions` must be sorted ascending.
    """
    density = sparse_density(positions, counts, kernel, radius)

    # local maxima among occupied bases: strictly greater than occupied
    # neighbors, above the density floor (strict both sides -> no plateau
    # double-calling)
    n = density.size
    left = np.concatenate(([np.inf], density[:-1]))
    right = np.concatenate((density[1:], [np.inf]))
    ismax = (density > min_density) & (density > left) & (density > right)
    peak_idx = np.nonzero(ismax)[0]
    if peak_idx.size == 0:
        return []

    peak_pos = positions[peak_idx]
    # merge maxima whose genomic positions are within merge_dist into clusters
    breaks = np.nonzero(np.diff(peak_pos) > merge_dist)[0]
    starts = np.concatenate(([0], breaks + 1))
    ends = np.concatenate((breaks, [peak_idx.size - 1]))

    # cumulative sum over occupied counts for O(1) cluster support
    csum = np.concatenate(([0.0], np.cumsum(counts)))

    results = []
    for a, b in zip(starts, ends):
        lo_pos = peak_pos[a]
        hi_pos = peak_pos[b]
        lo = np.searchsorted(positions, lo_pos)
        hi = np.searchsorted(positions, hi_pos, side='right')
        umi_support = float(csum[hi] - csum[lo])
        if umi_support < support_floor:
            continue
        # summit = highest raw UMI base within the merged window; tie -> density
        window = range(lo, hi)
        summit_j = max(window, key=lambda j: (counts[j], density[j]))
        kde_score = float(density[summit_j])
        results.append((int(positions[summit_j]), kde_score, umi_support))
    return results


def main():
    args = parse_args()
    if args.bandwidth <= 0:
        sys.exit("--bandwidth must be > 0")

    kernel, radius = gaussian_kernel(args.bandwidth)

    # resolve density floor: derive from bandwidth unless explicitly set
    if args.min_density < 0:
        min_density = _DENSITY_UMIS_PER_BW * kernel[radius]
    else:
        min_density = args.min_density

    data, total_umi = read_bed(args.input)
    support_floor = max(float(args.min_umi), args.min_umi_frac * total_umi)

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
                    min_density, args.merge_dist, support_floor):
                out.write("{}\t{}\t{}\t{:.6g}\t{:.0f}\n".format(
                    chrom, summit, strand, score, support))
                n_sites += 1

    sys.stderr.write(
        "kde_peaks: called {} RT priming sites across {} (chrom, strand) "
        "tracks; support floor = max({}, {:g} * {:.0f}) = {:.0f} UMIs; "
        "min_density = {:.6g}\n".format(
            n_sites, len(data), args.min_umi, args.min_umi_frac,
            total_umi, support_floor, min_density))


if __name__ == '__main__':
    main()

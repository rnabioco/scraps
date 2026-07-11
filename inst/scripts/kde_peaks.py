""" Call RT priming site summits from a strand-aware priming pileup.

Input is a strand-aware, single-base, UMI-deduplicated priming pileup bed
(chrom, start, end, count, strand). For each (chrom, strand) the per-base UMI
counts are smoothed with a Gaussian kernel (fixed bandwidth, in bp) to estimate
the RT priming density. Local maxima of the density above a threshold are
called as priming sites; maxima closer than --merge-dist are collapsed into one
cluster. Gaussian smoothing is used for cluster DETECTION and umi_support.

Summit placement within each cluster is controlled by --summit-method:

  downstream (default)  3'-end read pileups have a hard boundary at the cleavage
                        site and tail UPSTREAM (transcript orientation), so the
                        smoothed density maximum sits ~1-2 bp upstream of the
                        true cleavage site (in the tail's centre of mass). This
                        method instead takes the most-downstream (transcript-3')
                        base whose depth >= --summit-frac * cluster peak depth,
                        searched over the full occupied read pileup, snapping the
                        summit onto the cleavage boundary. Validated against
                        polyAdb: exact-hit rate rises from ~25% (kde) to ~65% and
                        the systematic upstream bias (~ -1.4 bp) is removed at
                        --summit-frac 0.75. The kde_score column is the raw UMI
                        depth at the summit under this method.

  kde                   legacy behavior: summit = highest (raw count, then
                        smoothed density) base within the local-maxima span;
                        kde_score is the smoothed density at the summit.

Because discovery operates on pseudobulk (UMI-deduplicated per cell, then summed
across all cells and, optionally, across pooled samples), a called peak is only
retained if its summed UMI support clears a floor:

    umi_support >= max(--min-umi, --min-umi-frac * T)

where T is the total UMI count across the whole input. The primary control is
the ABSOLUTE floor --min-umi. The depth-relative term --min-umi-frac * T is
DISABLED by default (--min-umi-frac 0): it scales linearly with sequencing depth
and pooling, so on deep or pooled data it imposes a punishing absolute bar (e.g.
1e-6 * 1e8 = 100 UMIs/peak) that silently drops usable mid-abundance peaks. Set
a small positive --min-umi-frac only if you deliberately want the bar to scale
with depth. The smoothed-density floor (--min-density) is the shape gate; by
default it is derived from the bandwidth (~2 UMIs concentrated within one
bandwidth) so that isolated single-UMI bumps cannot form a called maximum.

Density is evaluated only at occupied bases (sparse), which is equivalent to the
dense convolution at those positions but avoids allocating full-chromosome
arrays.

Output is a gzipped TSV: chrom, summit, strand, kde_score, umi_support
where summit is a 0-based single-base position and kde_score is the raw UMI
depth at the summit (downstream method) or the smoothed density (kde method).
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
                             "(default 10); this is the PRIMARY support floor")
    parser.add_argument('--min-umi-frac', type=float, default=0.0,
                        help="optional depth-relative UMI support floor as a "
                             "fraction of total UMIs T; effective floor is "
                             "max(--min-umi, --min-umi-frac * T). DISABLED by "
                             "default (0): the linear-in-depth term over-filters "
                             "usable peaks on deep/pooled data. Set a small "
                             "positive value only to make the bar scale with "
                             "depth.")
    parser.add_argument('--summit-method', choices=('downstream', 'kde'),
                        default='downstream',
                        help="how to place the summit within a called cluster. "
                             "'downstream' (default): most-downstream "
                             "(transcript-3') base with depth >= --summit-frac * "
                             "cluster peak depth, searched over the full occupied "
                             "read pileup -- snaps to the cleavage boundary "
                             "(score = raw UMI depth). 'kde' (legacy): max(count, "
                             "smoothed density) within the local-maxima span "
                             "(score = smoothed density), which is biased ~1-2 bp "
                             "upstream into the read tail.")
    parser.add_argument('--summit-frac', type=float, default=0.75,
                        help="for --summit-method downstream: depth threshold as "
                             "a fraction of the cluster peak depth (default 0.75, "
                             "empirically minimizes summit-vs-polyAdb offset).")
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
                   min_density, merge_dist, support_floor,
                   strand, summit_method, summit_frac):
    """Call peaks on one (chrom, strand) track.

    Returns list of (summit_pos, score, umi_support). Peaks below
    `support_floor` (summed UMI support over the merged cluster) are dropped.
    `positions` must be sorted ascending.

    Cluster DETECTION is unchanged (Gaussian density local maxima merged within
    merge_dist). Only the SUMMIT within each cluster depends on summit_method:

      'kde'        summit = max(raw count, density) within the local-maxima span
                   (legacy behavior). `score` = smoothed density at the summit.

      'downstream' 3'-end read pileups have a hard boundary at the cleavage site
                   and tail UPSTREAM (transcript orientation); the Gaussian
                   summit is pulled ~1-2 bp upstream into the tail's centre of
                   mass. Instead, expand the summit search to the full occupied
                   read cluster (bases connected within merge_dist around the
                   called maxima) and take the MOST-DOWNSTREAM (transcript-3')
                   base whose depth >= summit_frac * cluster-peak-depth. This
                   snaps the summit onto the cleavage boundary. `score` = raw UMI
                   depth at the summit.
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
        # local-maxima span: defines the called cluster + its UMI support (this
        # is unchanged, so the set of called sites and their support match the
        # legacy cluster-detection behavior)
        lo = int(np.searchsorted(positions, peak_pos[a]))
        hi = int(np.searchsorted(positions, peak_pos[b], side='right'))
        umi_support = float(csum[hi] - csum[lo])
        if umi_support < support_floor:
            continue

        if summit_method == 'downstream':
            # expand the SUMMIT SEARCH (not the cluster/support) to the full
            # occupied read pileup connected within merge_dist
            slo, shi = lo, hi
            while slo > 0 and (positions[slo] - positions[slo - 1]) <= merge_dist:
                slo -= 1
            while shi < n and (positions[shi] - positions[shi - 1]) <= merge_dist:
                shi += 1
            peak_depth = counts[slo:shi].max()
            thr = summit_frac * peak_depth
            # eligible bases at/above the fractional threshold; take the most
            # downstream in transcript orientation ('-' -> lowest genomic coord)
            elig = [j for j in range(slo, shi) if counts[j] >= thr]
            summit_j = min(elig) if strand == '-' else max(elig)
            score = float(counts[summit_j])
        else:  # 'kde' (legacy)
            window = range(lo, hi)
            summit_j = max(window, key=lambda j: (counts[j], density[j]))
            score = float(density[summit_j])

        results.append((int(positions[summit_j]), score, umi_support))
    return results


def main():
    args = parse_args()
    if args.bandwidth <= 0:
        sys.exit("--bandwidth must be > 0")
    if not (0.0 < args.summit_frac <= 1.0):
        sys.exit("--summit-frac must be in (0, 1]")

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
                    min_density, args.merge_dist, support_floor,
                    strand, args.summit_method, args.summit_frac):
                out.write("{}\t{}\t{}\t{:.6g}\t{:.0f}\n".format(
                    chrom, summit, strand, score, support))
                n_sites += 1

    sys.stderr.write(
        "kde_peaks: called {} RT priming sites across {} (chrom, strand) "
        "tracks; summit_method = {}{}; support floor = max({}, {:g} * {:.0f}) "
        "= {:.0f} UMIs; min_density = {:.6g}\n".format(
            n_sites, len(data), args.summit_method,
            (" (frac %g)" % args.summit_frac
             if args.summit_method == 'downstream' else ""),
            args.min_umi, args.min_umi_frac,
            total_umi, support_floor, min_density))


if __name__ == '__main__':
    main()

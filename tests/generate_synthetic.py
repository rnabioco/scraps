"""Generate synthetic *aligned* BAMs that mimic STARsolo output for scraps, for
each alignment mode (R1, R2, paired). Reads carry CB/UB tags (STAR already
1MM-directional error-corrects UMIs upstream) and are positioned in and around
real polyAdb SAF windows, including deliberate edge cases:

  - PCR duplicates: same (CB, UB) repeated at one position
  - distinct molecules: different UBs at one position
  - edit-distance-1 error UMIs at lower abundance (directional network stress)
  - the same (CB, UB) at multiple positions within one window (window- vs
    position-first probe)

The output BAMs are meant to be run through the REAL featureCounts (as the
assign_sites_* rules do) to obtain XT/XS tags, then through both the old
(umi_tools) and new (count_sites.py / pileup_sites.py) code paths for a
bit-exactness diff. See tests/validate_synthetic.sh.

Usage:
  python3 tests/generate_synthetic.py --saf ref/polyadb32.hg38.saf.gz \
      --outdir <dir> [--seed 42]
"""
import os
import gzip
import random
import argparse

import pysam

BASES = "ACGT"
READLEN = 80


def load_windows(saf, chrom_filter="chr1", n=30, min_gap=2000):
    windows = []
    with gzip.open(saf, "rt") as fh:
        next(fh)
        for line in fh:
            gid, chrom, start, end, strand = line.rstrip("\n").split("\t")
            if chrom != chrom_filter:
                continue
            windows.append((gid, chrom, int(start), int(end), strand))
    chosen, last = [], -10**9
    for w in windows:
        if w[2] - last > min_gap:
            chosen.append(w)
            last = w[3]
        if len(chosen) >= n:
            break
    return chosen


def rand_umi(rng, k=10):
    return "".join(rng.choice(BASES) for _ in range(k))


def mutate(rng, u):
    i = rng.randrange(len(u))
    return u[:i] + rng.choice([b for b in BASES if b != u[i]]) + u[i + 1:]


def header():
    return {"HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 250_000_000}]}


def se_read(qname, pos0, strand, cb, ub):
    a = pysam.AlignedSegment()
    a.query_name = qname
    a.query_sequence = "A" * READLEN
    a.flag = 16 if strand == "-" else 0
    a.reference_id = 0
    a.reference_start = pos0
    a.mapping_quality = 255
    a.cigartuples = [(0, READLEN)]
    a.query_qualities = pysam.qualitystring_to_array("I" * READLEN)
    a.set_tag("CB", cb, "Z")
    a.set_tag("UB", ub, "Z")
    a.set_tag("NH", 1, "i")
    a.set_tag("HI", 1, "i")
    return a


def pe_pair(qname, pos1, pos2, strand, cb, ub):
    """First-in-pair (read1) at pos1, mate (read2) at pos2. Both carry tags.
    For the paired bed path only read1 (flag 0x42) is piled at its 5' end."""
    r1 = se_read(qname, pos1, strand, cb, ub)
    r2 = se_read(qname, pos2, "+" if strand == "-" else "-", cb, ub)
    # read1 flags: paired(1) + proper(2) + mate-reverse depending + read1(64)
    r1.flag = 0x1 | 0x2 | 0x40 | (0x10 if strand == "-" else 0) | (0x20 if strand == "+" else 0)
    r2.flag = 0x1 | 0x2 | 0x80 | (0x10 if strand == "+" else 0) | (0x20 if strand == "-" else 0)
    r1.next_reference_id = 0
    r1.next_reference_start = pos2
    r2.next_reference_id = 0
    r2.next_reference_start = pos1
    return r1, r2


def build(mode, chosen, cells, rng):
    """Return sorted list of AlignedSegments for a given mode."""
    records = []
    rid = 0
    for gid, chrom, start, end, strand in chosen:
        for cell in cells:
            n_true = rng.randint(1, 5)
            true_umis = [rand_umi(rng) for _ in range(n_true)]
            for tu in true_umis:
                for _ in range(rng.randint(1, 5)):  # PCR dups
                    p = _pos_in_window(start, end, strand, rng)
                    records += _emit(mode, f"r{rid}", p, strand, cell, tu, rng)
                    rid += 1
                if rng.random() < 0.5:  # edit-distance-1 error copies
                    eu = mutate(rng, tu)
                    for _ in range(rng.randint(1, 2)):
                        p = _pos_in_window(start, end, strand, rng)
                        records += _emit(mode, f"r{rid}", p, strand, cell, eu, rng)
                        rid += 1
            # same (CB,UB) at two distinct positions in the window
            us = rand_umi(rng)
            for p in (_pos_in_window(start, end, strand, rng),
                      _pos_in_window(start, end, strand, rng)):
                records += _emit(mode, f"r{rid}", p, strand, cell, us, rng)
                rid += 1
    records.sort(key=lambda a: a.reference_start)
    return records


def _pos_in_window(start, end, strand, rng):
    if strand == "-":
        return rng.randint(start, end - 1)
    return max(rng.randint(start, end - 1) - (READLEN - 1), 0)


def _emit(mode, qname, p, strand, cell, ub, rng):
    if mode == "paired":
        # mate placed a short insert away; strand as PAS strand for read1
        mate = max(p + rng.randint(100, 400), 0)
        r1, r2 = pe_pair(qname, p, mate, strand, cell, ub)
        return [r1, r2]
    return [se_read(qname, p, strand, cell, ub)]


def main():
    ap = argparse.ArgumentParser(description="synthetic aligned BAMs for scraps validation")
    ap.add_argument("--saf", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.seed)
    chosen = load_windows(args.saf)
    cells = ["".join(rng.choice(BASES) for _ in range(16)) for _ in range(6)]

    for mode in ("R1", "R2", "paired"):
        recs = build(mode, chosen, cells, rng)
        out = os.path.join(args.outdir, f"synth_{mode}_aligned.bam")
        with pysam.AlignmentFile(out, "wb", header=header()) as bam:
            for a in recs:
                bam.write(a)
        pysam.index(out)
        print(f"{mode}: wrote {len(recs)} records -> {out}")


if __name__ == "__main__":
    main()

"""Generate synthetic *aligned* BAMs that mimic STARsolo (and long-read
aligner) output for scraps. Reads carry CB/UB tags (STAR already
1MM-directional error-corrects UMIs upstream) and are positioned in and around
real polyAdb SAF windows, including deliberate edge cases:

  - PCR duplicates: same (CB, UB) repeated at one position
  - distinct molecules: different UBs at one position
  - edit-distance-1 error UMIs at lower abundance (directional network stress)
  - the same (CB, UB) at multiple positions within one window (window- vs
    position-first probe)

Fixtures written to --outdir:
  synth_{R1,R2,paired}_aligned.bam           STARsolo-like short reads per mode
  synth_bam_{same,reverse,paired1}_aligned.bam
                                             pre-aligned long multi-block reads for
                                             the `bam` mode (sense / antisense /
                                             paired-end with mate 1 as priming mate)
  synth_shift_R2_aligned.bam + shift_windows.saf
                                             priming-offset expectations (R2)
  synth_shift_{bam,bamrev}_aligned.bam       the same for bam mode under -L
  synth_clip_aligned.bam                     filter_bam.py --end 3/5 expectations

The BAMs are meant to be run through the REAL featureCounts (as the
assign_sites_* rules do) to obtain XT/XS tags, then through both the old
(umi_tools) and new (count_sites.py / pileup_sites.py) code paths: count tables
must be bit-exact, beds directional. See tests/validate_synthetic.sh.

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


def se_read(qname, pos0, strand, cb, ub, soft5=0, splice=0):
    """Single read at alignment start pos0. Optional 5' soft-clip (soft5 bp) and
    an internal splice (splice bp N) exercise the real chemistry
    (--clip5pNbases) and spliced alignments that shifted coordinates on real
    data. The soft-clip is placed at the read's 5' end in transcript orientation
    (leading S on forward reads, trailing S on reverse)."""
    a = pysam.AlignedSegment()
    a.query_name = qname
    qlen = READLEN + soft5
    a.query_sequence = "A" * qlen
    a.flag = 16 if strand == "-" else 0
    a.reference_id = 0
    a.reference_start = pos0
    a.mapping_quality = 255
    if splice:
        half = READLEN // 2
        core = [(0, half), (3, splice), (0, READLEN - half)]  # M N M
    else:
        core = [(0, READLEN)]
    if soft5:
        if strand == "-":       # reverse read: 5' end is downstream -> trailing S
            a.cigartuples = core + [(4, soft5)]
        else:                   # forward read: 5' end upstream -> leading S
            a.cigartuples = [(4, soft5)] + core
    else:
        a.cigartuples = core
    a.query_qualities = pysam.qualitystring_to_array("I" * qlen)
    a.set_tag("CB", cb, "Z")
    a.set_tag("UB", ub, "Z")
    a.set_tag("NH", 1, "i")
    a.set_tag("HI", 1, "i")
    return a


def pe_pair(qname, pos1, pos2, strand, cb, ub, soft5=0, proper=True):
    """First-in-pair (read1) at pos1, mate (read2) at pos2. Both carry tags.
    For the paired bed path only read1 (flag 0x42) is piled at its 5' end.
    proper=False marks a non-proper pair (mate on a different contig position /
    unmatched), which the pileup must discard (matching umi_tools --paired)."""
    r1 = se_read(qname, pos1, strand, cb, ub, soft5=soft5)
    r2 = se_read(qname, pos2, "+" if strand == "-" else "-", cb, ub)
    proper_bit = 0x2 if proper else 0
    # read1 flags: paired(1) + proper? + mate-reverse depending + read1(64)
    r1.flag = 0x1 | proper_bit | 0x40 | (0x10 if strand == "-" else 0) | (0x20 if strand == "+" else 0)
    r2.flag = 0x1 | proper_bit | 0x80 | (0x10 if strand == "+" else 0) | (0x20 if strand == "-" else 0)
    r1.next_reference_id = 0
    r1.next_reference_start = pos2
    r2.next_reference_id = 0
    r2.next_reference_start = pos1
    return r1, r2


def build(mode, chosen, cells, rng):
    """Return sorted list of AlignedSegments for a given mode.

    Exercises the cases that real data exposed but the original clean synthetic
    data missed: variable 5' soft-clips (--clip5pNbases chemistry), spliced
    reads, PCR duplicates that share (cell, UMI, priming position) but differ in
    mate/fragmentation position, and (paired) non-proper / unpaired fragments
    that must be discarded.
    """
    records = []
    rid = 0
    for gid, chrom, start, end, strand in chosen:
        for cell in cells:
            n_true = rng.randint(1, 5)
            true_umis = [rand_umi(rng) for _ in range(n_true)]
            for tu in true_umis:
                for _ in range(rng.randint(1, 5)):  # PCR dups (varied soft-clip)
                    p = _pos_in_window(start, end, strand, rng)
                    records += _emit(mode, f"r{rid}", p, strand, cell, tu, rng)
                    rid += 1
                if rng.random() < 0.5:  # edit-distance-1 error copies
                    eu = mutate(rng, tu)
                    for _ in range(rng.randint(1, 2)):
                        p = _pos_in_window(start, end, strand, rng)
                        records += _emit(mode, f"r{rid}", p, strand, cell, eu, rng)
                        rid += 1
            # same (CB,UB) at two distinct priming positions in the window
            us = rand_umi(rng)
            for p in (_pos_in_window(start, end, strand, rng),
                      _pos_in_window(start, end, strand, rng)):
                records += _emit(mode, f"r{rid}", p, strand, cell, us, rng)
                rid += 1
            # PCR-duplicate stress: same (CB,UB) at the SAME priming position but
            # DIFFERENT mate (fragmentation) positions. Fragment dedup keeps both;
            # position dedup collapses to one. (paired only)
            if mode == "paired":
                uf = rand_umi(rng)
                pp = _pos_in_window(start, end, strand, rng)
                for _ in range(3):
                    mate = max(pp + rng.randint(100, 600), 0)
                    r1, r2 = pe_pair(f"r{rid}", pp, mate, strand, cell, uf,
                                     soft5=rng.choice([0, 48, 57]))
                    records += [r1, r2]
                    rid += 1
                # non-proper / unpaired fragment: must be discarded by the pileup
                r1, r2 = pe_pair(f"r{rid}", pp, max(pp + 5000, 0), strand, cell,
                                 rand_umi(rng), proper=False)
                records += [r1, r2]
                rid += 1
    records.sort(key=lambda a: a.reference_start)
    return records


def _pos_in_window(start, end, strand, rng):
    if strand == "-":
        return rng.randint(start, end - 1)
    return max(rng.randint(start, end - 1) - (READLEN - 1), 0)


def _emit(mode, qname, p, strand, cell, ub, rng):
    soft5 = rng.choice([0, 0, 48, 57])   # most reads clipped, some varied
    splice = rng.choice([0, 0, 0, 1000])  # occasional spliced alignment
    if mode == "paired":
        # mate placed a short insert away; strand as PAS strand for read1
        mate = max(p + rng.randint(100, 400), 0)
        r1, r2 = pe_pair(qname, p, mate, strand, cell, ub, soft5=soft5)
        return [r1, r2]
    # R1 and R2 are single-end alignments (starsolo_R1/_R2 take one FASTQ)
    return [se_read(qname, p, strand, cell, ub, soft5=soft5, splice=splice)]


def long_read(qname, ref_start0, strand, cb, ub, blocks, introns, mapq=60):
    """Multi-block spliced alignment (M N M N ... M), as a long-read aligner
    emits for ONT/PacBio cDNA. `blocks` are the M lengths, `introns` the N
    lengths between them. Reference span = sum(blocks) + sum(introns)."""
    a = pysam.AlignedSegment()
    a.query_name = qname
    qlen = sum(blocks)
    a.query_sequence = "G" * qlen
    a.flag = 16 if strand == "-" else 0
    a.reference_id = 0
    a.reference_start = ref_start0
    a.mapping_quality = mapq
    cig = []
    for i, b in enumerate(blocks):
        if i:
            cig.append((3, introns[i - 1]))
        cig.append((0, b))
    a.cigartuples = cig
    a.query_qualities = pysam.qualitystring_to_array("I" * qlen)
    a.set_tag("CB", cb, "Z")
    a.set_tag("UB", ub, "Z")
    a.set_tag("NH", 1, "i")
    a.set_tag("HI", 1, "i")
    return a


def _lr_geometry(rng, nmin=6, nmax=14):
    """Random block/intron layout exceeding featureCounts' non-`-L` M-block cap
    (10) often enough to matter. Span stays < ~12 kb."""
    n = rng.randint(nmin, nmax)
    blocks = [rng.randint(60, 220) for _ in range(n)]
    introns = [rng.randint(100, 800) for _ in range(n - 1)]
    return blocks, introns


def _place(t_1based, right_end, span):
    """0-based reference_start for a read whose counted end must sit at the
    1-based base t: counted end at the right end -> reference_end == t;
    at the left end -> reference_start == t-1."""
    return (t_1based - span) if right_end else (t_1based - 1)


def _counted_end_is_right(read_strand, orientation):
    """The counted end is the read 3' end for `same` and the 5' end for
    `reverse`. On the reference, a forward read's 3' end is its right end and
    a reverse read's 3' end is its left end (and vice versa for 5')."""
    return (read_strand == "+") == (orientation == "same")


def build_bam(chosen, cells, rng, orientation, mate=0):
    """Pre-aligned long-read BAM fixture for the scraps `bam` mode.

    orientation 'same'   : read on the PAS strand, its 3' end in the window
    orientation 'reverse': read on the opposite strand, its 5' end in the window
    mate 0               : single-end records
    mate 1|2             : proper pairs; the priming mate (that number) has the
                           geometry above, the other mate sits a short insert
                           away and must be dropped; a few non-proper pairs and
                           secondary/supplementary copies must be dropped too.
    PCR duplicates share (cell, UMI, counted end) but vary the OTHER end (read
    length / block layout), so position-keyed dedup collapses them while
    mapping-coordinate dedup (umi_tools) may not -> directional bed check.
    """
    records = []
    rid = 0
    windows = [w for w in chosen if w[2] > 60000]
    for gid, chrom, start, end, pas_strand in windows:
        read_strand = pas_strand if orientation == "same" else ("-" if pas_strand == "+" else "+")
        right = _counted_end_is_right(read_strand, orientation)
        for cell in cells:
            for tu in [rand_umi(rng) for _ in range(rng.randint(1, 4))]:
                t = rng.randint(start, end)  # counted end, 1-based, inside the window
                for _ in range(rng.randint(1, 4)):  # PCR duplicates
                    records += _emit_lr(f"r{rid}", t, right, read_strand, cell, tu, rng, mate)
                    rid += 1
                if rng.random() < 0.4:  # edit-distance-1 error copy
                    records += _emit_lr(f"r{rid}", t, right, read_strand, cell, mutate(rng, tu), rng, mate)
                    rid += 1
            if mate:
                # non-proper pair: must be discarded
                recs = _emit_lr(f"r{rid}", rng.randint(start, end), right, read_strand,
                                cell, rand_umi(rng), rng, mate, proper=False)
                records += recs
                rid += 1
        # secondary / supplementary copies of a real read, far away: must be dropped
        for bit in (0x100, 0x800):
            blocks, introns = _lr_geometry(rng)
            a = long_read(f"r{rid}", start + 50000, read_strand, cells[0], rand_umi(rng), blocks, introns)
            a.flag |= bit
            records.append(a)
            rid += 1
    records.sort(key=lambda a: (a.reference_start, a.query_name))
    return records


def _emit_lr(qname, t, right, read_strand, cell, ub, rng, mate, proper=True):
    blocks, introns = _lr_geometry(rng)
    span = sum(blocks) + sum(introns)
    a = long_read(qname, _place(t, right, span), read_strand, cell, ub, blocks, introns)
    if not mate:
        return [a]
    # paired: the priming mate is `mate`; the other mate is a 150 bp read on the
    # opposite strand, placed beyond the non-counted end
    other_strand = "-" if read_strand == "+" else "+"
    if right:
        other_pos = max(a.reference_start - rng.randint(200, 600) - 150, 0)
    else:
        other_pos = a.reference_start + span + rng.randint(200, 600)
    b = long_read(qname, other_pos, other_strand, cell, ub, [150], [])
    proper_bit = 0x2 if proper else 0
    me_bit, other_bit = (0x40, 0x80) if mate == 1 else (0x80, 0x40)
    a.flag = 0x1 | proper_bit | me_bit | (0x10 if read_strand == "-" else 0) | (0x20 if other_strand == "-" else 0)
    b.flag = 0x1 | proper_bit | other_bit | (0x10 if other_strand == "-" else 0) | (0x20 if read_strand == "-" else 0)
    a.next_reference_id = 0
    a.next_reference_start = b.reference_start
    b.next_reference_id = 0
    b.next_reference_start = a.reference_start
    return [a, b]


def clipped_read(qname, pos0, strand, cb, ub, lead="", trail=""):
    """80 bp aligned body of G's with explicit leading/trailing soft-clip
    sequences (reference orientation), for the filter_bam.py tests."""
    a = pysam.AlignedSegment()
    a.query_name = qname
    seq = lead + "G" * READLEN + trail
    a.query_sequence = seq
    a.flag = 16 if strand == "-" else 0
    a.reference_id = 0
    a.reference_start = pos0
    a.mapping_quality = 255
    cig = []
    if lead:
        cig.append((4, len(lead)))
    cig.append((0, READLEN))
    if trail:
        cig.append((4, len(trail)))
    a.cigartuples = cig
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    a.set_tag("CB", cb, "Z")
    a.set_tag("UB", ub, "Z")
    return a


def build_softclip(rng):
    """Reads whose qname encodes the expected filter_bam.py verdict for each
    counted end: '<keep3|drop3>-<keep5|drop5>_<case>_<rep>'. Thresholds (script
    defaults n=5, f=0.8): clip > 5 bp with >= 80% of the expected base, where
    end 3 = forward trailing clip of A / reverse leading clip of T, and
    end 5 = forward leading clip of T / reverse trailing clip of A."""
    A10, T10, MIX = "A" * 10, "T" * 10, "ACGTACGTAC"
    cases = [
        # (strand, lead, trail, verdict3, verdict5)
        ("+", "", A10, "keep3", "drop5"),      # poly(A) at fwd 3' end
        ("+", T10, "", "drop3", "keep5"),      # poly(T) at fwd 5' end
        ("+", T10, A10, "keep3", "keep5"),     # both
        ("+", "", "AAAA", "drop3", "drop5"),   # too short
        ("+", "", MIX, "drop3", "drop5"),      # low A fraction
        ("+", "", "", "drop3", "drop5"),       # no clip
        ("+", A10, "", "drop3", "drop5"),      # A's at the 5' end: wrong base for end 5
        ("-", T10, "", "keep3", "drop5"),      # rev read: 3' end is the leading clip, T
        ("-", "", A10, "drop3", "keep5"),      # rev read: 5' end is the trailing clip, A
        ("-", T10, A10, "keep3", "keep5"),
        ("-", "TTT", "", "drop3", "drop5"),
        ("-", MIX, "", "drop3", "drop5"),
        ("-", "", "", "drop3", "drop5"),
        ("-", "", T10, "drop3", "drop5"),      # T's at rev 5' (trailing): wrong base
    ]
    records = []
    for n, (strand, lead, trail, v3, v5) in enumerate(cases):
        for rep in range(3):
            pos0 = 1_000_000 + 1000 * (n * 3 + rep)
            records.append(clipped_read(f"{v3}-{v5}_{n}_{rep}", pos0, strand,
                                        "ACGT" * 4, rand_umi(rng), lead, trail))
    records.sort(key=lambda a: a.reference_start)
    return records


def build_shift_lr(chosen, cells, rng, orientation, shift=5):
    """Long multi-block reads at exact offsets from each SAF window for the
    `bam` mode (`-L`): same prefixes/expectations as build_shift (lrshiftin_,
    lrinside_, lredge_, lrfar_) plus `lrbody_`, whose body overlaps the window
    but whose counted end is 200 bp past it -> must be Unassigned with and
    without shift (proves --read2pos is honoured under -L). orientation
    'same' counts the 3' end (-s 1 --read2pos 3, shift downstream); 'reverse'
    counts the 5' end (-s 2 --read2pos 5, shift upstream)."""
    records = []
    rid = 0
    for gid, chrom, start, end, pas_strand in [w for w in chosen if w[2] > 60000]:
        read_strand = pas_strand if orientation == "same" else ("-" if pas_strand == "+" else "+")
        right = _counted_end_is_right(read_strand, orientation)
        if pas_strand == "+":
            first, last, step = start, end, 1
        else:
            first, last, step = end, start, -1
        cases = {
            "lrshiftin": first - step * shift,
            "lrinside": first + step * 2,
            "lredge": last,
            "lrfar": first - step * (shift + 1),
            "lrbody": last + step * 200,
        }
        for prefix, t in cases.items():
            for cell in cells[:2]:
                blocks, introns = _lr_geometry(rng, nmin=11, nmax=14)
                span = sum(blocks) + sum(introns)
                records.append(long_read(f"{prefix}_{rid}", _place(t, right, span),
                                         read_strand, cell, rand_umi(rng), blocks, introns))
                rid += 1
    records.sort(key=lambda a: a.reference_start)
    return records


def build_shift(chosen, cells, rng, shift=5):
    """R2-style single-end reads (no soft clips, like Ultima FASTQ cDNA reads)
    placed at exact offsets from each SAF window so that a featureCounts run
    with `--read2pos 3 --readShiftType downstream --readShiftSize <shift>` has
    a known expected assignment per read. The qname prefix encodes the case:

      shiftin_ : 3' base exactly `shift` bp UPSTREAM (transcript 5') of the
                 window's first transcript base -> Unassigned without shift,
                 Assigned with shift
      inside_  : 3' base inside the window, >= shift bp from its downstream
                 edge -> Assigned both ways
      edge_    : 3' base within the last `shift` bases of the window ->
                 Assigned without shift, Unassigned with shift (the true
                 cleavage site lies past the window)
      far_     : shift+1 bp upstream of the window -> Unassigned both ways

    Transcript direction is + = increasing coordinate, - = decreasing. The
    read's 3' base is reference_end-1 (forward) or reference_start (reverse),
    so a 'downstream' (read-3') shift is +N on forward and -N on reverse reads.
    Returns (records, windows) so the caller can also write a SAF restricted
    to these windows (the full polyAdb SAF has neighbouring sites that would
    otherwise make 'Unassigned' expectations unreliable).
    """
    def rec(prefix, rid, t_3prime_1based, strand, cell, ub):
        # forward: 1-based 3' base = reference_start0 + READLEN
        # reverse: 1-based 3' base = reference_start0 + 1
        pos0 = (t_3prime_1based - READLEN) if strand == "+" else (t_3prime_1based - 1)
        return se_read(f"{prefix}_{rid}", max(pos0, 0), strand, cell, ub)

    records = []
    rid = 0
    for gid, chrom, start, end, strand in chosen:  # start/end 1-based inclusive
        if strand == "+":
            first, last, step = start, end, 1
        else:
            first, last, step = end, start, -1
        cases = {
            "shiftin": first - step * shift,
            "inside": first + step * 2,
            "edge": last,
            "far": first - step * (shift + 1),
        }
        for prefix, t in cases.items():
            for cell in cells[:2]:
                records.append(rec(prefix, rid, t, strand, cell, rand_umi(rng)))
                rid += 1
    records.sort(key=lambda a: a.reference_start)
    return records, chosen


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

    # priming-offset fixture (Ultima FASTQ: site = R2 3' end + 5)
    recs, windows = build_shift(chosen, cells, rng, shift=5)
    out = os.path.join(args.outdir, "synth_shift_R2_aligned.bam")
    with pysam.AlignmentFile(out, "wb", header=header()) as bam:
        for a in recs:
            bam.write(a)
    pysam.index(out)
    saf_out = os.path.join(args.outdir, "shift_windows.saf")
    with open(saf_out, "w") as fh:
        fh.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for gid, chrom, start, end, strand in windows:
            fh.write(f"{gid}\t{chrom}\t{start}\t{end}\t{strand}\n")
    print(f"shift: wrote {len(recs)} records -> {out} (SAF subset: {saf_out})")

    # bam mode fixtures: pre-aligned long reads, both orientations, one paired
    for label, orientation, mate in (("bam_same", "same", 0),
                                     ("bam_reverse", "reverse", 0),
                                     ("bam_paired1", "reverse", 1)):
        recs = build_bam(chosen, cells, rng, orientation, mate=mate)
        out = os.path.join(args.outdir, f"synth_{label}_aligned.bam")
        with pysam.AlignmentFile(out, "wb", header=header()) as bam:
            for a in recs:
                bam.write(a)
        pysam.index(out)
        print(f"{label}: wrote {len(recs)} records -> {out}")

    # long-read shift fixtures (featureCounts -L + --read2pos + --readShiftType)
    for label, orientation in (("bam", "same"), ("bamrev", "reverse")):
        recs = build_shift_lr(chosen, cells, rng, orientation, shift=5)
        out = os.path.join(args.outdir, f"synth_shift_{label}_aligned.bam")
        with pysam.AlignmentFile(out, "wb", header=header()) as bam:
            for a in recs:
                bam.write(a)
        pysam.index(out)
        print(f"shift_{label}: wrote {len(recs)} records -> {out}")

    # soft-clip filter fixture (filter_bam.py --end 3 / --end 5)
    recs = build_softclip(rng)
    out = os.path.join(args.outdir, "synth_clip_aligned.bam")
    with pysam.AlignmentFile(out, "wb", header=header()) as bam:
        for a in recs:
            bam.write(a)
    pysam.index(out)
    print(f"clip: wrote {len(recs)} records -> {out}")


if __name__ == "__main__":
    main()

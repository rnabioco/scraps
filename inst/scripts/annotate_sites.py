""" Annotate de novo RT priming sites.

Takes KDE-called RT priming site summits (from kde_peaks.py) and classifies
each into one of three categories:

    known_pas               summit falls within the strand-aware, asymmetric
                            transcript-oriented match window (--match-up /
                            --match-down) of a validated polyAdb cleavage point
    likely_internal_priming downstream genomic A-content/A-run suggests priming
                            off a genomic poly(A) stretch (artifact)
    potential_novel_pas     no internal-priming signature and a canonical
                            poly(A) signal (PAS) motif is present upstream

Known-site matching is point-based against the polyAdb single-base cleavage
position, which is read from field 5 of the ';'- (hg38) or '_'- (mm10)
delimited SAF GeneID (gene;genbank;id;chrom;pos;strand;class). A SAF that does
not follow this encoding will not annotate any known sites. See docs/discovery.md.

IMPORTANT: validated poly(A) sites (polyAdb) remain the recommended reference
for downstream quantitative analyses. potential_novel_pas sites are candidates
that require orthogonal validation; likely_internal_priming sites should be
excluded, not quantified.

Outputs:
  --out-tsv  full annotation table (gzipped)
  --out-bed  BED6 of sites (name = status), gzipped
  --out-saf  de novo SAF reusable as POLYA_SITES, GeneID encoded as
             gene;genbank;id;chrom;summit;strand;class  (class = status)
"""

import argparse
import gzip
import sys

import pysam

# canonical poly(A) signal + common single-base variants (Beaudoing et al.)
PAS_MOTIFS = [
    "AATAAA", "ATTAAA", "AGTAAA", "TATAAA", "CATAAA", "GATAAA",
    "AATATA", "AATACA", "AATAGA", "AAAAAG", "ACTAAA", "AAGAAA", "AATGAA",
]

COMPLEMENT = str.maketrans("ACGTacgtN", "TGCAtgcaN")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--peaks', required=True,
                        help="KDE peaks TSV (gzipped) from kde_peaks.py")
    parser.add_argument('--polya', required=True,
                        help="polyAdb SAF reference (gzipped)")
    parser.add_argument('--fasta', required=True,
                        help="genome FASTA (with .fai index)")
    parser.add_argument('--out-tsv', required=True)
    parser.add_argument('--out-bed', required=True)
    parser.add_argument('--out-saf', required=True)
    parser.add_argument('--match-up', type=int, default=10,
                        help="max distance (bp) a summit may lie UPSTREAM "
                             "(transcript orientation) of a polyAdb cleavage "
                             "point to call a match (default 10)")
    parser.add_argument('--match-down', type=int, default=5,
                        help="max distance (bp) a summit may lie DOWNSTREAM "
                             "(transcript orientation) of a polyAdb cleavage "
                             "point to call a match (default 5)")
    parser.add_argument('--saf-slop-up', type=int, default=10)
    parser.add_argument('--saf-slop-down', type=int, default=5)
    parser.add_argument('--pas-start', type=int, default=-40,
                        help="upstream PAS window start (bp rel. to summit)")
    parser.add_argument('--pas-end', type=int, default=0,
                        help="upstream PAS window end (bp rel. to summit)")
    parser.add_argument('--a-start', type=int, default=1,
                        help="downstream A-content window start (bp)")
    parser.add_argument('--a-end', type=int, default=20,
                        help="downstream A-content window end (bp)")
    parser.add_argument('--a-content-max', type=float, default=0.7,
                        help="downstream A fraction to flag internal priming")
    parser.add_argument('--a-run-min', type=int, default=8,
                        help="contiguous downstream A-run length to flag "
                             "internal priming (default 8)")
    return parser.parse_args()


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def _split_geneid(geneid, sep, chrom, strand):
    """Split a GeneID into its 7 fields, robust to delimiters inside fields.

    The GeneID encodes gene;genbank;id;chrom;pos;strand;class. The ';' (human)
    delimiter never appears inside a field, so a plain split works. The '_'
    (mouse) delimiter DOES occur inside fields - RefSeq IDs like 'NR_152944'
    and scaffold chromosome names like 'chrUn_GL456...' contain underscores -
    so a naive split is unreliable.

    The trailing fields chrom/pos/strand/class are recovered unambiguously by
    anchoring on the SAF-column ``chrom``/``strand`` and reading from the right
    (class = last token, strand = next, pos = next, chrom = the SAF-column
    value). The leading gene/genbank/id metadata cannot be split reliably when
    it contains the delimiter, so it is returned as a single best-effort group
    in the gene field (genbank/id set to 'NA'); these fields are passthrough
    metadata only and do not affect matching, which uses pos/chrom/strand.

    Returns [gene, genbank, id, chrom, pos, strand, class] or None.
    """
    parts = geneid.split(sep)
    if sep == ';':
        return parts if len(parts) == 7 else None

    # exact 7-field with consistent chrom/strand: trust it
    if len(parts) == 7 and parts[3] == chrom and parts[5] == strand:
        return parts

    # '_' delimiter with embedded underscores: anchor from the right
    if len(parts) < 7:
        return None
    klass = parts[-1]
    strand_tok = parts[-2]
    pos_tok = parts[-3]
    if strand_tok != strand:
        return None
    chrom_tokens = chrom.split(sep)
    nct = len(chrom_tokens)
    chrom_start = len(parts) - 3 - nct
    if chrom_start < 1:
        return None
    if parts[chrom_start:chrom_start + nct] != chrom_tokens:
        return None
    head = parts[:chrom_start]  # gene/genbank/id group (delimiter-ambiguous)
    gene = sep.join(head) if head else 'NA'
    return [gene, 'NA', 'NA', chrom, pos_tok, strand, klass]


def load_polya(path):
    """Load polyAdb sites: {(chrom, strand): sorted [(pos1based, geneid)]}.

    The polyAdb SAF GeneID must use the scraps 7-field encoding
    (gene;genbank;id;chrom;pos;strand;class for human, '_'-delimited for mouse),
    where the single-base cleavage position is field 5 (0-based index 4) and the
    site class is field 7 (index 6). Rows whose GeneID is not 7-field are skipped
    with a warning; if no rows parse, the SAF format is almost certainly wrong
    and the user is alerted (every discovered site would otherwise be labeled
    novel). Helper converters in inst/scripts/polyadb/ produce conforming SAFs.
    """
    opener = gzip.open if path.endswith('.gz') else open
    sites = {}
    n_total = 0
    n_skipped = 0
    bad_examples = []
    with opener(path, 'rt') as fh:
        header = fh.readline()  # GeneID Chr Start End Strand
        for line in fh:
            if not line.strip():
                continue
            n_total += 1
            cols = line.rstrip('\n').split('\t')
            geneid = cols[0]
            chrom = cols[1]
            strand = cols[4]
            sep = ';' if ';' in geneid else '_'
            parts = _split_geneid(geneid, sep, chrom, strand)
            if parts is None:
                n_skipped += 1
                if len(bad_examples) < 3:
                    bad_examples.append(geneid)
                continue
            try:
                # coerce float/scientific-notation encodings (e.g. '7.7e+07'
                # from a float round-trip) back to the integer cleavage position
                pos = int(round(float(parts[4])))
            except (ValueError, IndexError):
                n_skipped += 1
                if len(bad_examples) < 3:
                    bad_examples.append(geneid)
                continue
            meta = {
                'gene': parts[0],
                'genbank': parts[1],
                'id': parts[2],
                'class': parts[6],
            }
            sites.setdefault((chrom, strand), []).append((pos, meta))
    for key in sites:
        sites[key].sort(key=lambda t: t[0])

    n_parsed = n_total - n_skipped
    if n_skipped:
        sys.stderr.write(
            "annotate_sites: WARNING - skipped {}/{} polyAdb rows with a "
            "non-7-field GeneID. Examples: {}\n".format(
                n_skipped, n_total, "; ".join(bad_examples)))
    if n_total and n_parsed == 0:
        sys.stderr.write(
            "annotate_sites: ERROR - no polyAdb rows parsed. The SAF GeneID "
            "must use the scraps 7-field encoding "
            "(gene;genbank;id;chrom;pos;strand;class). See "
            "inst/scripts/polyadb/ for converters. All sites will be novel.\n")
    return sites


def nearest_polya(sites, chrom, strand, pos1, match_up, match_down):
    """Match a summit to the nearest polyAdb cleavage point, strand-aware.

    The match window is asymmetric in transcript orientation: a summit matches
    if it lies within ``match_up`` bp upstream or ``match_down`` bp downstream
    of a polyAdb cleavage point.

    Returns (signed_dist, meta) where signed_dist is the transcript-oriented
    offset of the summit relative to the polyAdb point (negative = upstream,
    positive = downstream), or (None, None) if no site is in the window.
    """
    arr = sites.get((chrom, strand))
    if not arr:
        return None, None
    # binary search for first polyAdb point >= summit
    lo, hi = 0, len(arr)
    while lo < hi:
        mid = (lo + hi) // 2
        if arr[mid][0] < pos1:
            lo = mid + 1
        else:
            hi = mid

    # scan outward far enough to cover the widest side of the window
    reach = max(match_up, match_down)
    best_signed = None
    best_abs = None
    best_meta = None
    j = lo - 1
    while j >= 0 and (pos1 - arr[j][0]) <= reach:
        _consider = arr[j]
        signed = _signed_offset(strand, pos1, _consider[0])
        if -match_up <= signed <= match_down:
            if best_abs is None or abs(signed) < best_abs:
                best_abs, best_signed, best_meta = abs(signed), signed, _consider[1]
        j -= 1
    j = lo
    while j < len(arr) and (arr[j][0] - pos1) <= reach:
        _consider = arr[j]
        signed = _signed_offset(strand, pos1, _consider[0])
        if -match_up <= signed <= match_down:
            if best_abs is None or abs(signed) < best_abs:
                best_abs, best_signed, best_meta = abs(signed), signed, _consider[1]
        j += 1

    if best_meta is not None:
        return best_signed, best_meta
    return None, None


def _signed_offset(strand, summit_pos, polya_pos):
    """Transcript-oriented offset: negative=upstream, positive=downstream."""
    if strand == '-':
        return polya_pos - summit_pos
    return summit_pos - polya_pos


def fetch_seq(fasta, chrom, start0, end0):
    """Fetch [start0, end0) clamped to contig bounds, uppercased."""
    try:
        length = fasta.get_reference_length(chrom)
    except KeyError:
        return ""
    s = max(0, start0)
    e = min(length, end0)
    if e <= s:
        return ""
    return fasta.fetch(chrom, s, e).upper()


def max_a_run(seq):
    best = run = 0
    for base in seq:
        if base == 'A':
            run += 1
            best = max(best, run)
        else:
            run = 0
    return best


def seq_context(fasta, chrom, summit0, strand, args):
    """Return (pas_motif, pas_dist, a_content, a_run) in transcript orientation.

    summit0 is 0-based. Windows are expressed relative to the summit on the
    sense (transcript) strand: negative = upstream (5'), positive = downstream.
    """
    # downstream A-content window (sense strand, toward the poly(A) tail)
    if strand == '+':
        a_s, a_e = summit0 + args.a_start, summit0 + args.a_end + 1
        a_seq = fetch_seq(fasta, chrom, a_s, a_e)
    else:
        a_s, a_e = summit0 - args.a_end, summit0 - args.a_start + 1
        a_seq = revcomp(fetch_seq(fasta, chrom, a_s, a_e))

    a_content = (a_seq.count('A') / len(a_seq)) if a_seq else 0.0
    a_run = max_a_run(a_seq)

    # upstream PAS motif window (sense strand)
    if strand == '+':
        p_s, p_e = summit0 + args.pas_start, summit0 + args.pas_end + 1
        p_seq = fetch_seq(fasta, chrom, p_s, p_e)
    else:
        p_s, p_e = summit0 - args.pas_end, summit0 - args.pas_start + 1
        p_seq = revcomp(fetch_seq(fasta, chrom, p_s, p_e))

    pas_motif = "NA"
    pas_dist = "NA"
    best_off = None
    for motif in PAS_MOTIFS:
        idx = p_seq.find(motif)
        if idx >= 0:
            # distance upstream of summit (bp); window end aligns to summit
            off = len(p_seq) - idx
            if best_off is None or off < best_off:
                best_off = off
                pas_motif = motif
                pas_dist = off
    return pas_motif, pas_dist, a_content, a_run


def classify(matched, a_content, a_run, pas_motif, args):
    if matched:
        return "known_pas"
    internal = (a_content >= args.a_content_max) or (a_run >= args.a_run_min)
    if internal:
        return "likely_internal_priming"
    if pas_motif != "NA":
        return "potential_novel_pas"
    # no PAS signal and no strong internal-priming signature: conservative call
    return "likely_internal_priming"


def main():
    args = parse_args()
    sites = load_polya(args.polya)
    fasta = pysam.FastaFile(args.fasta)

    opener = gzip.open if args.peaks.endswith('.gz') else open
    counts = {"known_pas": 0, "likely_internal_priming": 0,
              "potential_novel_pas": 0}

    with opener(args.peaks, 'rt') as fh, \
            gzip.open(args.out_tsv, 'wt') as tsv, \
            gzip.open(args.out_bed, 'wt') as bed, \
            gzip.open(args.out_saf, 'wt') as saf:
        tsv.write("chrom\tsummit\tstrand\tkde_score\tumi_support\tstatus\t"
                  "gene\tid\tclass\tpolya_dist\tpas_motif\tpas_dist\t"
                  "a_content\ta_run\n")
        saf.write("GeneID\tChr\tStart\tEnd\tStrand\n")

        header = fh.readline()
        for line in fh:
            if not line.strip():
                continue
            chrom, summit, strand, kde_score, umi_support = \
                line.rstrip('\n').split('\t')
            summit0 = int(summit)
            pos1 = summit0 + 1  # 1-based to compare with polyAdb position

            dist, meta = nearest_polya(sites, chrom, strand, pos1,
                                       args.match_up, args.match_down)
            matched = meta is not None
            pas_motif, pas_dist, a_content, a_run = seq_context(
                fasta, chrom, summit0, strand, args)
            status = classify(matched, a_content, a_run, pas_motif, args)
            counts[status] += 1

            gene = meta['gene'] if matched else "NA"
            gid = meta['id'] if matched else "NA"
            genbank = meta['genbank'] if matched else "NA"
            klass = status
            polya_dist = dist if matched else "NA"

            tsv.write("\t".join(str(x) for x in [
                chrom, summit0, strand, kde_score, umi_support, status,
                gene, gid, klass, polya_dist, pas_motif, pas_dist,
                "{:.4f}".format(a_content), a_run]) + "\n")

            # BED6 (0-based, single base), name = status
            bed.write("\t".join(str(x) for x in [
                chrom, summit0, summit0 + 1, status, umi_support, strand
            ]) + "\n")

            # SAF: 1-based inclusive window with strand-aware slop
            if strand == '+':
                start1 = pos1 - args.saf_slop_up
                end1 = pos1 + args.saf_slop_down
            else:
                start1 = pos1 - args.saf_slop_down
                end1 = pos1 + args.saf_slop_up
            start1 = max(1, start1)
            geneid = ";".join([gene, genbank, gid, chrom, str(pos1),
                               strand, klass])
            saf.write("\t".join(str(x) for x in [
                geneid, chrom, start1, end1, strand]) + "\n")

    sys.stderr.write(
        "annotate_sites: {known_pas} known_pas, "
        "{likely_internal_priming} likely_internal_priming, "
        "{potential_novel_pas} potential_novel_pas\n".format(**counts))


if __name__ == '__main__':
    main()

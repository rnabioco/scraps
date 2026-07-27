#!/usr/bin/env python3
""" Validate a scraps SAF file (and optionally compare to a reference).

Checks:
  1. FORMAT  - header, 7-field GeneID, GeneID chrom/pos/strand consistency with
               the SAF columns, strand-specific window, sort order.
  2. CONTENT - row-level comparison against a reference SAF at shared
               (chrom, encoded position, strand) positions (optional).

The GeneID delimiter (';' human / '_' mouse) is auto-detected per file, or set
with --delim. Inputs may be gzipped.

Usage:
  python3 validate_saf.py GENERATED.saf[.gz] [REFERENCE.saf[.gz]]
"""

import argparse
import gzip
from collections import defaultdict


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def detect_delim(gene_id, override=None):
    if override:
        return override
    return ";" if ";" in gene_id else "_"


def split_geneid(gene_id, delim, chrom, strand):
    """Robustly split a 7-field GeneID, tolerant of '_' inside fields.

    The '_' (mouse) delimiter occurs inside RefSeq IDs ('NR_152944') and some
    scaffold chromosome names, so trailing fields are recovered by anchoring on
    the SAF-column chrom/strand and reading from the right. Returns
    [gene, genbank, id, chrom, pos, strand, class] or None.
    """
    parts = gene_id.split(delim)
    if delim == ";":
        return parts if len(parts) == 7 else None
    if len(parts) == 7 and parts[3] == chrom and parts[5] == strand:
        return parts
    if len(parts) < 7:
        return None
    klass, strand_tok, pos_tok = parts[-1], parts[-2], parts[-3]
    if strand_tok != strand:
        return None
    chrom_tokens = chrom.split(delim)
    nct = len(chrom_tokens)
    chrom_start = len(parts) - 3 - nct
    if chrom_start < 1 or parts[chrom_start:chrom_start + nct] != chrom_tokens:
        return None
    head = parts[:chrom_start]
    gene = delim.join(head) if head else "NA"
    return [gene, "NA", "NA", chrom, pos_tok, strand, klass]


def load_saf(path, delim=None):
    records = []
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n")
        assert header == "GeneID\tChr\tStart\tEnd\tStrand", \
            "Bad header in {}: {!r}".format(path, header)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            assert len(f) == 5, "Expected 5 columns, got {}: {!r}".format(
                len(f), line)
            records.append({
                "gene_id": f[0],
                "chrom": f[1],
                "start": int(f[2]),
                "end": int(f[3]),
                "strand": f[4],
            })
    return records


def validate_format(records, label, delim=None):
    errors = []
    prev_chrom, prev_start = None, -1
    for i, r in enumerate(records):
        d = detect_delim(r["gene_id"], delim)
        g = split_geneid(r["gene_id"], d, r["chrom"], r["strand"])
        if g is None:
            errors.append(
                "[{}] record {}: GeneID not parseable as 7 fields: {!r}".format(
                    label, i, r["gene_id"]))
            continue
        _, _, _, gid_chrom, pos_str, gid_strand, _ = g
        if gid_chrom != r["chrom"]:
            errors.append("[{}] record {}: GeneID chrom={!r} != column {!r}".format(
                label, i, gid_chrom, r["chrom"]))
        if gid_strand != r["strand"]:
            errors.append("[{}] record {}: GeneID strand={!r} != column {!r}".format(
                label, i, gid_strand, r["strand"]))
        try:
            # tolerate float/scientific-notation encodings (e.g. '7.7e+07')
            pos = int(round(float(pos_str)))
        except ValueError:
            errors.append("[{}] record {}: GeneID pos not integer: {!r}".format(
                label, i, pos_str))
            continue
        window = r["end"] - r["start"]
        if window < 15:
            errors.append("[{}] record {}: window too small ({}bp)".format(
                label, i, window))
        if not (r["start"] <= pos <= r["end"]):
            errors.append(
                "[{}] record {}: GeneID pos={} outside window [{}, {}]".format(
                    label, i, pos, r["start"], r["end"]))
        if prev_chrom is not None:
            if r["chrom"] < prev_chrom:
                errors.append("[{}] record {}: chrom out of order: {!r} after {!r}".format(
                    label, i, r["chrom"], prev_chrom))
            elif r["chrom"] == prev_chrom and r["start"] < prev_start:
                errors.append("[{}] record {}: start out of order in {}: {} < {}".format(
                    label, i, r["chrom"], r["start"], prev_start))
        prev_chrom, prev_start = r["chrom"], r["start"]
    return errors


def index_positions(path, delim=None):
    d = defaultdict(set)
    with open_text(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            dd = detect_delim(f[0], delim)
            g = split_geneid(f[0], dd, f[1], f[4])
            if g is not None:
                d[(f[1], g[4], f[4])].add(f[0])  # chrom, enc_pos, strand
    return d


def compare_content(generated, reference, delim=None):
    gen = index_positions(generated, delim)
    ref = index_positions(reference, delim)
    gk, rk = set(gen), set(ref)
    shared = gk & rk
    exact = mismatch = 0
    examples = []
    for key in shared:
        if gen[key] == ref[key]:
            exact += 1
        else:
            mismatch += 1
            if len(examples) < 5:
                examples.append((key, sorted(gen[key]), sorted(ref[key])))
    return {
        "gen_rows": sum(len(v) for v in gen.values()),
        "ref_rows": sum(len(v) for v in ref.values()),
        "gen_positions": len(gk), "ref_positions": len(rk),
        "shared": len(shared), "only_gen": len(gk - rk), "only_ref": len(rk - gk),
        "exact": exact, "mismatch": mismatch, "examples": examples,
    }


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("generated")
    ap.add_argument("reference", nargs="?", default=None)
    ap.add_argument("--delim", default=None, help="force ';' or '_'")
    args = ap.parse_args()

    print("=== Format: {} ===".format(args.generated))
    records = load_saf(args.generated, args.delim)
    errors = validate_format(records, "generated", args.delim)
    if errors:
        print("FAIL - {} format error(s):".format(len(errors)))
        for e in errors[:20]:
            print("  " + e)
    else:
        print("PASS - {:,} records, all format checks clean".format(len(records)))

    if not args.reference:
        return

    print("\n=== Content: {} vs {} ===".format(args.generated, args.reference))
    r = compare_content(args.generated, args.reference, args.delim)
    print("Generated : {:>8,} rows ({:,} positions)".format(
        r["gen_rows"], r["gen_positions"]))
    print("Reference : {:>8,} rows ({:,} positions)".format(
        r["ref_rows"], r["ref_positions"]))
    print("Shared positions  : {:,}".format(r["shared"]))
    if r["shared"]:
        print("  Exact match     : {:,} ({:.2f}%)".format(
            r["exact"], 100 * r["exact"] / r["shared"]))
        print("  Mismatch        : {:,}".format(r["mismatch"]))
    print("Only in generated : {:,}".format(r["only_gen"]))
    print("Only in reference : {:,}".format(r["only_ref"]))
    if r["examples"]:
        print("\nMismatch examples:")
        for key, gg, rg in r["examples"]:
            print("  pos={}".format(key))
            for g in gg:
                print("    gen: {}".format(g))
            for g in rg:
                print("    ref: {}".format(g))


if __name__ == "__main__":
    main()

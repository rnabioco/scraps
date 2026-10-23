#!/usr/bin/env python3
""" Convert polyAdb 4 PAS tables to the scraps SAF format.

Handles both polyAdb 4 table types, for human and mouse:

  main : *.PAS.main.tsv  (rich annotation; GeneSymbol + PAS_type_RefSeq_label)
  max  : *.PAS.max.tsv   (per-site highest-expressed gene only; no PAS type)

The PAS_ID column encodes the site as ``chr:strand:pos`` (1-based, target build).
Columns are resolved by header name. The mode (main/max) is auto-detected from
the header, or set with --mode. The output is a 7-field-GeneID SAF:

    gene_symbol{D}refseq{D}ensembl{D}chrom{D}pos{D}strand{D}pas_type

polyAdb 4 tables do not provide RefSeq or Ensembl gene IDs, so those fields are
'NA'. The max tables also lack a PAS type, so pas_type is 'NA' there. PAS type
labels from v4 main tables are passed through VERBATIM and therefore differ from
the polyAdb 3.2 vocabulary (e.g. v4 "3'UTR"/"Intergenic" vs 3.2
"3'UTR(M)"/"intergenic").

liftOver is supported for future targets via --source-build/--target-build;
inputs are hg38/mm10 by default, so liftOver is skipped when the builds match.

Examples:
  python3 pasv4_to_saf.py --species human --pas HumanPas_v4/hg38.PAS.main.tsv \
      --source-build hg38 --target-build hg38 --out ref/polyadb4.hg38.main.saf.gz
  python3 pasv4_to_saf.py --species mouse --pas MousePas_v4/mm10.PAS.max.tsv \
      --source-build mm10 --target-build mm10 --out ref/polyadb4.mm10.max.saf.gz
"""

import argparse
import sys

import saf_common as sc


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", choices=["human", "mouse"], required=True)
    ap.add_argument("--pas", required=True, help="input *.PAS.{main,max}.tsv")
    ap.add_argument("--out", required=True,
                    help="output SAF (.gz writes gzipped)")
    ap.add_argument("--mode", choices=["main", "max"], default=None,
                    help="table type (default: auto-detect from header)")
    ap.add_argument("--source-build", required=True,
                    help="genome build of the input (e.g. hg38, mm10)")
    ap.add_argument("--target-build", required=True,
                    help="desired build; equal to --source-build skips liftOver")
    ap.add_argument("--liftover", default=None,
                    help="path to liftOver binary (default: `liftOver` on PATH, "
                         "else a platform-matched UCSC download)")
    ap.add_argument("--chain", default=None,
                    help="path to source->target chain (default: auto-download)")
    ap.add_argument("--delim", default=None, choices=sc.VALID_DELIMS,
                    help="override GeneID delimiter (default by --species)")
    ap.add_argument("--tmpdir", default=None,
                    help="scratch directory for liftOver files (default: a fresh "
                         "mkdtemp, removed afterwards; a user-supplied directory "
                         "is never removed)")
    ap.add_argument("--keep-tmp", action="store_true",
                    help="do not remove the mkdtemp scratch directory")
    return ap.parse_args()


def parse_pas_id(pas_id, where=""):
    """'chr1:+:940182' -> (chrom, strand, pos_int)."""
    try:
        chrom, strand, pos_str = pas_id.rsplit(":", 2)
    except ValueError:
        sys.exit("{}: PAS_ID not of the form chr:strand:pos: {!r}".format(where, pas_id))
    return chrom, sc.check_strand(strand, where), sc.parse_pos(pos_str, where)


def detect_mode(header_line):
    cols = header_line.rstrip("\n").split("\t")
    return "main" if "PAS_type_RefSeq_label" in cols else "max"


def parse_pas(path, mode):
    """Yield records from a v4 table using header-name column resolution."""
    with sc.open_text(path) as fh:
        header = fh.readline()
        if mode is None:
            mode = detect_mode(header)
        wanted = {"pas_id": "PAS_ID", "sym": "GeneSymbol"}
        if mode == "main":
            wanted["pas_type"] = "PAS_type_RefSeq_label"
        idx = sc.resolve_columns(header, wanted)
        for lineno, line in enumerate(fh, start=2):
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            where = "{} line {}".format(path, lineno)
            try:
                chrom, strand, pos = parse_pas_id(f[idx["pas_id"]], where)
                rec = {
                    "chrom": chrom,
                    "pos": pos,
                    "strand": strand,
                    "sym": f[idx["sym"]],
                    "pas_type": f[idx["pas_type"]] if mode == "main" else "NA",
                }
            except IndexError:
                sys.exit("{}: short/ragged line ({} columns)".format(where, len(f)))
            yield rec


def rows_no_liftover(records, delim):
    for r in records:
        start, end = sc.site_window(r["pos"], r["strand"])
        gene_id = sc.encode_geneid(
            r["sym"], "NA", "NA",
            r["chrom"], r["pos"], r["strand"], r["pas_type"], delim)
        yield (r["chrom"], start, end, r["strand"], r["pos"], gene_id)


def rows_with_liftover(records, delim, args, tmpdir):
    records = list(records)
    point_bed = str(tmpdir / "src_point.bed")
    with open(point_bed, "w") as fp:
        for i, r in enumerate(records):
            fp.write("{}\t{}\t{}\t{}\t0\t{}\n".format(
                r["chrom"], r["pos"] - 1, r["pos"], i, r["strand"]))

    liftover_bin = sc.ensure_liftover(tmpdir, args.liftover)
    chain = sc.ensure_chain(tmpdir, args.source_build, args.target_build,
                            args.chain)
    lifted_bed = str(tmpdir / "tgt_point.bed")
    unmapped = str(tmpdir / "unmapped_point.bed")
    n_unmapped = sc.run_liftover(liftover_bin, chain, point_bed, lifted_bed, unmapped)
    lifted = sc.load_lifted(lifted_bed)
    sys.stderr.write("  lifted {}/{} cleavage positions ({} unmapped)\n".format(
        len(lifted), len(records), n_unmapped))

    for i, r in enumerate(records):
        pt = lifted.get(str(i))
        if pt is None:
            continue
        tgt_chrom = pt[0]
        tgt_pos = pt[1] + 1
        tgt_strand = pt[3]
        start, end = sc.site_window(tgt_pos, tgt_strand)
        gene_id = sc.encode_geneid(
            r["sym"], "NA", "NA",
            tgt_chrom, tgt_pos, tgt_strand, r["pas_type"], delim)
        yield (tgt_chrom, start, end, tgt_strand, tgt_pos, gene_id)


def main():
    args = parse_args()
    delim = sc.check_delim(args.delim or sc.SPECIES_DELIM[args.species])
    need_tmp = args.source_build != args.target_build
    tmpdir, created = (sc.make_tmpdir(args.tmpdir, "pasv4_to_saf")
                       if need_tmp else (None, False))

    # determine mode (for the log) up front
    mode = args.mode
    if mode is None:
        with sc.open_text(args.pas) as fh:
            mode = detect_mode(fh.readline())
    sys.stderr.write("parsing {} (mode={}) ...\n".format(args.pas, mode))
    records = parse_pas(args.pas, mode)

    try:
        if not need_tmp:
            sys.stderr.write("  source build == target build; skipping liftOver\n")
            raw = rows_no_liftover(records, delim)
        else:
            sys.stderr.write("  lifting {} -> {}\n".format(
                args.source_build, args.target_build))
            raw = rows_with_liftover(records, delim, args, tmpdir)

        rows = sc.dedup_sort_rows(raw)
        sc.write_saf(rows, args.out)
        sys.stderr.write("wrote {} records -> {}\n".format(len(rows), args.out))
    finally:
        if need_tmp:
            sc.cleanup_tmpdir(tmpdir, created, args.keep_tmp)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
""" Convert polyAdb 3.2 PAS tables to the scraps SAF format.

Input is a polyAdb 3.2 ``*.PAS.txt`` table (e.g. human.PAS.txt, mouse.PAS.txt).
Columns are resolved by header name, so human and mouse files (which differ in
column count) both work.

The cleavage position is taken from the input and, when ``--source-build`` and
``--target-build`` differ, remapped with UCSC liftOver. The output is a
7-field-GeneID SAF (see saf_common.py):

    gene_symbol{D}refseq{D}ensembl{D}chrom{D}pos{D}strand{D}pas_type

with ';' (human) or '_' (mouse) delimiter.

Examples:
  # human, hg19 -> hg38 (auto-download liftOver + chain)
  python3 pas32_to_saf.py --species human --pas human.PAS.txt \
      --source-build hg19 --target-build hg38 \
      --out ref/polyadb32.hg38.saf.gz

  # mouse, already mm10 -> mm10 (no liftover)
  python3 pas32_to_saf.py --species mouse --pas mouse.PAS.txt \
      --source-build mm10 --target-build mm10 \
      --out ref/polyadb32.mm10.saf.gz
"""

import argparse
import shutil
import sys
from pathlib import Path

import saf_common as sc


COLUMNS = {
    "chrom": "Chromosome",
    "pos": "Position",
    "strand": "Strand",
    "ensembl": "Ensemble ID",      # polyAdb's spelling
    "refseq": "RefSeq Gene ID",
    "sym": "Gene Symbol",
    "pas_type": "PAS type",
}


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", choices=["human", "mouse"], required=True,
                    help="sets the GeneID delimiter (; human / _ mouse)")
    ap.add_argument("--pas", required=True, help="input *.PAS.txt (3.2)")
    ap.add_argument("--out", required=True,
                    help="output SAF (.gz writes gzipped)")
    ap.add_argument("--source-build", required=True,
                    help="genome build of the input (e.g. hg19, mm10)")
    ap.add_argument("--target-build", required=True,
                    help="desired genome build (e.g. hg38); equal to "
                         "--source-build skips liftOver")
    ap.add_argument("--liftover", default=None,
                    help="path to liftOver binary (default: auto-download)")
    ap.add_argument("--chain", default=None,
                    help="path to source->target chain (default: auto-download)")
    ap.add_argument("--delim", default=None,
                    help="override GeneID delimiter (default by --species)")
    ap.add_argument("--tmpdir", default="/tmp/pas32_to_saf")
    ap.add_argument("--keep-tmp", action="store_true")
    return ap.parse_args()


def parse_pas(path):
    """Yield records from a 3.2 PAS table using header-name column resolution."""
    with sc.open_text(path) as fh:
        idx = sc.resolve_columns(fh.readline(), COLUMNS)
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            yield {
                "chrom": f[idx["chrom"]],
                "pos": int(f[idx["pos"]]),
                "strand": f[idx["strand"]],
                "ensembl": f[idx["ensembl"]],
                "refseq": f[idx["refseq"]],
                "sym": f[idx["sym"]],
                "pas_type": f[idx["pas_type"]],
            }


def rows_no_liftover(records, delim):
    for r in records:
        start, end = sc.site_window(r["pos"], r["strand"])
        gene_id = sc.encode_geneid(
            r["sym"], r["refseq"], r["ensembl"],
            r["chrom"], r["pos"], r["strand"], r["pas_type"], delim)
        yield (r["chrom"], start, end, r["strand"], r["pos"], gene_id)


def rows_with_liftover(records, delim, args, tmpdir):
    """Lift the point cleavage position; rebuild the window in target coords."""
    records = list(records)
    point_bed = str(tmpdir / "src_point.bed")
    with open(point_bed, "w") as fp:
        for i, r in enumerate(records):
            # BED is 0-based half-open; single base [pos-1, pos)
            fp.write("{}\t{}\t{}\t{}\t0\t{}\n".format(
                r["chrom"], r["pos"] - 1, r["pos"], i, r["strand"]))

    liftover_bin = sc.ensure_liftover(tmpdir, args.liftover)
    chain = sc.ensure_chain(tmpdir, args.source_build, args.target_build,
                            args.chain)
    lifted_bed = str(tmpdir / "tgt_point.bed")
    unmapped = str(tmpdir / "unmapped_point.bed")
    sc.run_liftover(liftover_bin, chain, point_bed, lifted_bed, unmapped)
    lifted = sc.load_lifted(lifted_bed)
    sys.stderr.write("  lifted {}/{} cleavage positions\n".format(
        len(lifted), len(records)))

    for i, r in enumerate(records):
        pt = lifted.get(str(i))
        if pt is None:
            continue  # failed liftover
        tgt_chrom = pt[0]
        tgt_pos = pt[1] + 1  # 0-based BED start -> 1-based position
        tgt_strand = pt[3]
        start, end = sc.site_window(tgt_pos, tgt_strand)
        gene_id = sc.encode_geneid(
            r["sym"], r["refseq"], r["ensembl"],
            tgt_chrom, tgt_pos, tgt_strand, r["pas_type"], delim)
        yield (tgt_chrom, start, end, tgt_strand, tgt_pos, gene_id)


def main():
    args = parse_args()
    delim = args.delim or sc.SPECIES_DELIM[args.species]
    tmpdir = Path(args.tmpdir)
    tmpdir.mkdir(parents=True, exist_ok=True)

    sys.stderr.write("parsing {} ...\n".format(args.pas))
    records = parse_pas(args.pas)

    if args.source_build == args.target_build:
        sys.stderr.write("  source build == target build; skipping liftOver\n")
        raw = rows_no_liftover(records, delim)
    else:
        sys.stderr.write("  lifting {} -> {}\n".format(
            args.source_build, args.target_build))
        raw = rows_with_liftover(records, delim, args, tmpdir)

    rows = sc.dedup_sort_rows(raw)
    sc.write_saf(rows, args.out)
    sys.stderr.write("wrote {} records -> {}\n".format(len(rows), args.out))

    if not args.keep_tmp and tmpdir.exists():
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()

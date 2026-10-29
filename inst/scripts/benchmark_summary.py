"""Summarize Snakemake benchmark files per sample, broken down by rule.

Every scraps rule writes `{results}/benchmarks/<rule-stem>_<sample>_<mode>.tsv`
(discovery rules: `<rule>_<unit>.tsv`) via Snakemake's `benchmark:` directive.
Each file holds one row per benchmark repeat with the job's wall time `s`
(measured where the job ran, so SLURM queue wait is excluded), peak memory
`max_rss` (MB) and CPU time. Snakemake samples memory with Linux-only process
statistics, so memory columns are NA on macOS and for very short jobs. This script reads them all and prints, per sample
(or discovery unit), the wall time of each rule and the total.

Usage:
  python3 inst/scripts/benchmark_summary.py ../results/benchmarks [--tsv out.tsv]
"""

import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict

# rule stems in the order they run; longest-first matching keeps
# "assign_sites" from being read as "assign" + sample "sites_..."
RULE_STEMS = ["cutadapt", "starsolo", "filter", "assign_sites", "count",
              "stranded_bed", "pool_beds", "kde_peaks", "annotate_sites"]
UNIT_RULES = {"pool_beds", "kde_peaks", "annotate_sites"}
MODES = ("R1", "R2", "paired", "bam")


def parse_name(path):
    """benchmark file name -> (rule, group, mode); group is sample or unit."""
    base = os.path.basename(path)[:-len(".tsv")]
    for stem in sorted(RULE_STEMS, key=len, reverse=True):
        if base.startswith(stem + "_"):
            rest = base[len(stem) + 1:]
            if stem in UNIT_RULES:
                return stem, rest, ""
            m = re.match(r"^(.*)_({})$".format("|".join(MODES)), rest)
            if m:
                return stem, m.group(1), m.group(2)
            return stem, rest, ""
    return None


def read_benchmark(path):
    """Mean wall seconds and max max_rss over the file's rows (repeats)."""
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    secs = [float(r["s"]) for r in rows if r.get("s") not in (None, "", "NA")]
    rss = [float(r["max_rss"]) for r in rows if r.get("max_rss") not in (None, "", "NA", "-")]
    if not secs:
        return None
    return sum(secs) / len(secs), (max(rss) if rss else float("nan"))


def fmt(seconds):
    seconds = int(round(seconds))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return "{:d}:{:02d}:{:02d}".format(h, m, s)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("benchdir", help="{results}/benchmarks directory")
    ap.add_argument("--tsv", help="also write one row per benchmark file to this TSV")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.benchdir, "*.tsv")))
    if not files:
        sys.exit("no benchmark files in " + args.benchdir)

    rows = []
    skipped = []
    for path in files:
        parsed = parse_name(path)
        bench = read_benchmark(path)
        if parsed is None or bench is None:
            skipped.append(os.path.basename(path))
            continue
        rule, group, mode = parsed
        rows.append({"group": group, "rule": rule, "mode": mode,
                     "seconds": bench[0], "max_rss_mb": bench[1],
                     "file": os.path.basename(path)})

    by_group = defaultdict(list)
    for r in rows:
        by_group[r["group"]].append(r)
    order = {stem: i for i, stem in enumerate(RULE_STEMS)}
    for group in sorted(by_group):
        items = sorted(by_group[group], key=lambda r: (r["mode"], order[r["rule"]]))
        total = sum(r["seconds"] for r in items)
        print("{}  (total {})".format(group, fmt(total)))
        for r in items:
            label = r["rule"] + ("_" + r["mode"] if r["mode"] else "")
            pct = 100 * r["seconds"] / total if total else 0
            rss = r["max_rss_mb"]
            rss_str = "NA" if rss != rss else "{:.0f} MB".format(rss)  # NaN check
            print("  {:<22} {:>10} {:6.1f}%   max_rss {:>10}".format(
                label, fmt(r["seconds"]), pct, rss_str))
        print()
    if skipped:
        print("skipped (unrecognized name or empty): " + ", ".join(skipped))

    if args.tsv:
        with open(args.tsv, "w") as fh:
            fh.write("sample_or_unit\trule\tmode\tseconds\tmax_rss_mb\tfile\n")
            for r in rows:
                rss = r["max_rss_mb"]
                fh.write("{}\t{}\t{}\t{:.2f}\t{}\t{}\n".format(
                    r["group"], r["rule"], r["mode"], r["seconds"],
                    "NA" if rss != rss else "{:.1f}".format(rss), r["file"]))
        print("rows written to " + args.tsv)


if __name__ == "__main__":
    main()

""" Shared helpers for converting polyAdb releases to the scraps SAF format.

The scraps SAF format used as ``POLYA_SITES`` is a tab-separated table with a
header ``GeneID  Chr  Start  End  Strand`` where the GeneID encodes seven
fields:

    gene_symbol{D}refseq_gene_id{D}ensembl_id{D}chrom{D}pos{D}strand{D}pas_type

with delimiter ``{D}`` = ';' for human references and '_' for mouse references
(matching the bundled ref/polyadb32.{hg38,mm10}.saf.gz). ``pos`` is the 1-based
single-base cleavage position (field 5) and ``pas_type`` is the site class
(field 7). Missing values are encoded as 'NA'.

The SAF Start/End columns describe a strand-specific 15 bp window around the
cleavage position (transcript-oriented -10 / +5 slop):

    + strand: Start = pos - 10, End = pos + 5
    - strand: Start = pos - 5,  End = pos + 10

This module also provides an optional UCSC liftOver wrapper so converters can
remap coordinates to a target genome build.
"""

import gzip
import os
import stat
import subprocess
import sys
import urllib.request
from pathlib import Path

# transcript-oriented slop applied around the cleavage position
SLOP_UP = 10
SLOP_DOWN = 5

SAF_HEADER = "GeneID\tChr\tStart\tEnd\tStrand"

# delimiter convention per species (matches bundled references)
SPECIES_DELIM = {"human": ";", "mouse": "_"}


def normalise(val):
    """Return 'NA' for empty/na-like values, else the stripped value."""
    if val is None:
        return "NA"
    v = val.strip()
    if not v or v.lower() == "na":
        return "NA"
    return v


def encode_geneid(sym, refseq, ensembl, chrom, pos, strand, pas_type, delim):
    """Build the 7-field scraps GeneID string."""
    fields = [
        normalise(sym),
        normalise(refseq),
        normalise(ensembl),
        chrom,
        str(pos),
        strand,
        normalise(pas_type),
    ]
    return delim.join(fields)


def site_window(pos, strand):
    """Strand-specific SAF (Start, End) around a 1-based cleavage position.

    Mirrors the bundled polyAdb references: a 15 bp window with -10 / +5
    transcript-oriented slop. Returned coordinates use the same convention as
    those references (Start = pos - up_slop, End = pos + down_slop on '+').
    """
    if strand == "+":
        return pos - SLOP_UP, pos + SLOP_DOWN
    return pos - SLOP_DOWN, pos + SLOP_UP


def dedup_sort_rows(rows):
    """Deduplicate identical (chrom,pos,strand,GeneID) rows and sort.

    ``rows`` is an iterable of (chrom, start, end, strand, pos, gene_id).
    Returns a list of (chrom, start, end, strand, gene_id) sorted by
    (chrom lexicographic, start numeric).
    """
    seen = set()
    out = []
    for chrom, start, end, strand, pos, gene_id in rows:
        key = (chrom, pos, strand, gene_id)
        if key in seen:
            continue
        seen.add(key)
        out.append((chrom, start, end, strand, gene_id))
    out.sort(key=lambda r: (r[0], r[1]))
    return out


def write_saf(rows, path):
    """Write SAF rows. Output is gzipped if ``path`` ends with '.gz'.

    ``rows`` is an iterable of (chrom, start, end, strand, gene_id).
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as out:
        out.write(SAF_HEADER + "\n")
        for chrom, start, end, strand, gene_id in rows:
            out.write("{}\t{}\t{}\t{}\t{}\n".format(
                gene_id, chrom, start, end, strand))


def open_text(path):
    """Open a possibly-gzipped text file for reading."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def resolve_columns(header_line, wanted):
    """Map wanted column names to indices using a tab-separated header.

    ``wanted`` is a dict {key: header_name}. Returns {key: index}. Raises
    KeyError listing any header names not found (handles releases whose column
    order or count differs, e.g. human vs mouse polyAdb files).
    """
    cols = header_line.rstrip("\n").split("\t")
    index = {name: i for i, name in enumerate(cols)}
    out = {}
    missing = []
    for key, name in wanted.items():
        if name in index:
            out[key] = index[name]
        else:
            missing.append(name)
    if missing:
        raise KeyError(
            "columns not found in header: {} (available: {})".format(
                missing, cols))
    return out


# --------------------------------------------------------------------------
# liftOver support
# --------------------------------------------------------------------------

_LIFTOVER_BASE = "https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/liftOver"
_CHAIN_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/{frm}/liftOver/"
    "{frm}To{To}.over.chain.gz"
)


def _download(url, dest, label):
    if Path(dest).exists():
        sys.stderr.write("  {}: using existing {}\n".format(label, dest))
        return
    sys.stderr.write("  downloading {} from {}\n".format(label, url))
    urllib.request.urlretrieve(url, dest)


def ensure_liftover(tmpdir, liftover=None):
    """Return a path to a liftOver binary, downloading if not provided."""
    if liftover:
        return liftover
    lo = Path(tmpdir) / "liftOver"
    _download(_LIFTOVER_BASE, lo, "liftOver binary")
    lo.chmod(lo.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return str(lo)


def ensure_chain(tmpdir, source_build, target_build, chain=None):
    """Return a path to a decompressed chain file, downloading if needed."""
    if chain:
        if str(chain).endswith(".gz"):
            out = Path(tmpdir) / (Path(chain).stem)
            with gzip.open(chain, "rb") as src, open(out, "wb") as dst:
                dst.write(src.read())
            return str(out)
        return chain
    to_cap = target_build[0].upper() + target_build[1:]
    url = _CHAIN_URL.format(frm=source_build, To=to_cap)
    chain_gz = Path(tmpdir) / "{}To{}.over.chain.gz".format(source_build, to_cap)
    chain_txt = Path(tmpdir) / "{}To{}.over.chain".format(source_build, to_cap)
    _download(url, chain_gz, "chain file")
    if not chain_txt.exists():
        with gzip.open(chain_gz, "rb") as src, open(chain_txt, "wb") as dst:
            dst.write(src.read())
    return str(chain_txt)


def run_liftover(liftover_bin, chain, bed_in, bed_out, unmapped):
    """Run UCSC liftOver. liftOver exits 1 when some records are unmapped."""
    cmd = [liftover_bin, bed_in, chain, bed_out, unmapped]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode not in (0, 1):
        sys.exit("liftOver failed:\n{}".format(result.stderr))


def load_lifted(path):
    """Load a lifted BED into {name: (chrom, start, end, strand)}."""
    d = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            d[f[3]] = (f[0], int(f[1]), int(f[2]), f[5])
    return d

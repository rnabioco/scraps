""" Shared helpers for converting polyAdb releases to the scraps SAF format.

The scraps SAF format used as ``POLYA_SITES`` is a tab-separated table with a
header ``GeneID  Chr  Start  End  Strand`` where the GeneID encodes seven
fields:

    gene_symbol{D}refseq_gene_id{D}ensembl_id{D}chrom{D}pos{D}strand{D}pas_type

with delimiter ``{D}`` = ';' for human references and '_' for mouse references
(matching the bundled ref/polyadb32.{hg38,mm10}.saf.gz). ``pos`` is the 1-based
single-base cleavage position (field 5) and ``pas_type`` is the site class
(field 7). Missing values are encoded as 'NA'.

The SAF Start/End columns describe a strand-specific 16-base window
(End - Start = 15) around the cleavage position (transcript-oriented -10 / +5
slop):

    + strand: Start = pos - 10, End = pos + 5
    - strand: Start = pos - 5,  End = pos + 10

This module also provides an optional UCSC liftOver wrapper so converters can
remap coordinates to a target genome build. The liftOver binary is taken from
``--liftover``, else from PATH, else downloaded for the current platform.
"""

import gzip
import os
import platform
import shutil
import stat
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

# transcript-oriented slop applied around the cleavage position
SLOP_UP = 10
SLOP_DOWN = 5

SAF_HEADER = "GeneID\tChr\tStart\tEnd\tStrand"

# delimiter convention per species (matches bundled references)
SPECIES_DELIM = {"human": ";", "mouse": "_"}
VALID_DELIMS = (";", "_")
VALID_STRANDS = ("+", "-")


def check_delim(delim):
    """Return `delim` if it is one of the supported GeneID delimiters."""
    if delim not in VALID_DELIMS:
        sys.exit("--delim must be ';' (human) or '_' (mouse), got {!r}".format(delim))
    return delim


def check_strand(strand, context=""):
    if strand not in VALID_STRANDS:
        sys.exit("invalid strand {!r}{}; expected '+' or '-'".format(
            strand, (" at " + context) if context else ""))
    return strand


def parse_pos(value, context=""):
    """Parse a 1-based position, tolerating float / scientific notation
    (e.g. '7.7e+07' from a spreadsheet round-trip) like validate_saf.py."""
    try:
        return int(round(float(value)))
    except (TypeError, ValueError):
        sys.exit("cannot parse position {!r}{}".format(
            value, (" at " + context) if context else ""))


def normalise(val):
    """Return 'NA' for empty/na-like values, else the stripped value."""
    if val is None:
        return "NA"
    v = val.strip()
    if not v or v.lower() == "na":
        return "NA"
    return v


def encode_geneid(sym, refseq, ensembl, chrom, pos, strand, pas_type, delim):
    """Build the 7-field scraps GeneID string.

    Fails if any field already contains the delimiter: such a GeneID could not
    be split back into 7 fields by scraps' parsers (annotate_sites.py,
    validate_saf.py, the R helpers)."""
    fields = [
        normalise(sym),
        normalise(refseq),
        normalise(ensembl),
        chrom,
        str(pos),
        strand,
        normalise(pas_type),
    ]
    if delim == ";":
        bad = [f for f in fields if delim in f]
        if bad:
            sys.exit("GeneID field(s) {} contain the delimiter {!r} (site {}:{}:{}); "
                     "choose another --delim".format(bad, delim, chrom, pos, strand))
    # '_' legitimately occurs inside RefSeq ids and scaffold names; parsers
    # anchor on the SAF Chr/Strand columns to recover those, so only the fields
    # that would break that anchoring are checked.
    elif delim in strand or delim in normalise(pas_type):
        sys.exit("strand/pas_type must not contain '_' (site {}:{}:{})".format(
            chrom, pos, strand))
    return delim.join(fields)


def site_window(pos, strand):
    """Strand-specific SAF (Start, End) around a 1-based cleavage position.

    Mirrors the bundled polyAdb references: a 16-base window (End - Start = 15)
    with -10 / +5 transcript-oriented slop. Returned coordinates use the same
    convention as those references (Start = pos - up_slop, End = pos + down_slop
    on '+').
    """
    check_strand(strand)
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
    """Write SAF rows. Output is gzipped if ``path`` ends with '.gz'
    (with mtime=0 so identical inputs give identical bytes).

    ``rows`` is an iterable of (chrom, start, end, strand, gene_id).
    """
    path = str(path)
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    if path.endswith(".gz"):
        raw = open(path, "wb")
        out = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    else:
        raw = None
        out = open(path, "wb")
    try:
        out.write((SAF_HEADER + "\n").encode())
        for chrom, start, end, strand, gene_id in rows:
            out.write("{}\t{}\t{}\t{}\t{}\n".format(
                gene_id, chrom, start, end, strand).encode())
    finally:
        out.close()
        if raw is not None:
            raw.close()


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
# temp dir handling
# --------------------------------------------------------------------------

def make_tmpdir(user_tmpdir, prefix):
    """Return (Path, created). A user-supplied directory is created if missing
    but never deleted by cleanup_tmpdir; without one, a fresh mkdtemp is used."""
    if user_tmpdir:
        p = Path(user_tmpdir)
        p.mkdir(parents=True, exist_ok=True)
        return p, False
    return Path(tempfile.mkdtemp(prefix=prefix + ".")), True


def cleanup_tmpdir(tmpdir, created, keep):
    """Delete the temp dir only if this process created it and --keep-tmp is
    not set. User-supplied --tmpdir directories are left alone."""
    if keep:
        sys.stderr.write("  keeping temp dir {}\n".format(tmpdir))
        return
    if created and Path(tmpdir).exists():
        shutil.rmtree(tmpdir, ignore_errors=True)


# --------------------------------------------------------------------------
# liftOver support
# --------------------------------------------------------------------------

_UCSC_EXE = "https://hgdownload.soe.ucsc.edu/admin/exe/{}/liftOver"
_UCSC_PLATFORMS = {
    ("Darwin", "arm64"): "macOSX.arm64",
    ("Darwin", "x86_64"): "macOSX.x86_64",
    ("Linux", "x86_64"): "linux.x86_64",
}
_CHAIN_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/{frm}/liftOver/"
    "{frm}To{To}.over.chain.gz"
)


def liftover_url():
    """UCSC liftOver download URL for the running platform, or exit."""
    system, machine = platform.system(), platform.machine()
    if machine in ("amd64", "AMD64"):
        machine = "x86_64"
    key = (system, machine)
    if key not in _UCSC_PLATFORMS:
        sys.exit("no UCSC liftOver build for {}/{}; install liftOver (e.g. conda "
                 "ucsc-liftover) so it is on PATH, or pass --liftover".format(*key))
    return _UCSC_EXE.format(_UCSC_PLATFORMS[key])


def _download(url, dest, label):
    """Download to dest atomically (dest.part -> dest), reusing a complete
    existing dest."""
    dest = Path(dest)
    if dest.exists():
        sys.stderr.write("  {}: using existing {}\n".format(label, dest))
        return
    sys.stderr.write("  downloading {} from {}\n".format(label, url))
    part = dest.with_name(dest.name + ".part")
    try:
        urllib.request.urlretrieve(url, part)
    except Exception as e:  # noqa: BLE001 - report and exit cleanly
        if part.exists():
            part.unlink()
        sys.exit("download of {} failed: {}".format(url, e))
    os.replace(part, dest)


def ensure_liftover(tmpdir, liftover=None):
    """Return a path to a liftOver binary: --liftover, else PATH, else a
    platform-matched UCSC download into tmpdir."""
    if liftover:
        return liftover
    found = shutil.which("liftOver")
    if found:
        sys.stderr.write("  liftOver binary: using {} from PATH\n".format(found))
        return found
    lo = Path(tmpdir) / "liftOver"
    _download(liftover_url(), lo, "liftOver binary")
    lo.chmod(lo.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return str(lo)


def ensure_chain(tmpdir, source_build, target_build, chain=None):
    """Return a path to a decompressed chain file, downloading if needed."""
    if chain:
        if str(chain).endswith(".gz"):
            out = Path(tmpdir) / (Path(chain).stem)
            with gzip.open(chain, "rb") as src, open(out, "wb") as dst:
                shutil.copyfileobj(src, dst)
            return str(out)
        return chain
    to_cap = target_build[0].upper() + target_build[1:]
    url = _CHAIN_URL.format(frm=source_build, To=to_cap)
    chain_gz = Path(tmpdir) / "{}To{}.over.chain.gz".format(source_build, to_cap)
    chain_txt = Path(tmpdir) / "{}To{}.over.chain".format(source_build, to_cap)
    _download(url, chain_gz, "chain file")
    if not chain_txt.exists():
        with gzip.open(chain_gz, "rb") as src, open(chain_txt, "wb") as dst:
            shutil.copyfileobj(src, dst)
    return str(chain_txt)


def run_liftover(liftover_bin, chain, bed_in, bed_out, unmapped):
    """Run UCSC liftOver. liftOver exits 1 when some records are unmapped;
    returns the number of unmapped records (from the unmapped file)."""
    cmd = [liftover_bin, bed_in, chain, bed_out, unmapped]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
    except OSError as e:
        sys.exit("could not run liftOver binary {}: {}".format(liftover_bin, e))
    if result.returncode not in (0, 1):
        sys.exit("liftOver failed:\n{}".format(result.stderr))
    n_unmapped = 0
    if Path(unmapped).exists():
        with open(unmapped) as fh:
            n_unmapped = sum(1 for line in fh if line.strip() and not line.startswith("#"))
    return n_unmapped


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

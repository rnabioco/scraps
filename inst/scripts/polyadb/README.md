# polyAdb → scraps SAF converters

Helper scripts that convert [polyAdb](https://exon.apps.wistar.org/polya_db/)
releases into the scraps SAF format used as `POLYA_SITES`. All scripts use only
the Python 3 standard library (plus UCSC `liftOver`, auto-downloaded when
needed).

| Script | Input | Notes |
|---|---|---|
| `pas32_to_saf.py` | polyAdb 3.2 `*.PAS.txt` (human, mouse) | liftOver to a target build |
| `pasv4_to_saf.py` | polyAdb 4 `*.PAS.{main,max}.tsv` (human, mouse) | main/max, optional liftOver |
| `validate_saf.py` | a scraps SAF | format checks + optional comparison to a reference |
| `saf_common.py` | — | shared encoding/window/liftOver helpers (imported by the above) |

## The scraps SAF format

Tab-separated with header `GeneID  Chr  Start  End  Strand`. The **GeneID has
7 fields**:

```
gene_symbol{D}refseq_gene_id{D}ensembl_id{D}chrom{D}pos{D}strand{D}pas_type
```

- `{D}` is the delimiter: **`;` for human**, **`_` for mouse** (matches the
  bundled `ref/polyadb32.{hg38,mm10}.saf.gz`).
- `pos` (field 5) is the **1-based single-base cleavage position**.
- `pas_type` (field 7) is the **site class**.
- Missing values are `NA`.

`Start`/`End` describe a strand-specific 15 bp window (transcript-oriented
-10 / +5 slop) around the cleavage position:

```
+ strand: Start = pos - 10, End = pos + 5
- strand: Start = pos - 5,  End = pos + 10
```

Rows are deduplicated on `(chrom, pos, strand, GeneID)` and sorted by
`(chrom lexicographic, start numeric)`.

> **Note on the `_` (mouse) delimiter.** RefSeq IDs (e.g. `NR_152944`) and some
> scaffold chromosome names (e.g. `chrUn_GL456...`) contain underscores, so a
> naive split of a mouse GeneID can yield more than 7 tokens. scraps parses
> these robustly by anchoring on the SAF `Chr`/`Strand` columns and reading the
> trailing `chrom/pos/strand/class` fields from the right
> (`inst/scripts/annotate_sites.py`, `validate_saf.py`). The leading
> gene/refseq/ensembl group is passthrough metadata and does not affect site
> matching, which uses `pos`/`chrom`/`strand` only.

## Usage

### polyAdb 3.2 (`*.PAS.txt`)

polyAdb 3.2 human coordinates are hg19 and require liftOver to hg38; mouse is
already mm10.

```bash
# human: hg19 -> hg38 (auto-downloads liftOver binary + chain)
python3 pas32_to_saf.py --species human --pas human.PAS.txt \
  --source-build hg19 --target-build hg38 \
  --out ../../../ref/polyadb32.hg38.saf.gz

# mouse: already mm10 (no liftOver)
python3 pas32_to_saf.py --species mouse --pas mouse.PAS.txt \
  --source-build mm10 --target-build mm10 \
  --out ../../../ref/polyadb32.mm10.saf.gz
```

### polyAdb 4 (`*.PAS.{main,max}.tsv`)

v4 tables are already in the target build (hg38 / mm10). The `main` table
carries gene symbol + PAS type; the `max` table carries only the
highest-expressed gene per site (no PAS type). Mode is auto-detected from the
header (override with `--mode`).

```bash
python3 pasv4_to_saf.py --species human --pas HumanPas_v4/hg38.PAS.main.tsv \
  --source-build hg38 --target-build hg38 --out ref/polyadb4.hg38.main.saf.gz

python3 pasv4_to_saf.py --species mouse --pas MousePas_v4/mm10.PAS.max.tsv \
  --source-build mm10 --target-build mm10 --out ref/polyadb4.mm10.max.saf.gz
```

liftOver is available for v4 too (future builds): set `--source-build` and
`--target-build` to differing builds, e.g. `--source-build hg38
--target-build hs1`.

### liftOver options (both converters)

- `--liftover PATH` — use a local liftOver binary (else auto-download).
- `--chain PATH` — use a local `source→target` chain (`.gz` accepted; else
  auto-download by build name).
- `--tmpdir PATH`, `--keep-tmp` — control intermediate files.

If the source and target builds are equal, liftOver is skipped.

### Validation

```bash
# format only
python3 validate_saf.py polyadb32.hg38.saf.gz

# format + content comparison against a reference
python3 validate_saf.py generated.saf.gz ../../../ref/polyadb32.hg38.saf.gz
```

`.gz` inputs are read transparently.

## PAS type / class vocabulary

PAS type labels are passed through **verbatim** from each release, so the class
field differs between releases:

| | example class values |
|---|---|
| polyAdb 3.2 | `3'UTR(M)`, `3'UTR(L)`, `Intron`, `CDS`, `intergenic`, `Pseudogene`, `LncRNA(FANTOM5)` |
| polyAdb 4 main | `3'UTR`, `Intron`, `CDS`, `Intergenic`, `5'UTR`, `Downstream`, `Upstream`, `3'-most exon`, `Single exon` |
| polyAdb 4 max | `NA` (not provided) |

Downstream code that filters by class (e.g. `parse_saf(types=...)`,
discovery annotation) must account for these differences when mixing releases.

## Quantitation caveat

These SAFs are intended as the validated `POLYA_SITES` reference for scraps
quantitation. For the de novo discovery module's relationship to validated
sites, see [`docs/discovery.md`](../../../docs/discovery.md).

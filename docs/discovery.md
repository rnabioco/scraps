# De novo RT priming site discovery

> **Use validated poly(A) sites for quantitative analyses.** scraps quantifies
> poly(A) usage against curated, validated sites (polyAdb) by default, and this
> remains the recommended reference for all downstream **quantitative** work
> (per-cell counting, PSI, differential poly(A) usage). The discovery module
> described here is for **hypothesis generation and novel-site nomination**. Its
> output is a set of *candidate* sites that require orthogonal validation, and
> includes sites flagged as likely artifacts that should be **excluded**, not
> quantified.

## Motivation

scraps maps **reverse-transcription (RT) priming sites at single-base
resolution**. By default it quantifies these priming events only at validated
polyAdb cleavage/polyadenylation (C/PA) sites, because internal priming
(RT priming off genomically-encoded A-rich stretches rather than a true poly(A)
tail) is a major artifact of oligo-dT priming and is genuinely hard to identify
from sequence alone.

The default approach is safe but underuses a key strength of scraps: its
single-base resolution makes it possible to *discover* priming sites directly
from the data. This module calls RT priming sites de novo and then annotates
each one, separating validated/candidate C/PA sites from likely artifacts.

## What it does

1. **Pseudobulk, strand-aware, single-base priming signal.** For the chosen
   alignment mode (`DISCOVERY.read`), reads are UMI-deduplicated per cell
   (`umi_tools dedup`) and collapsed to single-base priming positions per
   strand (`bedtools genomecov`). Cell identity is intentionally collapsed —
   discovery operates on pseudobulk.

   The emitted strand is the **transcript strand**, not the read's own aligned
   strand. `bedtools genomecov` reports each position on the strand the read
   aligned to, so the strand is converted per alignment mode to match the
   validated `featureCounts` strandedness used for quantitation
   (`rules/count.snake`): **R2** is same-strand (`-s 1`), so the read-aligned
   strand already equals the transcript strand; **R1** and **paired** are
   reverse-stranded (`-s 2`), so the read-aligned strand is the *opposite* of
   the transcript strand and is flipped when the stranded bed is written
   (`rules/discovery.snake`, `stranded_bed`). All downstream steps (KDE per
   `(chrom, strand)`, polyAdb matching, sequence-context scanning, SAF slop)
   are strand-aware and consume this transcript strand.

2. **Optional sample pooling.** Samples can be pooled into named groups
   (`DISCOVERY.groups`) so that priming counts are summed across a set of
   samples before site calling. Ungrouped samples are processed individually.

3. **Kernel density estimation (KDE).** Per `(chromosome, strand)` track, the
   UMI-weighted per-base signal is smoothed with a Gaussian kernel
   (`kde_bandwidth`, in bp). Local maxima of the density are called as RT
   priming sites; maxima within `peak_merge_dist` are merged, and the reported
   summit is the highest-UMI base within the merged region.
   (`inst/scripts/kde_peaks.py`)

   Because discovery is **pseudobulk** — UMIs are deduplicated per cell, then
   summed across *all* cells and (with `groups`) across pooled samples — a raw
   UMI count of 1–2 at a base is background, not signal. A merged peak is kept
   only if its summed support clears

   ```
   umi_support >= max(min_umi, min_umi_frac * T)
   ```

   where `T` is the total UMI count across the whole (pooled) input. The
   **primary control** on peak count is the **absolute floor `min_umi`**. The
   depth-relative term `min_umi_frac * T` is **off by default** (`min_umi_frac:
   0`): because it scales *linearly* with total depth, on deep or pooled data it
   imposes a very high absolute bar (e.g. `1e-6 * 1e8 = 100` UMIs/peak) that
   silently drops usable mid-abundance peaks. Set a small positive
   `min_umi_frac` only if you deliberately want the bar to scale with depth.
   `min_density` is the **shape gate**: by default it is derived from
   `kde_bandwidth` (≈ 2 UMIs concentrated within one bandwidth) so that isolated
   single-UMI bumps cannot form a called maximum; set a positive value to
   override. (Density is evaluated sparsely, only at occupied bases.)

4. **Annotation into three categories.** Each summit is classified
   (`inst/scripts/annotate_sites.py`):

   | Status | Meaning | Use |
   |:--|:--|:--|
   | `known_pas` | Summit falls within the strand-aware `match_window` of a validated polyAdb cleavage point | Validated |
   | `likely_internal_priming` | Downstream genomic A-content / A-run indicates priming off a genomic poly(A) stretch | **Exclude (artifact)** |
   | `potential_novel_pas` | No internal-priming signature **and** a canonical poly(A) signal (PAS) motif present upstream | Candidate — validate before use |

## How sites are classified

- **Known vs novel:** the summit (1-based) is compared, strand-aware, to the
  polyAdb single-base cleavage positions. Matching is **point-based** against
  the cleavage position only — the polyAdb SAF `Start`/`End` interval is *not*
  used for matching. The match window is **asymmetric and transcript-oriented**
  (`match_window` = `[upstream, downstream]`, default `[10, 5]`): a summit is
  labeled `known_pas` if it lies within `upstream` bp upstream or `downstream`
  bp downstream of a polyAdb cleavage point, in the direction of transcription.
  The asymmetry reflects the skew of RT priming read distributions, and the
  default `[10, 5]` mirrors the slop of the provided polyAdb SAF references. A
  matched site inherits the polyAdb gene/id/class metadata, and the reported
  `polya_dist` is the **signed transcript-oriented** distance from the summit to
  the polyAdb point (negative = upstream, positive = downstream). Summits that
  fall outside the window proceed to the novel-vs-internal-priming branch below.

  > **Dependent on the scraps polyAdb SAF format** (see the spec below).
  > Known-site annotation reads the single-base cleavage position from
  > **field 5** of the `GeneID`. A SAF that does not follow this encoding will
  > fail to annotate known sites — every summit is treated as novel.
  > `annotate_sites.py` validates the GeneID field count and warns/errors if no
  > rows parse, so this failure mode is surfaced rather than silent. The bundled
  > `ref/polyadb32.{hg38,mm10}.saf.gz` already conform; converters for polyAdb
  > 3.2 and 4 are in `inst/scripts/polyadb/`.

- **Internal-priming signature (novel sites):** the genomic sequence
  *downstream* of the summit on the sense strand (window `a_content_window`,
  default 1–20 bp) is scored for A-content. A site is flagged
  `likely_internal_priming` if the A fraction is ≥ `a_content_max` (default 0.7)
  **or** it contains a contiguous A-run ≥ `--a-run-min` (default 8). This is the
  same biology that motivates using validated sites in the first place: a true
  poly(A) tail is added post-transcriptionally and is *not* templated, so a
  genomic A-stretch immediately downstream is the hallmark of an artifact.

- **PAS motif (novel sites):** the sequence *upstream* of the summit
  (window `pas_window`, default −40–0 bp, sense strand) is scanned for the
  canonical poly(A) signal `AATAAA` and common single-base variants. A novel
  site with no internal-priming signature and a detectable PAS is labeled
  `potential_novel_pas`. A novel site with neither a PAS nor a strong
  internal-priming signature is conservatively labeled
  `likely_internal_priming`.

## Reference SAF format

The polyAdb reference supplied as `POLYA_SITES` (and the de novo SAF emitted by
discovery) is a tab-separated table with header `GeneID  Chr  Start  End
Strand`. The **GeneID encodes seven fields**:

```
gene_symbol{D}refseq_gene_id{D}ensembl_id{D}chrom{D}pos{D}strand{D}pas_type
```

| Field | Index | Meaning |
|:--|:--|:--|
| `gene_symbol` | 1 | gene symbol (or `NA`) |
| `refseq_gene_id` | 2 | RefSeq/Entrez gene id (or `NA`) |
| `ensembl_id` | 3 | Ensembl gene id (or `NA`) |
| `chrom` | 4 | chromosome |
| `pos` | **5** | **1-based single-base cleavage position** (used for matching) |
| `strand` | 6 | `+` / `-` |
| `pas_type` | **7** | **site class** (used as the discovery `class`) |

- Delimiter `{D}` is **`;` for human** and **`_` for mouse** references (matching
  the bundled `ref/polyadb32.{hg38,mm10}.saf.gz`).
- `Start`/`End` form a strand-specific 15 bp window around `pos`
  (transcript-oriented -10 / +5 slop): `+` → `[pos-10, pos+5]`,
  `-` → `[pos-5, pos+10]`.

**Mouse `_` delimiter caveat:** RefSeq IDs (e.g. `NR_152944`) and scaffold
chromosome names (e.g. `chrUn_GL456...`) contain underscores, so a naive split
of a mouse GeneID over-counts fields. scraps parses these robustly by anchoring
on the SAF `Chr`/`Strand` columns and reading the trailing
`chrom/pos/strand/class` fields from the right; the leading gene/refseq/ensembl
group is passthrough metadata only.

**PAS type vocabulary differs by release** (labels are passed through verbatim
by the converters): polyAdb 3.2 uses values like `3'UTR(M)`, `3'UTR(L)`,
`Intron`, `intergenic`, `Pseudogene`, `LncRNA(FANTOM5)`; polyAdb 4 `main` uses
`3'UTR`, `Intron`, `Intergenic`, `5'UTR`, `Downstream`, `Upstream`,
`3'-most exon`, `Single exon`; polyAdb 4 `max` has no PAS type (`NA`). Code that
filters by class (e.g. `parse_saf(types=...)`) must account for this when mixing
releases.

**Generating references:** converters for polyAdb 3.2 (`*.PAS.txt`, with
liftOver) and polyAdb 4 (`*.PAS.{main,max}.tsv`) live in
[`inst/scripts/polyadb/`](../inst/scripts/polyadb/); see its README for usage.

### Caveats and limitations

- **Internal priming cannot be perfectly called by sequence alone.** The
  A-content heuristic will both miss some artifacts (false negatives) and flag
  some genuine sites near A-rich 3' UTRs (false positives). Thresholds are
  exposed in config and should be tuned per dataset.
- **PAS presence is suggestive, not definitive.** Many real sites use
  non-canonical signals; some artifact sites have an incidental upstream PAS.
- **`potential_novel_pas` are candidates, not validated sites.** Confirm with
  orthogonal evidence (3'-seq/QuantSeq, polyAdb updates, conservation, RNAse
  protection, etc.) before treating them as real.
- Discovery is **pseudobulk**: it does not provide per-cell single-base
  resolution. Per-cell quantitation still goes through the SAF-window counting
  path (`results/counts/`).
- **Peak count is controlled primarily by the absolute `min_umi`.** Raise
  `min_umi` for a smaller, higher-confidence set; lower it to nominate more
  candidates. The depth-relative `min_umi_frac` is **off by default** (`0`)
  because its linear-in-depth bar over-filters usable peaks on deep/pooled data;
  enable it (e.g. `1e-7`) only if you specifically want the floor to scale with
  sequencing depth. `min_density` mainly suppresses single-UMI artifacts.

## Enabling discovery

In `config.yaml`:

```yaml
GENOME_FASTA: "ref/genome.fa"     # must match the alignment reference; needs .fai

DISCOVERY:
  enabled: true
  read: R2                        # R1, R2, or paired
  groups:                         # optional; omit/empty for per-sample
    groupA: [chromiumv2_test, dropseq_test]
  kde_bandwidth: 10
  min_umi: 10                     # primary absolute support floor
  min_umi_frac: 0.0               # optional depth-relative term, off by default
                                  # (support >= max(min_umi, min_umi_frac * total_UMIs))
  min_density: 0.0                # 0 / omit -> derived from kde_bandwidth
  peak_merge_dist: 24
  match_window: [10, 5]           # [upstream, downstream] bp, transcript-oriented
  saf_slop: [10, 5]              # de novo SAF window; matches provided polyAdb slop
  pas_window: [-40, 0]
  a_content_window: [1, 20]
  a_content_max: 0.7
```

`match_window` and `saf_slop` are both transcript-oriented
`[upstream, downstream]` pairs. `match_window` controls how a discovered summit
is matched to a known polyAdb cleavage point (asymmetric, see above); `saf_slop`
controls the window written around each summit in the output de novo SAF. The
defaults `[10, 5]` reproduce the -10/+5 slop of the bundled polyAdb references,
so a de novo SAF re-used as `POLYA_SITES` produces feature windows consistent
with the originals.

`GENOME_FASTA` is **required** when `DISCOVERY.enabled` is true (the pipeline
will error otherwise) and must be indexed:

```bash
samtools faidx ref/genome.fa
```

Then run as usual; discovery targets are added automatically:

```bash
snakemake -npr --configfile config.yaml          # validate
snakemake --configfile config.yaml --keep-going  # run
```

## Outputs

For each discovery unit (group name or ungrouped sample) under
`results/discovery/`:

- **`<unit>_sites.tsv.gz`** — full annotation table:
  `chrom, summit, strand, kde_score, umi_support, status, gene, id, class,
  polya_dist, pas_motif, pas_dist, a_content, a_run`. `polya_dist` is the
  **signed transcript-oriented** distance (bp) from the summit to the matched
  polyAdb cleavage point (negative = upstream, positive = downstream;
  `NA` for unmatched/novel sites).
- **`<unit>_sites.bed.gz`** — BED6 of summits (`name` = status) for browsing.
- **`<unit>_denovo.saf.gz`** — a SAF using the same `GeneID` encoding as
  polyAdb (`gene;genbank;id;chrom;summit;strand;class`, where `class` is the
  status), windowed by `saf_slop` (transcript-oriented `[upstream, downstream]`,
  default `[10, 5]`).

## Re-quantifying against discovered sites (advanced)

Because `<unit>_denovo.saf.gz` uses the standard scraps SAF encoding, it can be
supplied as `POLYA_SITES` to re-run the normal per-cell quantitation against the
discovered sites — no code changes required:

```yaml
POLYA_SITES: "results/discovery/groupA_denovo.saf.gz"
```

> **Before doing this, filter the SAF.** Drop `likely_internal_priming`
> features, and treat `potential_novel_pas` features as unvalidated. For any
> publication-grade quantitative analysis, prefer the curated polyAdb reference.
> Re-quantifying against de novo sites is appropriate for exploratory work and
> for evaluating nominated novel sites, not as a default replacement for
> validated sites.

## Implementation

| Component | Location |
|:--|:--|
| Snakemake rules | `rules/discovery.snake` |
| KDE peak calling | `inst/scripts/kde_peaks.py` |
| Annotation + SAF | `inst/scripts/annotate_sites.py` |
| Config | `config.yaml` (`GENOME_FASTA`, `DISCOVERY`) |

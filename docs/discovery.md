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

2. **Optional sample pooling.** Samples can be pooled into named groups
   (`DISCOVERY.groups`) so that priming counts are summed across a set of
   samples before site calling. Ungrouped samples are processed individually.

3. **Kernel density estimation (KDE).** Per `(chromosome, strand)` track, the
   UMI-weighted per-base signal is smoothed with a Gaussian kernel
   (`kde_bandwidth`, in bp). Local maxima of the density above `min_density`
   are called as RT priming sites; maxima within `peak_merge_dist` are merged,
   and the reported summit is the highest-UMI base within the merged region.
   (`inst/scripts/kde_peaks.py`)

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

  > **Dependent on the scraps polyAdb SAF format.** Known-site annotation reads
  > the single-base cleavage position from **field 5** of the `;`- (hg38) or
  > `_`- (mm10) delimited `GeneID`
  > (`gene;genbank;id;chrom;pos;strand;class`). A SAF that does not follow this
  > encoding (or stores the position in a different field) will silently fail to
  > annotate any known sites — every summit would be treated as novel. The
  > bundled `ref/polyadb32.{hg38,mm10}.saf.gz` already conform. Conversion
  > helpers for polyAdb 3.2 and 4.0 are planned in future commits; until then,
  > custom references must match this encoding.

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
  min_density: 0.0
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

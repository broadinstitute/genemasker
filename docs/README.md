# genemasker Documentation

## Overview

`genemasker` generates Regenie-compatible gene/group mask files from VEP-style annotation data for rare variant burden testing.

The installed CLI entrypoint is:

```bash
genemasker
```

Package version in this repository: `1.0`.

## What It Produces

For output prefix `--out results/study1`, the pipeline writes:

- `results/study1.regenie.setlist.tsv`
- `results/study1.regenie.annotations.<mask>.tsv` (one file per mask)
- `results/study1.regenie.mask.<mask>.tsv` (one file per mask)
- `results/study1.results.tsv.gz` (full tabular results)
- `results/study1.genemasker.log` (or `results/study1.chunk<N>.genemasker.log` when `--chunk` is used)

It can also write intermediate artifacts (for reuse) and summary files such as:

- `*.rankscore_miss.tsv`
- `*.rankscore_maf_corr.tsv`
- `*.pc_maf_corr.tsv`
- `*.ic_maf_corr.tsv`
- `*.pca_explained_variance.tsv`
- `*.combined_score_corr.tsv`
- `*.damaging_prop.tsv`
- `*.impute_pipeline.pkl`, `*.pca_pipeline.pkl`, `*.ica_pipeline.pkl` (when `--save-all`)

## Input Requirements

### 1) Annotation file (`--annot`)

`genemasker` expects a tab-delimited VEP-like annotation file with the columns defined in `genemasker/definitions.py`.

- If the file name ends with `.bgz`, gzip compression is assumed.
- Otherwise pandas auto-detect (`infer`) is used.

Core required columns include:

- `#Uploaded_variation`, `Feature`, `Feature_type`, `Gene`, `PICK`, `Consequence`, `IMPACT`, `DOMAINS`, `LoF`
- Prediction/categorical fields used by masks (examples: `SIFT_pred`, `Polyphen2_HVAR_pred`, `clinvar_clnsig`)
- Rankscore fields used for imputation/PCA/ICA (many `*_rankscore` columns)
- `CADD_phred_hg19`, `REVEL_score`

### 2) Variant frequency/stat file (`--stat`)

A tab-delimited file joined to annotation rows by variant ID. You must provide:

- `--stat`
- `--stat-id-col`
- At least one of: `--stat-maf-col`, `--stat-mac-col`

In practice, masks and downstream correlation steps depend on MAF, so provide `--stat-maf-col` for standard use.

### 3) Optional user extension files

- `--user-definitions`: Python file loaded as a module for custom definitions.
- `--user-defined-filters`: Python file with additional mask/filter functions.

## Pipeline Modes

`genemasker` supports three execution patterns.

### Mode A: Full run from annotation

Reads annotation + stat data, computes combined damage scores, applies masks, and writes Regenie files.

### Mode B: Resume from scored parquet

Use `--generate-from-scored` to skip score computation and regenerate filters/group files.

### Mode C: Resume from filtered parquet

Use `--generate-from-filtered` to skip directly to group file generation.

## Argument Validation Rules

- `--out` is always required.
- If **neither** `--generate-from-scored` nor `--generate-from-filtered` is used, you must provide:
  - full projection reuse set:
    - `--impute-pipeline`, `--pca-pipeline`, `--ica-pipeline`, `--rankscore-maf-corr`, `--pc-maf-corr`, `--ic-maf-corr`, `--rankscore-miss`
  - and/or stat-file set:
    - `--stat`, `--stat-id-col`, and at least one of `--stat-maf-col`, `--stat-mac-col`
- `--run-masks-file` and `--run-masks` are mutually exclusive.
- If `--chunk` is supplied, group file generation is skipped (chunk-level processing only).

## Quick Start

```bash
genemasker \
  --annot data/vep.annot.tsv.bgz \
  --stat data/variant_stats.tsv \
  --stat-id-col variant_id \
  --stat-maf-col maf \
  --stat-mac-col mac \
  --out results/study1
```

## Resume Examples

From scored parquet files:

```bash
genemasker \
  --generate-from-scored 'results/study1_tmp/*.scored.parquet' \
  --out results/study1_rerun
```

From filtered parquet files:

```bash
genemasker \
  --generate-from-filtered 'results/study1_tmp/*.filters.parquet' \
  --out results/study1_rerun
```

## Custom Mask Selection

Run only selected masks by name:

```bash
genemasker \
  --annot data/vep.annot.tsv.bgz \
  --stat data/variant_stats.tsv \
  --stat-id-col variant_id \
  --stat-maf-col maf \
  --run-masks new_damaging_og25,x37348876_m8 \
  --out results/study1_subset
```

Or from a file (`one mask name per line`):

```bash
genemasker \
  --annot data/vep.annot.tsv.bgz \
  --stat data/variant_stats.tsv \
  --stat-id-col variant_id \
  --stat-maf-col maf \
  --run-masks-file config/masks.txt \
  --out results/study1_subset
```

## Chromosome Recoding Example

```bash
genemasker \
  --annot data/vep.annot.tsv.bgz \
  --stat data/variant_stats.tsv \
  --stat-id-col variant_id \
  --stat-maf-col maf \
  --recode-chrs '{"^X$":"23","^Y$":"24","^MT$":"26"}' \
  --out results/study1_recode
```

## Transcripts and Conserved Domains

Use all transcript-level annotations and restrict to conserved domains:

```bash
genemasker \
  --annot data/vep.annot.tsv.bgz \
  --stat data/variant_stats.tsv \
  --stat-id-col variant_id \
  --stat-maf-col maf \
  --include-transcripts \
  --conserved-domains-only \
  --out results/study1_tx_domains
```

## Chunked Execution Example

```bash
genemasker \
  --annot data/vep.annot.tsv.bgz \
  --stat data/variant_stats.tsv \
  --stat-id-col variant_id \
  --stat-maf-col maf \
  --chunk-size 200000 \
  --chunk 3 \
  --out results/study1
```

## Intermediate Artifact Reuse Example

```bash
genemasker \
  --annot data/vep.annot.tsv.bgz \
  --impute-pipeline results/study1.impute_pipeline.pkl \
  --pca-pipeline results/study1.pca_pipeline.pkl \
  --ica-pipeline results/study1.ica_pipeline.pkl \
  --rankscore-miss results/study1.rankscore_miss.tsv \
  --rankscore-maf-corr results/study1.rankscore_maf_corr.tsv \
  --pc-maf-corr results/study1.pc_maf_corr.tsv \
  --ic-maf-corr results/study1.ic_maf_corr.tsv \
  --out results/study1_reuse
```

## Built-in Masks

See full list in `docs/masks.md`.

## Full CLI

See `docs/cli-reference.md` for every option and examples.

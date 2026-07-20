# Mask Production Workflows

This guide describes each supported place to begin a `genemasker` run. Choose
the latest checkpoint that still contains the information needed for the change
you want to make.

| Starting point | Use when | Required input | Work performed |
| --- | --- | --- | --- |
| Raw annotation | No prior model/checkpoint is available, or scoring must be retrained | `--annot` and a stat file | Read/filter annotation, train score models, score, filter, write Regenie files |
| Projection | A previous run trained the score models and a new annotation should use that same model | `--annot`, stat file, seven saved projection artifacts | Read/filter annotation, project models, filter, write Regenie files |
| Scored parquet | Scores are already present but masks need changing | `--generate-from-scored` glob | Apply filters, write Regenie files |
| Filtered parquet | Mask Boolean columns are already present and only final files must be regenerated | `--generate-from-filtered` glob | Write Regenie files |

Quote glob patterns so the shell passes the pattern to `genemasker` unchanged.

## 1. Full run from raw annotation data

Use this workflow for a new dataset or whenever imputation, PCA, ICA, and the
combined damage scores should be fitted again. `--save-all` is recommended when
you may later reuse the trained model or scored parquet files.

```bash
genemasker \
  --annot data/study1.vep.tsv.bgz \
  --stat data/study1.variant_stats.tsv \
  --stat-id-col variant_id \
  --stat-maf-col maf \
  --stat-mac-col mac \
  --save-all \
  --out results/study1
```

With no `--run-masks` or `--run-masks-file`, every built-in mask is calculated.
To calculate only named masks, give a comma-separated list:

```bash
genemasker \
  --annot data/study1.vep.tsv.bgz \
  --stat data/study1.variant_stats.tsv \
  --stat-id-col variant_id \
  --stat-maf-col maf \
  --run-masks new_damaging_og25,x37348876_m8 \
  --save-all \
  --out results/study1_selected
```

The reusable artifacts from a `--save-all` run are:

```text
results/study1.impute_pipeline.pkl
results/study1.pca_pipeline.pkl
results/study1.ica_pipeline.pkl
results/study1.rankscore_miss.tsv
results/study1.rankscore_maf_corr.tsv
results/study1.pc_maf_corr.tsv
results/study1.ic_maf_corr.tsv
results/study1_tmp/*.scored.parquet
results/study1_tmp/*.filters.parquet
```

## 2. Project a prior model onto a new annotation file

Projection reuses the model and correlation choices from a prior run while
scoring a new annotation dataset. Supply the annotation file and all seven
projection artifacts together. Supply the stat file too when the masks or
correlation calculations need `MAF`/`MAC` (the usual case).

The new annotation must be compatible with the prior model: it needs the same
rank-score features retained by the original run and the same annotation field
semantics. Projection is not a way to add or remove rank-score inputs.

```bash
genemasker \
  --annot data/study2.vep.tsv.bgz \
  --stat data/study2.variant_stats.tsv \
  --stat-id-col variant_id \
  --stat-maf-col maf \
  --impute-pipeline results/study1.impute_pipeline.pkl \
  --pca-pipeline results/study1.pca_pipeline.pkl \
  --ica-pipeline results/study1.ica_pipeline.pkl \
  --rankscore-miss results/study1.rankscore_miss.tsv \
  --rankscore-maf-corr results/study1.rankscore_maf_corr.tsv \
  --pc-maf-corr results/study1.pc_maf_corr.tsv \
  --ic-maf-corr results/study1.ic_maf_corr.tsv \
  --save-all \
  --out results/study2_projected
```

`--save-all` on the source run is what creates the three `.pkl` pipeline files.
The four TSV summaries are written during scoring, but retaining them with the
source results is necessary for a reproducible projection.

## 3. Start from scored parquet files

Scored parquet contains the retained annotation columns plus combined damage
scores, but not necessarily the requested mask Boolean columns. This is the
right restart point for changing mask selection, adding a new filter file, or
changing mask logic without recomputing scores.

```bash
genemasker \
  --generate-from-scored 'results/study1_tmp/*.scored.parquet' \
  --run-masks new_damaging_og25,x37348876_m8 \
  --out results/study1_remasked
```

Use `--save-all` on this remasking command if the scored parquet must remain
available afterwards; otherwise each input `*.scored.parquet` is removed after
its `*.filters.parquet` replacement is written.

```bash
genemasker \
  --generate-from-scored 'results/study1_tmp/*.scored.parquet' \
  --user-defined-filters custom/study_filters.py \
  --run-masks study_lof_or_high_cadd \
  --save-all \
  --out results/study1_custom
```

## 4. Start from filtered parquet files

Filtered parquet already contains the mask Boolean columns. This workflow only
creates the Regenie files and `results.tsv.gz`; it does not rerun filter
functions. Use it to recreate outputs from an existing filtered checkpoint.

```bash
genemasker \
  --generate-from-filtered 'results/study1_tmp/*.filters.parquet' \
  --run-masks new_damaging_og25,x37348876_m8 \
  --out results/study1_regenerated
```

The requested mask names must correspond to columns in the parquet and to
registered filters in the current process. For a custom mask, load the same
custom filter file even though its function will not be evaluated:

```bash
genemasker \
  --generate-from-filtered 'results/study1_custom_tmp/*.filters.parquet' \
  --user-defined-filters custom/study_filters.py \
  --run-masks study_lof_or_high_cadd \
  --out results/study1_custom_regenerated
```

## 5. Chunked processing

`--chunk` selects one checkpoint path by replacing `*` in the supplied glob and
skips final group-file generation. It is intended for processing one chunk at a
time, followed by a non-chunked filtered-parquet aggregation run.

For example, remask checkpoint number 3:

```bash
genemasker \
  --generate-from-scored 'results/study1_tmp/study1.vep.tsv.bgz-*.scored.parquet' \
  --chunk 3 \
  --run-masks new_damaging_og25 \
  --out results/study1_chunked
```

After every chunk has produced a filtered parquet file, create the final output
without `--chunk`:

```bash
genemasker \
  --generate-from-filtered 'results/study1_chunked_tmp/*.filters.parquet' \
  --run-masks new_damaging_og25 \
  --out results/study1_chunked_final
```

## Common controls

Use `--include-transcripts` to group by gene/transcript instead of the VEP
`PICK` row, and `--conserved-domains-only` to retain only annotations with a
conserved-domain value. `--recode-chrs` applies regex-based chromosome recoding
when final files are written:

```bash
genemasker \
  --generate-from-filtered 'results/study1_tmp/*.filters.parquet' \
  --recode-chrs '{"^X$":"23","^Y$":"24","^MT$":"26"}' \
  --out results/study1_recoded
```

Use exactly one of `--run-masks` and `--run-masks-file`. A masks file contains
one function name per line. Names not registered by a built-in or loaded custom
filter are not selected, so check the run log and expected output files when
using a new name.

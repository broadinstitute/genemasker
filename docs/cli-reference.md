# CLI Reference

This reference documents all command-line options exposed in `genemasker/config.py`.

## Command

```bash
genemasker [OPTIONS] --out <prefix>
```

## Global

### `--help`

Show help.

Example:

```bash
genemasker --help
```

### `--version`

Show version and exit.

Example:

```bash
genemasker --version
```

### `--out` (required)

Output file prefix.

Example:

```bash
genemasker --generate-from-filtered 'results/run_tmp/*.filters.parquet' --out results/run
```

## Input/Resume Source Options

### `--annot`

Annotation input file path.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col variant_id --stat-maf-col maf --out results/run
```

### `--generate-from-scored`

Glob for precomputed `*.scored.parquet` inputs.

Example:

```bash
genemasker --generate-from-scored 'results/run_tmp/*.scored.parquet' --out results/run_regen
```

### `--generate-from-filtered`

Glob for precomputed `*.filters.parquet` inputs.

Example:

```bash
genemasker --generate-from-filtered 'results/run_tmp/*.filters.parquet' --out results/run_regen
```

## Stat/MAF Integration

### `--stat`

Stat (MAF/MAC) TSV path.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col variant_id --stat-maf-col maf --out results/run
```

### `--stat-id-col`

Column in stat file matching `#Uploaded_variation`.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --out results/run
```

### `--stat-maf-col`

MAF column in stat file.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --out results/run
```

### `--stat-mac-col`

MAC column in stat file.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --stat-mac-col MAC --out results/run
```

## Model/Pipeline Reuse

### `--impute-pipeline`

Load imputation pipeline from prior run.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --impute-pipeline results/base.impute_pipeline.pkl --pca-pipeline results/base.pca_pipeline.pkl --ica-pipeline results/base.ica_pipeline.pkl --rankscore-miss results/base.rankscore_miss.tsv --rankscore-maf-corr results/base.rankscore_maf_corr.tsv --pc-maf-corr results/base.pc_maf_corr.tsv --ic-maf-corr results/base.ic_maf_corr.tsv --out results/reuse
```

### `--pca-pipeline`

Load PCA pipeline from prior run.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --impute-pipeline results/base.impute_pipeline.pkl --pca-pipeline results/base.pca_pipeline.pkl --ica-pipeline results/base.ica_pipeline.pkl --rankscore-miss results/base.rankscore_miss.tsv --rankscore-maf-corr results/base.rankscore_maf_corr.tsv --pc-maf-corr results/base.pc_maf_corr.tsv --ic-maf-corr results/base.ic_maf_corr.tsv --out results/reuse
```

### `--ica-pipeline`

Load ICA pipeline from prior run.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --impute-pipeline results/base.impute_pipeline.pkl --pca-pipeline results/base.pca_pipeline.pkl --ica-pipeline results/base.ica_pipeline.pkl --rankscore-miss results/base.rankscore_miss.tsv --rankscore-maf-corr results/base.rankscore_maf_corr.tsv --pc-maf-corr results/base.pc_maf_corr.tsv --ic-maf-corr results/base.ic_maf_corr.tsv --out results/reuse
```

### `--rankscore-miss`

Provide precomputed rankscore missingness file.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --impute-pipeline results/base.impute_pipeline.pkl --pca-pipeline results/base.pca_pipeline.pkl --ica-pipeline results/base.ica_pipeline.pkl --rankscore-miss results/base.rankscore_miss.tsv --rankscore-maf-corr results/base.rankscore_maf_corr.tsv --pc-maf-corr results/base.pc_maf_corr.tsv --ic-maf-corr results/base.ic_maf_corr.tsv --out results/reuse
```

### `--rankscore-maf-corr`

Provide precomputed rankscore vs MAF correlations.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --impute-pipeline results/base.impute_pipeline.pkl --pca-pipeline results/base.pca_pipeline.pkl --ica-pipeline results/base.ica_pipeline.pkl --rankscore-miss results/base.rankscore_miss.tsv --rankscore-maf-corr results/base.rankscore_maf_corr.tsv --pc-maf-corr results/base.pc_maf_corr.tsv --ic-maf-corr results/base.ic_maf_corr.tsv --out results/reuse
```

### `--pc-maf-corr`

Provide precomputed PC vs MAF correlations.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --impute-pipeline results/base.impute_pipeline.pkl --pca-pipeline results/base.pca_pipeline.pkl --ica-pipeline results/base.ica_pipeline.pkl --rankscore-miss results/base.rankscore_miss.tsv --rankscore-maf-corr results/base.rankscore_maf_corr.tsv --pc-maf-corr results/base.pc_maf_corr.tsv --ic-maf-corr results/base.ic_maf_corr.tsv --out results/reuse
```

### `--ic-maf-corr`

Provide precomputed IC vs MAF correlations.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --impute-pipeline results/base.impute_pipeline.pkl --pca-pipeline results/base.pca_pipeline.pkl --ica-pipeline results/base.ica_pipeline.pkl --rankscore-miss results/base.rankscore_miss.tsv --rankscore-maf-corr results/base.rankscore_maf_corr.tsv --pc-maf-corr results/base.pc_maf_corr.tsv --ic-maf-corr results/base.ic_maf_corr.tsv --out results/reuse
```

## Training/Computation Controls

### `--chunk-size`

Rows per annotation chunk.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --chunk-size 200000 --out results/run
```

### `--chunk`

Run a specific chunk index; skips group-file generation.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --chunk-size 200000 --chunk 2 --out results/run
```

### `--pca-fit-method`

PCA fitting method: `standard` or `incremental`.

Example (`standard`):

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --pca-fit-method standard --out results/run
```

Example (`incremental`):

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --pca-fit-method incremental --out results/run
```

### `--pca-training-frac`

Fraction of variants used for PCA training.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --pca-training-frac 0.25 --out results/run
```

### `--impute-training-frac`

Fraction of variants used for imputer training.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --impute-training-frac 0.25 --out results/run
```

### `--ica-training-frac`

Fraction of variants used for ICA training.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --ica-training-frac 0.25 --out results/run
```

### `--impute-standardized`

Standardize inputs before imputer training.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --impute-standardized --out results/run
```

### `--pca-standardized`

Standardize inputs before PCA training.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --pca-standardized --out results/run
```

### `--impute-tol`

Iterative imputer convergence tolerance.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --impute-tol 0.0005 --out results/run
```

## Filtering and Annotation Scope

### `--conserved-domains-only`

Keep only rows with non-null conserved domain annotations.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --conserved-domains-only --out results/run
```

### `--include-transcripts`

Use transcript-level grouping instead of `PICK == 1` canonical rows.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --include-transcripts --out results/run
```

### `--recode-chrs`

JSON dictionary for regex-based chromosome recoding.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --recode-chrs '{"^X$":"23","^Y$":"24"}' --out results/run
```

## Mask Selection

### `--run-masks`

Comma-separated mask function names.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --run-masks new_damaging_ic25,x31383942_m4 --out results/run
```

### `--run-masks-file`

File of mask names (one per line).

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --run-masks-file config/masks.txt --out results/run
```

## Custom Extensions

### `--user-definitions`

Load custom definitions module.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --user-definitions custom/definitions.py --out results/run
```

### `--user-defined-filters`

Load custom filters module (function names must not conflict with built-ins).

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --user-defined-filters custom/filters.py --out results/run
```

## Output/Runtime Behavior

### `--skip-calc-perc-damaging`

Skip `*.damaging_prop.tsv` computation.

Example:

```bash
genemasker --generate-from-scored 'results/run_tmp/*.scored.parquet' --skip-calc-perc-damaging --out results/run
```

### `--save-all`

Retain temporary intermediate parquet files and save model pipelines.

Example:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --save-all --out results/run
```

## Common Option Combinations

Full run from raw inputs:

```bash
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --stat-mac-col MAC --out results/run
```

Run from existing filtered chunks and write Regenie outputs only:

```bash
genemasker --generate-from-filtered 'results/base_tmp/*.filters.parquet' --run-masks x37348876_m8 --out results/run
```

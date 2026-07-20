# genemasker Documentation

`genemasker` produces Regenie-compatible gene/group mask files for rare-variant
analysis from VEP-style annotations.  It can train damage-scoring models from an
annotation file, project a previously trained model onto a new annotation file,
or resume from saved parquet checkpoints.

The installed command is:

```bash
genemasker
```

## Start here

- [Mask production workflows](workflows.md) explains the raw-data, projection,
  scored-parquet, filtered-parquet, and chunked workflows, including what each
  one reads and writes.
- [Custom masks and annotation definitions](custom-masks.md) explains how to
  supply Python files through `--user-defined-filters` and `--user-definitions`.
- [Built-in masks](masks.md) lists the masks included with the package.
- [CLI reference](cli-reference.md) describes every command-line option.

## Inputs

### Annotation file

`--annot` is a tab-delimited VEP-style annotation file. Its required columns
and dtypes are defined in `genemasker/definitions.py`; the reader requests all
of those columns. Files ending in `.bgz` are read as gzip-compressed files.

At a minimum, the standard definitions include variant ID (`#Uploaded_variation`),
gene and transcript fields, consequence/impact fields, LoF and predictor fields,
many rank-score columns, `CADD_phred_hg19`, and `REVEL_score`. The variant ID is
expected to be colon-delimited (for example, `chr1:12345:A:G`) so the final
Regenie files can be ordered by chromosome and position.

By default, the `PICK == "1"` VEP row is used for each variant. `--include-transcripts`
instead retains transcript annotations and makes each gene/transcript combination
a separate group. `--conserved-domains-only` keeps only rows with a non-null
`DOMAINS` value.

### Stat file

When calculating scores from an annotation file, provide a tab-delimited stat
file with `--stat`, `--stat-id-col`, and at least one of `--stat-maf-col` or
`--stat-mac-col`. The selected fields are renamed internally to `MAF` and `MAC`.
Most built-in masks use `MAF`, and the score-correlation steps also use it, so a
normal run should supply `--stat-maf-col`.

## Final outputs

For `--out results/study1`, a completed non-chunked run writes:

- `results/study1.regenie.setlist.tsv` — gene/group definitions and their variants.
- `results/study1.regenie.annotations.<mask>.tsv` — variants assigned to each mask.
- `results/study1.regenie.mask.<mask>.tsv` — Regenie mask declaration for each mask.
- `results/study1.results.tsv.gz` — the final annotation table, including mask columns.
- `results/study1.genemasker.log` — run log.

Score summaries are also written as applicable: `rankscore_miss.tsv`,
`rankscore_maf_corr.tsv`, `pc_maf_corr.tsv`, `ic_maf_corr.tsv`,
`pca_explained_variance.tsv`, `combined_score_corr.tsv`, and `damaging_prop.tsv`.

Use `--save-all` when the run must be reused. It preserves intermediate parquet
files and writes the fitted imputation, PCA, and ICA pipelines; without it, the
pipeline removes several prior-stage parquet files while advancing.

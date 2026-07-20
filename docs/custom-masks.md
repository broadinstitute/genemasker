# Custom Masks and Annotation Definitions

`genemasker` can load two user-supplied Python files before processing begins:

- `--user-definitions` extends the annotation columns read from the VEP-style
  file.
- `--user-defined-filters` registers additional mask functions.

The definition file is loaded first, then the filter file. Use both on a raw or
projection run when a custom mask depends on an extra annotation column. When
remasking scored parquet, only the filter file is needed if that column is
already in the parquet. When generating from filtered parquet, load the filter
file again so `genemasker` recognizes the custom mask name.

## 1. Add input annotation fields with a definitions file

The annotation reader only requests fields present in the shared `annot_cols`
mapping. A definitions file can extend that mapping before the annotation file
is read. It can also add an annotation missing-value convention.

Create `custom/study_definitions.py`:

```python
from genemasker.definitions import annot_cols, annot_na_values

# Column name in the VEP-style TSV: pandas dtype used when it is read.
annot_cols["study_damage_score"] = "float64"

# Treat the VEP-style '-' placeholder as missing for this field.
annot_na_values["study_damage_score"] = "-"
```

Every column added to `annot_cols` must exist in the annotation input. Use a
pandas-compatible dtype such as `"str"` or `"float64"`. A custom field ending
in `rankscore` becomes an input to the imputation/PCA/ICA workflow, so only use
that suffix when the new field is intentionally part of score modelling.

Name the file something other than `definitions.py`; a distinctive name such as
`study_definitions.py` avoids module-name collisions during dynamic loading.

## 2. Define and register a mask

A custom filter file imports the built-in helper functions and `mask` decorator
from the module named `filters`. Decorating a function registers it in the
package-wide mask list. The function receives a pandas DataFrame and must
return a Boolean Series aligned to that DataFrame.

Create `custom/study_filters.py`:

```python
from filters import CADD_phred_20, LoF_HC, mask, missense


def study_damage_score_080(row):
    return row["study_damage_score"] >= 0.80


@mask
def study_lof_or_high_cadd(df):
    """LoF HC variants, or missense variants with CADD >= 20."""
    return df.apply(
        lambda row: LoF_HC(row) or (missense(row) and CADD_phred_20(row)),
        axis=1,
    )


@mask
def study_damage_080(df):
    return df.apply(lambda row: study_damage_score_080(row), axis=1)
```

Use a unique function name for every decorated mask. It must not overlap a
built-in filter function name; the loader raises an error when it detects such
a collision. Only functions decorated with `@mask` become selectable masks.
Helper functions may have any unique name and do not need the decorator.

Do not name this file `filters.py`. The loader exposes the built-in module as
`filters` for imports, so a distinct filename such as `study_filters.py` is
needed for the example import to reliably refer to built-in helpers.

## 3. Run a custom mask from raw data

This run reads the added `study_damage_score` column, registers the custom
masks, calculates only `study_damage_080`, and retains checkpoints for reuse.

```bash
genemasker \
  --annot data/study1.vep.tsv.bgz \
  --stat data/study1.variant_stats.tsv \
  --stat-id-col variant_id \
  --stat-maf-col maf \
  --user-definitions custom/study_definitions.py \
  --user-defined-filters custom/study_filters.py \
  --run-masks study_damage_080 \
  --save-all \
  --out results/study1_custom
```

To produce both custom and built-in masks, list all desired names:

```bash
genemasker \
  --annot data/study1.vep.tsv.bgz \
  --stat data/study1.variant_stats.tsv \
  --stat-id-col variant_id \
  --stat-maf-col maf \
  --user-definitions custom/study_definitions.py \
  --user-defined-filters custom/study_filters.py \
  --run-masks study_damage_080,new_damaging_og25 \
  --out results/study1_custom_and_builtin
```

## 4. Custom masks from checkpoints

To add a mask that uses only fields already present in a scored parquet file,
load the filter file and start from scored parquet. The definitions file is not
needed because no annotation TSV is read.

```bash
genemasker \
  --generate-from-scored 'results/study1_custom_tmp/*.scored.parquet' \
  --user-defined-filters custom/study_filters.py \
  --run-masks study_damage_080 \
  --save-all \
  --out results/study1_custom_remasked
```

If the filtered parquet already contains `study_damage_080`, only final Regenie
files need to be regenerated. Still load the filter file so the requested name
is registered:

```bash
genemasker \
  --generate-from-filtered 'results/study1_custom_tmp/*.filters.parquet' \
  --user-defined-filters custom/study_filters.py \
  --run-masks study_damage_080 \
  --out results/study1_custom_regenerated
```

## Verification checklist

After a custom run, confirm that the log lists the supplied extension paths and
that the requested output files exist:

```text
results/study1_custom.regenie.annotations.study_damage_080.tsv
results/study1_custom.regenie.mask.study_damage_080.tsv
```

The mask column is also present in `results/study1_custom.results.tsv.gz` and
in the saved `*.filters.parquet` files. Before using outputs for analysis,
inspect the count of `True` values and a small set of included variants to make
sure the custom field's encoding and missing-value rule match the intended mask.

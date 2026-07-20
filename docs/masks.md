# Built-in Masks

The built-in baseline masks are curated into three selectable sets:

- `all` (default): all 16 baseline masks (and two pLoF masks for interpretability).
- `low-frequency`: the 8 masks in the **Baseline low-frequency** strategy (and low-frequency pLoF for interpretability).
- `rare`: the 8 masks in the **Baseline rare** strategy (and rare pLoF for interpretability).

Select a set with `--mask-set`. `--mask-set` cannot be combined with
`--run-masks` or `--run-masks-file`; use an explicit list when selecting a
custom or one-off subset.

```bash
# All curated baseline masks (the default)
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --out results/all

# Low-frequency masks (MAF < 1%)
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --mask-set low-frequency --out results/low_frequency

# Rare masks (generally MAF < 0.1%)
genemasker --annot data/vep.annot.tsv.bgz --stat data/stat.tsv --stat-id-col ID --stat-maf-col MAF --mask-set rare --out results/rare
```

## Low-frequency masks

| Mask ID (in group files) | Summary |
| --- | --- |
| `x24507775_m6_0_01` | pLoF variants or damaging missense predicted by PolyPhen2; MAF < 1%. |
| `x29177435_m1` | Exonic variants; MAF < 1%. |
| `x29378355_m1_0_01` | pLoF or missense variants with 0.5% <= MAF < 1%. |
| `x32141622_m7_0_01` | pLoF or splice-site variants; MAF < 1%. |
| `x32853339_m1` | High-impact or (likely) pathogenic ClinVar variants; MAF < 1%. |
| `x36327219_m3` | High-impact, indels, or missense called damaging by five algorithms; MAF < 1%. |
| `new_damaging_og25_0_01` | pLoF or missense predicted damaging by at least 25% of algorithms; MAF < 1%. |
| `new_damaging_og50_0_01` | pLoF or missense predicted damaging by at least 50% of algorithms; MAF < 1%. |
| `LoF_HC_0_01` | High-confidence loss-of-function variants; MAF < 1%. |

## Rare masks

| Mask ID (in group files) | Summary |
| --- | --- |
| `x30269813_m4` | pLoF, indels, or damaging missense predicted by PolyPhen2; MAF < 0.1%. |
| `x31118516_m5_0_001` | pLoF or missense called damaging by five algorithms; MAF < 0.1%. |
| `x31383942_m10` | pLoF or (likely) pathogenic ClinVar variants; MAF < 0.1%. |
| `x31383942_m4` | pLoF or REVEL-damaging missense variants; MAF < 0.1%. |
| `x34183866_m1` | pLoF, indels, or PolyPhen2-damaging missense variants; MAF < 0.01%. |
| `x34216101_m3_0_001` | PolyPhen2-damaging missense variants; MAF < 0.1%. |
| `x36411364_m4_0_001` | pLoF or REVEL-damaging missense variants; MAF < 0.1%. |
| `x37348876_m8` | pLoF or CADD-damaging variants; MAF < 0.1%. |
| `LoF_HC_0_001` | High-confidence loss-of-function variants; MAF < 0.1%. |

## Explicit selection

`--run-masks` accepts a comma-separated list of built-in or user-defined mask
function names. Every requested name is validated before output is written.

```bash
genemasker \
  --annot data/vep.annot.tsv.bgz \
  --stat data/stat.tsv \
  --stat-id-col ID \
  --stat-maf-col MAF \
  --run-masks LoF_HC_0_001,x37348876_m8 \
  --out results/custom
```

`--run-masks-file` provides the same selection in a file with one mask name per
line.

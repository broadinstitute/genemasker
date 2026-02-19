# Built-in Masks

These are built-in filter function names you can pass to `--run-masks` or include in `--run-masks-file`.

## Mask Names

- `new_damaging_ic25`
- `new_damaging_og25`
- `new_damaging_og25_0_01`
- `new_damaging_og50`
- `new_damaging_og50_0_01`
- `x23633568_m1`
- `x24507775_m6_0_01`
- `x29177435_m1`
- `x29378355_m1_0_01`
- `x30269813_m4`
- `x30828346_m1`
- `x31118516_m5_0_001`
- `x31383942_m10`
- `x31383942_m4`
- `x32141622_m4`
- `x32141622_m7`
- `x32141622_m7_0_01`
- `x32853339_m1`
- `x34183866_m1`
- `x34216101_m3_0_001`
- `x36327219_m3`
- `x36411364_m4_0_001`
- `x37348876_m8`

## `--run-masks` Example

```bash
genemasker \
  --annot data/vep.annot.tsv.bgz \
  --stat data/stat.tsv \
  --stat-id-col ID \
  --stat-maf-col MAF \
  --run-masks new_damaging_ic25,new_damaging_og50,x37348876_m8 \
  --out results/run
```

## `--run-masks-file` Example

Create `config/masks.txt`:

```text
new_damaging_ic25
new_damaging_og50
x37348876_m8
```

Run:

```bash
genemasker \
  --annot data/vep.annot.tsv.bgz \
  --stat data/stat.tsv \
  --stat-id-col ID \
  --stat-maf-col MAF \
  --run-masks-file config/masks.txt \
  --out results/run
```

genemasker
**********

This application produces masked group/gene variant inclusion files for rare variant analysis.

|DOI|

.. |DOI| image:: https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21456131-blue.svg
   :target: https://doi.org/10.5281/zenodo.21456131
   :alt: DOI: 10.5281/zenodo.21456131

Installation
************

This application has been tested using the following version of Python and required modules.
   - python=3.9.21
   - numpy=2.0.2
   - pandas=2.2.3
   - scikit-learn=1.6.1
   - scipy=1.13.1
   - psutil=6.1.1
   - dask=2024.8.0

It can be installed in a conda environment via the provided yml file as below.

.. code-block:: bash

   conda env create -f env/environment.yml
   source activate genemasker

With your environment activated, you are ready to install genemasker.

.. code-block:: bash

   python setup.py install

Another option is to run genemasker via Docker (see the included Dockerfile).

Getting Started
***************

.. code-block:: bash

   source activate genemasker

Verify that genemasker is functional using the following command to display help.

.. code-block:: bash

   genemasker --help

Documentation Overview
**********************

``genemasker`` generates Regenie-compatible gene/group mask files from
VEP-style annotation data for rare-variant burden testing. It can run from raw
annotation and frequency data, reuse fitted score models to project onto a new
annotation file, or resume from scored or filtered parquet checkpoints.

For detailed usage, see the documentation in ``docs/``:

- `Mask production workflows <docs/workflows.md>`_ — raw-data, projection,
  scored-parquet, filtered-parquet, and chunked workflows.
- `Custom masks and annotation definitions <docs/custom-masks.md>`_ — creating
  custom Python definition and filter files.
- `Built-in masks <docs/masks.md>`_ — available built-in mask names.
- `CLI reference <docs/cli-reference.md>`_ — every command-line option.

What It Produces
****************

For output prefix ``--out results/study1``, a completed non-chunked run writes:

- ``results/study1.regenie.setlist.tsv``
- ``results/study1.regenie.annotations.<mask>.tsv`` (one file per mask)
- ``results/study1.regenie.mask.<mask>.tsv`` (one file per mask)
- ``results/study1.results.tsv.gz`` (full tabular results)
- ``results/study1.genemasker.log`` (or
  ``results/study1.chunk<N>.genemasker.log`` when ``--chunk`` is used)

It also writes score summaries as applicable: ``rankscore_miss.tsv``,
``rankscore_maf_corr.tsv``, ``pc_maf_corr.tsv``, ``ic_maf_corr.tsv``,
``pca_explained_variance.tsv``, ``combined_score_corr.tsv``, and
``damaging_prop.tsv``. With ``--save-all``, it also retains intermediate parquet
files and writes the fitted imputation, PCA, and ICA pipeline files.

Input Requirements
******************

Annotation file (``--annot``)
=============================

Provide a tab-delimited VEP-style annotation file containing the columns and
dtypes defined in ``genemasker/definitions.py``. Files ending in ``.bgz`` are
read as gzip-compressed files; other compression is inferred by pandas.

The standard definitions include ``#Uploaded_variation``, ``Feature``,
``Feature_type``, ``Gene``, ``PICK``, ``Consequence``, ``IMPACT``, ``DOMAINS``,
``LoF``, predictor/categorical fields, rank-score fields,
``CADD_phred_hg19``, and ``REVEL_score``. By default, only VEP rows where
``PICK == "1"`` are retained. Use ``--include-transcripts`` to group by gene and
transcript instead, or ``--conserved-domains-only`` to retain only annotations
with a non-null ``DOMAINS`` value.

Variant frequency/stat file (``--stat``)
========================================

For a run that calculates scores, provide a tab-delimited stat file with
``--stat``, ``--stat-id-col``, and at least one of ``--stat-maf-col`` or
``--stat-mac-col``. The selected fields are renamed internally to ``MAF`` and
``MAC``. Most built-in masks and the score-correlation steps use MAF, so normal
runs should include ``--stat-maf-col``.

Optional custom extension files
===============================

Use ``--user-definitions`` to load a Python module with custom annotation
definitions, and ``--user-defined-filters`` to load a Python module with
additional mask/filter functions. See `Custom masks and annotation definitions
<docs/custom-masks.md>`_ for the required file structure and examples.

Pipeline Modes
**************

- **Full run from annotation:** reads annotation and stat data, calculates
  damage scores, applies masks, and writes Regenie files.
- **Projection onto new annotation:** reuses the imputation/PCA/ICA pipelines
  and score-correlation artifacts from a prior run while scoring a new
  annotation file.
- **Resume from scored parquet:** use ``--generate-from-scored`` to recompute
  filters and group files without recalculating scores.
- **Resume from filtered parquet:** use ``--generate-from-filtered`` to write
  Regenie files from mask-bearing parquet data without recalculating filters.

Quick Start
***********

.. code-block:: bash

   genemasker \
     --annot data/vep.annot.tsv.bgz \
     --stat data/variant_stats.tsv \
     --stat-id-col variant_id \
     --stat-maf-col maf \
     --stat-mac-col mac \
     --out results/study1

Resume from scored parquet:

.. code-block:: bash

   genemasker \
     --generate-from-scored 'results/study1_tmp/*.scored.parquet' \
     --out results/study1_rerun

Resume from filtered parquet:

.. code-block:: bash

   genemasker \
     --generate-from-filtered 'results/study1_tmp/*.filters.parquet' \
     --out results/study1_rerun

Reusable Model Projection
*************************

To score a new annotation file with a prior run's model, provide ``--annot``,
the stat-file arguments, and all seven projection artifacts. The new annotation
must have compatible rank-score fields and annotation semantics.

.. code-block:: bash

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
     --out results/study2_projected

Run the source analysis with ``--save-all`` to preserve the fitted ``.pkl``
pipeline files and intermediate checkpoint parquet files.

Mask Selection and Output Controls
**********************************

By default, all 18 curated baseline masks are run. Select only the
low-frequency or rare baseline strategy with ``--mask-set``:

.. code-block:: bash

   genemasker \
     --annot data/vep.annot.tsv.bgz \
     --stat data/variant_stats.tsv \
     --stat-id-col variant_id \
     --stat-maf-col maf \
     --mask-set rare \
     --out results/study1_rare

Run an explicit subset, including custom masks, by name:

.. code-block:: bash

   genemasker \
     --annot data/vep.annot.tsv.bgz \
     --stat data/variant_stats.tsv \
     --stat-id-col variant_id \
     --stat-maf-col maf \
     --run-masks new_damaging_og25_0_01,x37348876_m8 \
     --out results/study1_subset

Or provide one mask name per line with ``--run-masks-file config/masks.txt``.
``--mask-set``, ``--run-masks``, and ``--run-masks-file`` are mutually
exclusive; when none is used, all curated baseline masks are calculated.

Use ``--recode-chrs`` to apply regex-based chromosome recoding in final output:

.. code-block:: bash

   genemasker \
     --generate-from-filtered 'results/study1_tmp/*.filters.parquet' \
     --recode-chrs '{"^X$":"23","^Y$":"24","^MT$":"26"}' \
     --out results/study1_recoded

# Equivalence Testing

The main script used to perform equivalence testing is [equiv_test_vectorized.py](https://github.com/juliaolivieri/equivalence_testing/blob/main/scripts/equiv_test_vectorized.py). It assumes input file has been normalized (in this pipeline, it is normalized using DESeq2). 

Current project status and the recent workflow changes are summarized in:

1. [CURRENT_STATE.md](/Users/jolivie1/Desktop/Research/equivalence_testing/CURRENT_STATE.md:1)

## Quick start with test data

The [`test_data`](https://github.com/juliaolivieri/equivalence_testing/tree/main/test_data) folder contains several small example datasets. For a first run, use the "intuitive" dataset:

1. [`intuitive_data.csv`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/intuitive_data.csv)
1. [`intuitive_meta.csv`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/intuitive_meta.csv)
1. [`run_intuitive_test.sh`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/run_intuitive_test.sh)

From the [`test_data`](https://github.com/juliaolivieri/equivalence_testing/tree/main/test_data) folder, run:

`$ bash run_intuitive_test.sh`

This creates the following files:

Intermediary:

1. [`intuitive_condition_deseq2_results.csv`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/intuitive_condition_deseq2_results.csv) (DESeq2 results for input data)
1. [`intuitive_deseq2_normalized_counts.csv`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/intuitive_deseq2_normalized_counts.csv) (DESeq2 normalized matrix used in the equivalence step)

Result:
1. [`intuitive_condition_1.0_results.csv`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/intuitive_condition_1.0_results.csv) (Output file with differential, equivalence, and low-information calls)
1. [`intuitive_condition_1.0_deseq2_pvals.png`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/intuitive_condition_1.0_deseq2_pvals.png) (Plot showing correlation between DESeq2 difference p values and current difference p values)
1. [`intuitive_condition_1.0_deseq2_foldchange.png`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/intuitive_condition_1.0_deseq2_foldchange.png) (Plot showing correlation between DESeq2 fold change and current fold change)
1. [`intuitive_condition_1.0_readdepth.png`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/intuitive_condition_1.0_readdepth.png) (Plot showing how the fraction of different, equivalent, inconclusive, and low-information genes changes with sequencing depth)
1. [`intuitive_condition_1.0_volcano.png`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/intuitive_condition_1.0_volcano.png) (Volcano plot of results colored by category)

This dataset is meant to be easier to interpret by eye than the original toy example. It keeps sample totals fairly similar and includes several stable genes, so DESeq2 normalization changes the rows only modestly. That makes it easier to see how the final categories arise after normalization and log transformation.

With the current defaults (`delta = 1`, `low_count_threshold = 5`, `min_nonzero_per_group = 2`), this dataset should show:

1. one gene called `different`
1. one gene called `inconclusive`
1. one gene called `low_information`
1. several genes called `equivalent`

The older files ([`data.csv`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/data.csv), [`meta.csv`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/meta.csv), and [`run_test.sh`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/run_test.sh)) are still available, but they are less intuitive for a first pass because normalization and `log2(x+1)` distort the small-count rows more strongly.

## Input arguments

1. `--savename`: Name to save output as
1. `--outpath`: Folder to save output to
1. `--infile`: .csv file with rows corresponding to genes and columns corresponding to samples.
1. `--meta`: .csv file with rows corresponding to samples and columns corresponding to metadata about samples. All sample names in the first column should match column names in `infile` (though `meta` can contain a subset of the samples in `infile`). `condition` must be a column in the file that has two categories.  
1. `--condition`: Column of the `meta` file that will be used to split the samples into two groups. The script requires exactly two groups, each with at least two samples.
1. `--delta`: Input a quantitative value that defines the practical-effect boundary on the log2 scale. Default: `delta = 1`. This corresponds to a 2x fold change. The differential test asks whether the absolute group difference is significantly greater than `delta`, and the equivalence test asks whether the absolute group difference is significantly less than `delta`.
1. `--low_count_threshold`: Default `5`. A gene is flagged as low-information if both group means in the normalized pre-log matrix are below this threshold.
1. `--min_nonzero_per_group`: Default `2`. A gene is flagged as low-information if either group has fewer than this many nonzero samples.

## Batch analysis workflow

For GEO-scale runs, the preferred workflow is now manifest-driven.

The curated screening table lives at:

1. [`download_geo_datasets/geo_bulk_rnaseq_screening.csv`](https://github.com/juliaolivieri/equivalence_testing/blob/main/download_geo_datasets/geo_bulk_rnaseq_screening.csv)

An analysis-ready manifest lives at:

1. [`download_geo_datasets/analysis_manifest.csv`](https://github.com/juliaolivieri/equivalence_testing/blob/main/download_geo_datasets/analysis_manifest.csv)

The manifest is one row per planned analysis. It starts from the eligible screened GSEs and adds the local information needed to actually run the pipeline:

1. local count-matrix path
1. local metadata path
1. condition column
1. control and treatment labels
1. whether the row is selected for analysis
1. whether the row is analysis-ready

Two helper scripts support this workflow:

1. [`scripts/bootstrap_analysis_manifest.py`](https://github.com/juliaolivieri/equivalence_testing/blob/main/scripts/bootstrap_analysis_manifest.py): build or refresh an analysis manifest from the GEO screening CSV
1. [`scripts/run_batch_equiv.py`](https://github.com/juliaolivieri/equivalence_testing/blob/main/scripts/run_batch_equiv.py): run all selected and analysis-ready rows from the manifest
1. [`scripts/prepare_geo_counts.py`](https://github.com/juliaolivieri/equivalence_testing/blob/main/scripts/prepare_geo_counts.py): extract a counts-only matrix and matching metadata table from a GEO supplement plus a simple sample-condition sheet

Typical workflow:

1. Start from [`download_geo_datasets/analysis_manifest.csv`](https://github.com/juliaolivieri/equivalence_testing/blob/main/download_geo_datasets/analysis_manifest.csv)
1. For each GSE you want to run, prepare a counts-only matrix and metadata table, then fill in `local_counts_path`, `local_meta_path`, `condition_column`, and any comparison notes
1. Set `selected_for_analysis = True` and `analysis_ready = True`
1. Run the batch script

For simple GEO count files, preparation can be done with:

```bash
python3 scripts/prepare_geo_counts.py \
  --input <geo_count_matrix> \
  --sample-sheet <sample_condition_csv> \
  --output-counts <prepared_counts.csv> \
  --output-meta <prepared_meta.csv> \
  --gene-column gene_id
```

Example:

```bash
python3 scripts/run_batch_equiv.py \
  --manifest download_geo_datasets/analysis_manifest.csv \
  --runs-dir runs/batch_equiv
```

The runner writes one folder per analysis under `runs/batch_equiv/<analysis_id>/` and a combined dataset-level summary table at:

1. `runs/batch_equiv/batch_summary.csv`

The first real GEO example already prepared in this workflow is `GSE327455`. Its current outputs live in:

1. [runs/batch_equiv/GSE327455](/Users/jolivie1/Desktop/Research/equivalence_testing/runs/batch_equiv/GSE327455:1)
1. [runs/batch_equiv/batch_summary.csv](/Users/jolivie1/Desktop/Research/equivalence_testing/runs/batch_equiv/batch_summary.csv:1)

The manifest includes a local smoke-test row, `TEST_INTUITIVE`, that points to the small intuitive toy dataset. This is useful for validating the batch pipeline before preparing real GEO datasets.

If your Python environment is not the default one on your machine, pass the explicit interpreter:

```bash
python3 scripts/run_batch_equiv.py \
  --manifest download_geo_datasets/analysis_manifest.csv \
  --runs-dir runs/batch_equiv \
  --python-executable /opt/anaconda3/bin/python3
```

The older [`data/data_files.csv`](https://github.com/juliaolivieri/equivalence_testing/blob/main/data/data_files.csv) and [`scripts/write_submission_script.py`](https://github.com/juliaolivieri/equivalence_testing/blob/main/scripts/write_submission_script.py) workflow is still present for historical runs, but it is no longer the recommended path for multi-dataset GEO analysis.

## Output

1. `<savename>_<condition>_<delta>_results.csv`: contains one row per gene. Includes adjusted differential and equivalence p values, raw test calls, downstream reporting category, and a low-information flag.

## Statistical interpretation

For each gene, the pipeline runs two complementary Welch-style tests on the log2-transformed normalized expression values:

1. Differential test: `H0: |mu1 - mu2| <= delta` vs `H1: |mu1 - mu2| > delta`
1. Equivalence test: `H0: |mu1 - mu2| >= delta` vs `H1: |mu1 - mu2| < delta`

This means `delta` is built into both tests directly, rather than applying a standard difference test followed by an effect-size cutoff.

## Output categories

The result file reports two related category columns:

1. `test_category`: the direct statistical call from the two tests (`different`, `equivalent`, or `inconclusive`)
1. `category`: the downstream reporting category

If `low_information = True`, then `category` is set to `low_information` even if `test_category` is `equivalent` or `inconclusive`. This keeps low-count genes out of the downstream `different` / `equivalent` / `inconclusive` totals while preserving the underlying test result.

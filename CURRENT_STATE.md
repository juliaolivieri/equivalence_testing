# Current State

This file records the current project state after the recent cleanup of the statistical workflow, toy datasets, and GEO batch-analysis setup.

## Statistical workflow changes

The current vectorized pipeline is implemented in:

1. [scripts/equiv_test_funcs_vectorized.py](/Users/jolivie1/Desktop/Research/equivalence_testing/scripts/equiv_test_funcs_vectorized.py:1)
1. [scripts/equiv_test_vectorized.py](/Users/jolivie1/Desktop/Research/equivalence_testing/scripts/equiv_test_vectorized.py:1)

Recent changes:

1. The differential p value was corrected so it tests evidence for `|effect| > delta`, rather than evidence for any nonzero effect.
1. A separate low-information flag was added using:
   `low_count_threshold` (default `5`) and `min_nonzero_per_group` (default `2`).
1. The output now distinguishes:
   `test_category` = raw statistical call and
   `category` = downstream reporting category.
1. If `low_information = True`, then `category` becomes `low_information` so downstream totals satisfy:
   `different + equivalent + inconclusive + low_information = total genes`.

## Toy datasets

Three small datasets are currently useful:

1. [test_data/data.csv](/Users/jolivie1/Desktop/Research/equivalence_testing/test_data/data.csv:1)
   Original toy dataset. Kept for reference, but less intuitive after normalization and `log2(x+1)`.
1. [test_data/showcase_data.csv](/Users/jolivie1/Desktop/Research/equivalence_testing/test_data/showcase_data.csv:1)
   Small designed dataset that spans all downstream buckets.
1. [test_data/intuitive_data.csv](/Users/jolivie1/Desktop/Research/equivalence_testing/test_data/intuitive_data.csv:1)
   Recommended first test dataset. Sample totals are fairly similar and normalization changes the rows only modestly.

Associated scripts:

1. [test_data/run_test.sh](/Users/jolivie1/Desktop/Research/equivalence_testing/test_data/run_test.sh:1)
1. [test_data/run_showcase_test.sh](/Users/jolivie1/Desktop/Research/equivalence_testing/test_data/run_showcase_test.sh:1)
1. [test_data/run_intuitive_test.sh](/Users/jolivie1/Desktop/Research/equivalence_testing/test_data/run_intuitive_test.sh:1)

## GEO screening and manifest workflow

The GEO screening handoff lives in:

1. [download_geo_datasets/geo_bulk_rnaseq_screening.csv](/Users/jolivie1/Desktop/Research/equivalence_testing/download_geo_datasets/geo_bulk_rnaseq_screening.csv:1)

This file is the source registry for real-data work. It currently contains eligible, candidate, and exclude rows from the earlier GEO screening pass.

The analysis manifest lives in:

1. [download_geo_datasets/analysis_manifest.csv](/Users/jolivie1/Desktop/Research/equivalence_testing/download_geo_datasets/analysis_manifest.csv:1)

This file is one row per planned analysis, not just one row per screened GEO series. It stores:

1. local counts path
1. local metadata path
1. condition column
1. control/treatment labels
1. whether the row is selected for analysis
1. whether the row is analysis-ready
1. current statistical settings (`delta`, low-information thresholds)

Helper scripts:

1. [scripts/bootstrap_analysis_manifest.py](/Users/jolivie1/Desktop/Research/equivalence_testing/scripts/bootstrap_analysis_manifest.py:1)
   Build or refresh the analysis manifest from the screening CSV.
1. [scripts/prepare_geo_counts.py](/Users/jolivie1/Desktop/Research/equivalence_testing/scripts/prepare_geo_counts.py:1)
   Reduce a GEO supplement to the counts-only matrix used by the runner and write a matching metadata file from a simple `sample,condition` sheet.
1. [scripts/run_batch_equiv.py](/Users/jolivie1/Desktop/Research/equivalence_testing/scripts/run_batch_equiv.py:1)
   Run all manifest rows where `selected_for_analysis = True` and `analysis_ready = True`.

## First real GEO run

The first real eligible GEO dataset prepared and run is:

1. `GSE327455`

Prepared input files:

1. [test_data/prepared_geo/GSE327455/counts.csv](/Users/jolivie1/Desktop/Research/equivalence_testing/test_data/prepared_geo/GSE327455/counts.csv:1)
1. [test_data/prepared_geo/GSE327455/meta.csv](/Users/jolivie1/Desktop/Research/equivalence_testing/test_data/prepared_geo/GSE327455/meta.csv:1)
1. [download_geo_datasets/prep_configs/GSE327455_sample_conditions.csv](/Users/jolivie1/Desktop/Research/equivalence_testing/download_geo_datasets/prep_configs/GSE327455_sample_conditions.csv:1)

The manifest row for `GSE327455` in [download_geo_datasets/analysis_manifest.csv](/Users/jolivie1/Desktop/Research/equivalence_testing/download_geo_datasets/analysis_manifest.csv:1) is now filled in and marked analysis-ready.

Current output location for the real-data run:

1. [runs/batch_equiv/GSE327455](/Users/jolivie1/Desktop/Research/equivalence_testing/runs/batch_equiv/GSE327455:1)
1. [runs/batch_equiv/batch_summary.csv](/Users/jolivie1/Desktop/Research/equivalence_testing/runs/batch_equiv/batch_summary.csv:1)

Main per-gene result table:

1. [runs/batch_equiv/GSE327455/GSE327455_condition_1.0_results.csv](/Users/jolivie1/Desktop/Research/equivalence_testing/runs/batch_equiv/GSE327455/GSE327455_condition_1.0_results.csv:1)

Dataset-level outcome for `GSE327455`:

1. genes tested: `23152`
1. different: `173`
1. equivalent: `7634`
1. inconclusive: `8227`
1. low_information: `7118`

## What the batch summary currently records

The combined batch summary currently includes:

1. `analysis_id`
1. `gse_accession`
1. `description`
1. `species`
1. `comparison_name`
1. `comparison_description`
1. `condition_column`
1. `n_group1`
1. `n_group2`
1. `total_n`
1. `num_genes_tested`
1. `num_different`
1. `num_equivalent`
1. `num_inconclusive`
1. `num_low_information`
1. category fractions
1. `mean_library_size`
1. `median_library_size`
1. `delta`
1. low-information thresholds
1. `status`
1. `error_message`
1. `result_path`

So replicate counts and read-depth summaries are already being tracked at the dataset level.

## Known cleanup items

1. The `result_path` value inside the existing [runs/batch_equiv/batch_summary.csv](/Users/jolivie1/Desktop/Research/equivalence_testing/runs/batch_equiv/batch_summary.csv:1) row for `GSE327455` still reflects the pre-move location from before the run outputs were moved into `runs/batch_equiv`.
1. Some temporary validation files remain in `test_data/`:
   `analysis_manifest_test.csv`,
   `analysis_manifest_gse327455_test.csv`,
   `batch_runs/`,
   `batch_runs_real/`,
   and the downloaded raw GEO supplement `GSE327455_gene_count.txt.gz`.
1. Prepared GEO inputs are currently staged under `test_data/prepared_geo/` because this environment allowed writes there. If desired, they can be moved later to a more permanent location once write permissions are settled.

## Recommended next step

Repeat the same preparation pattern for the next highest-confidence eligible GEO series:

1. create a `sample,condition` file under `download_geo_datasets/prep_configs/`
1. run [scripts/prepare_geo_counts.py](/Users/jolivie1/Desktop/Research/equivalence_testing/scripts/prepare_geo_counts.py:1)
1. fill the row in [download_geo_datasets/analysis_manifest.csv](/Users/jolivie1/Desktop/Research/equivalence_testing/download_geo_datasets/analysis_manifest.csv:1)
1. run [scripts/run_batch_equiv.py](/Users/jolivie1/Desktop/Research/equivalence_testing/scripts/run_batch_equiv.py:1)

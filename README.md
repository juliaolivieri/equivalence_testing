# Equivalence Testing

The main script used to perform equivalence testing is [equiv_test_vectorized.py](https://github.com/juliaolivieri/equivalence_testing/blob/main/scripts/equiv_test_vectorized.py). It assumes input file has been normalized (in this pipeline, it is normalized using DESeq2). 

## Quick start with test data

The [`test_data`](https://github.com/juliaolivieri/equivalence_testing/tree/main/test_data) contains example input files ([`data.csv`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/data.csv) and [`meta.csv`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/meta.csv)), and an example submission script [`run_test.sh`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/run_test.sh). 

From the [`test_data`](https://github.com/juliaolivieri/equivalence_testing/tree/main/test_data) folder, run the script:

`$ bash run_test.sh`

This should create the following output files:

Intermediary:

1. [`test_condition_deseq2_results.csv`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/test_condition_deseq2_results.csv) (DESeq2 results for input data)
1. [`test_deseq2_normalized_counts.csv`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/test_deseq2_normalized_counts.csv) (DESeq2 normalized matrix to be used in next step)

Result:
1. [`test_condition_1.0_results.csv`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/test_condition_1.0_results.csv) (Output file with different and equivalent p values)
1. [`test_condition_1.0_deseq2_pvals.png`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/test_condition_1.0_deseq2_pvals.png) (Plot showing correlation between DESeq2 difference p values and current difference p values)
1. [`test_condition_1.0_deseq2_foldchange.png`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/test_condition_1.0_deseq2_foldchange.png) (Plot showing correlation between DESeq2 fold change and current fold change)
1. [`test_condition_1.0_readdepth.png`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/test_condition_1.0_readdepth.png) (Plot showing how the fraction of different, equivalent, and inconclusive genes changes with sequencing depth)
1. [`test_condition_1.0_volcano.png`](https://github.com/juliaolivieri/equivalence_testing/blob/main/test_data/test_condition_1.0_volcano.png) (Volcano plot of results colored by different, equivalent, inconclusive)

This test data has six samples and four genes. The results should show that two genes are significantly equivalent between groups, one is significantly different, and one is inconclusive.

## Input arguments

1. `--savename`: Name to save output as
1. `--outpath`: Folder to save output to
1. `--infile`: .csv file with rows corresponding to genes and columns corresponding to samples.
1. `--meta`: .csv file with rows corresponding to samples and columns corresponding to metadata about samples. All sample names in the first column should match column names in `infile` (though `meta` can contain a subset of the samples in `infile`). `condition` must be a column in the file that has two categories.  
1. `--condition`: Column of the `meta` file that will be used to split the samples into two groups.
1. `--delta`: Input a quantitative value to use as the cutoff for equivalence testing (so, the p value corresponds to the two groups having log2 average values within delta of each other). Default: `delta = 1`. This corresponds to 2x fold change.

## Running instructions

The most updated way of running is to put the details of the dataset you want to run in [`data/data_files.csv`](https://github.com/juliaolivieri/equivalence_testing/blob/main/data/data_files.csv).

Then, edit [`scripts/write_submission_script.py`](https://github.com/juliaolivieri/equivalence_testing/blob/main/scripts/write_submission_script.py) to set `identifier` equal to the id you want to run for. Run this code:

`$ python write_submission_script.py`

This will create the file [`scripts/run_auto.sh`](https://github.com/juliaolivieri/equivalence_testing/blob/main/scripts/run_auto.sh). This will run DESeq2 and the equivalence testing script. Run it as follows:

`$ bash run_auto.sh`

Output will be created in `/Users/jolivie1/Desktop/Research/equivalence_testing_output/scripts/output/equiv_test_vectorized/`.

## Output

1. `<savename>_<condition>_<delta>_results.csv`: contains one row per gene. For each gene, includes information such as the adjusted difference and adjusted equivalence p value.






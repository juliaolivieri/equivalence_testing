# Equivalence Testing

## Input arguments

1. `--savename`: Name to save output as
1. `--outpath`: Folder to save output to
1. `--infile`: .csv file with rows corresponding to samples and columns corresponding to genes. The last column should be labeled `cell_type`, and assign each sample to the corresponding group. 
1. `--delta`: Input a quantitative value to use as the cutoff for equivalence testing (so, the p value corresponds to the two groups having average values within delta of each other)
1. `--log_scale`: Including this flag log-scales the values before equivalence testing
1. `--normalize`: Including this flag normalizes the input by row count

## Running instructions

You can run the script `run_perform_test.sh` after replacing any command-line arguments you want to change:

```
bash run_perform_test.sh
```

Or you can run the command itself. For example:

```
python perform_test_tabular_data.py \
        --savename monocyte_raw \
        --outpath results/ \
        --infile HLCA_normal_1000_classical-monocyte_non-classical-monocyte_1000_100_raw.csv \
        --delta 1 \
        --log_scale

```



## Output

1. `<savename>_results.csv`: contains one row per gene. For each gene, includes information such as the adjusted difference and adjusted equivalence p value.
1. `<savename>_<delta>_comparison.png`: Plot of diff vs equiv p value
1. `<savename>_<delta>_equiv_pval.png`: histogram of equiv p value
1. `<savename>_<delta>_diff_pval.png`: histogram of diff p value





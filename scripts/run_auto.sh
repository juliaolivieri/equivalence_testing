#!/usr/bin/env bash

date

python equiv_test_vectorized.py --savename HLCA_column_alph_1000 --outpath ~/Desktop/Research/equivalence_testing_output/scripts/output/equiv_test_vectorized/ --infile ~/Desktop/Research/equivalence_testing_output/scripts/output/deseq_normalization/HLCA_column_alph_1000_deseq2_normalized_counts.csv  --meta ~/Desktop/Research/equivalence_testing_output/notebooks/output/reformat_HLCA/HLCA_normal_1000_ciliated-columnar-cell-of-tracheobronchial-tree_CD8-positive-alpha-beta-T-cell_all_all_raw_2000cells_meta.csv --condition cell_type --delta 1 --deseq2_results ~/Desktop/Research/equivalence_testing_output/scripts/output/deseq_normalization/HLCA_column_alph_1000_cell_type_deseq2_results.csv

date
echo 'complete'


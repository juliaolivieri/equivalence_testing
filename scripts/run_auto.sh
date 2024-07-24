#!/usr/bin/env bash

date

Rscript deseq_normalization.R ~/Desktop/Research/equivalence_testing_output/notebooks/output/reformat_HLCA/HLCA_normal_1000_classical-monocyte_non-classical-monocyte_all_all_raw_counts.csv ~/Desktop/Research/equivalence_testing_output/notebooks/output/reformat_HLCA/HLCA_normal_1000_classical-monocyte_non-classical-monocyte_all_all_raw_meta.csv cell_type ~/Desktop/Research/equivalence_testing_output/scripts/output/deseq_normalization/HLCA_monocytes

date
echo 'normalized'

python perform_test_tabular_data.py --savename HLCA_monocytes_cell_type_1_0.1 --outpath ~/Desktop/Research/equivalence_testing_output/scripts/output/perform_test_tabular_data/ --infile ~/Desktop/Research/equivalence_testing_output/scripts/output/deseq_normalization/HLCA_monocytes_deseq2_normalized_counts.csv --delta 1 --epsilon 0.1 --log_scale --meta ~/Desktop/Research/equivalence_testing_output/notebooks/output/reformat_HLCA/HLCA_normal_1000_classical-monocyte_non-classical-monocyte_all_all_raw_meta.csv --condition cell_type

date
echo 'complete'


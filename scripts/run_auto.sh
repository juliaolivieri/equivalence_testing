#!/usr/bin/env bash

date

Rscript deseq_normalization.R ~/Desktop/Research/Datasets/GSE206932_merged.counts.bulk.csv ~/Desktop/Research/Datasets/GSE206932_meta.csv age ~/Desktop/Research/equivalence_testing_output/scripts/output/deseq_normalization/GSE206932

date
echo 'normalized'

python perform_test_tabular_data.py --savename GSE206932_age_1_0.1 --outpath ~/Desktop/Research/equivalence_testing_output/scripts/output/perform_test_tabular_data/ --infile ~/Desktop/Research/equivalence_testing_output/scripts/output/deseq_normalization/GSE206932_deseq2_normalized_counts.csv --delta 1 --epsilon 0.1 --log_scale --meta ~/Desktop/Research/Datasets/GSE206932_meta.csv --condition age

date
echo 'complete'


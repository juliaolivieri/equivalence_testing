#!/usr/bin/env bash

date

Rscript deseq_normalization.R ~/Desktop/Research/Datasets/GSE117891_all_6148.umi.count.matrix.csv ~/Desktop/Research/Datasets/GSE117891_meta_twocat.csv sample_type ~/Desktop/Research/equivalence_testing_output/scripts/output/deseq_normalization/GSE117891

date
echo 'normalized'

python perform_test_tabular_data.py --savename GSE117891_sample_type_1 --outpath ~/Desktop/Research/equivalence_testing_output/scripts/output/perform_test_tabular_data/ --infile ~/Desktop/Research/Datasets/GSE117891_all_6148.umi.count.matrix.csv --delta 1 --log_scale --meta ~/Desktop/Research/Datasets/GSE117891_meta_twocat.csv --condition sample_type

date
echo 'complete'


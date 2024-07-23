#!/usr/bin/env bash

date

Rscript deseq_normalization.R ~/Desktop/Research/Datasets/GSE268000_Kallisto_raw.counts.csv ~/Desktop/Research/Datasets/GSE268000_meta_twocat.csv treatment ~/Desktop/Research/equivalence_testing_output/scripts/output/deseq_normalization/GSE268000

date
echo 'normalized'

python perform_test_tabular_data.py --savename GSE268000_treatment_1 --outpath ~/Desktop/Research/equivalence_testing_output/scripts/output/perform_test_tabular_data/ --infile ~/Desktop/Research/Datasets/GSE268000_Kallisto_raw.counts.csv --delta 1 --log_scale --meta ~/Desktop/Research/Datasets/GSE268000_meta_twocat.csv --condition treatment

date
echo 'complete'


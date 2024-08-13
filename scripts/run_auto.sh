#!/usr/bin/env bash

date

Rscript deseq_normalization.R ~/Desktop/Research/Datasets/GSE206932_merged.counts.bulk.csv ~/Desktop/Research/Datasets/GSE206932_meta.csv age ~/Desktop/Research/equivalence_testing_output/scripts/output/deseq_normalization/GSE206932_age

date
echo 'normalized'

python equiv_test_vectorized.py --savename GSE206932_age --outpath ~/Desktop/Research/equivalence_testing_output/scripts/output/equiv_test_vectorized/ --infile ~/Desktop/Research/equivalence_testing_output/scripts/output/deseq_normalization/GSE206932_age_deseq2_normalized_counts.csv  --meta ~/Desktop/Research/Datasets/GSE206932_meta.csv --condition age --delta 1

date
echo 'complete'


#!/usr/bin/env bash

date

Rscript deseq_normalization.R /Users/jolivie1/Desktop/Research/oncogenetic_fusions/fc_output/all.csv /Users/jolivie1/Desktop/Research/oncogenetic_fusions/metadata/full_meta.csv condition output/deseq_normalization/onco

date
echo 'normalized'

python equiv_test_vectorized.py --savename onco --outpath output/equiv_test_vectorized/ --infile output/deseq_normalization/onco_deseq2_normalized_counts.csv  --meta /Users/jolivie1/Desktop/Research/oncogenetic_fusions/metadata/full_meta.csv --condition condition --delta 1 --deseq2_results output/deseq_normalization/onco_condition_deseq2_results.csv

date
echo 'complete'


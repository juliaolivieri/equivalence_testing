#!/usr/bin/env bash

python perform_test_tabular_data.py \
	--savename HLCA_normal_1000_CD4-positive-alpha-beta-T-cell_CD8-positive-alpha-beta-T-cell \
	--outpath /exports/home/jolivieri/equivalence_testing_output/scripts/output/perform_test_tabular_data \
	--infile /exports/home/jolivieri/equivalence_testing_output/notebooks/output/sample_tabular_data/HLCA_normal_1000_CD4-positive-alpha-beta-T-cell_CD8-positive-alpha-beta-T-cell_1000_100.csv \
	--delta 0.05

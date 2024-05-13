#!/usr/bin/env bash

#SAVENAME="glioma_matrix_edited_2types"
#INFILE="/exports/home/jolivieri/data/test/${SAVENAME}.csv"

#SAVENAME="glioma_test"
#INFILE="/exports/home/jolivieri/data/test/${SAVENAME}.csv"


# SAVENAME=HLCA_normal_1000_CD4-positive-alpha-beta-T-cell_CD8-positive-alpha-beta-T-cell
# INFILE=/exports/home/jolivieri/equivalence_testing_output/notebooks/output/sample_tabular_data/${SAVENAME}_1000_100.csv

SAVENAME="monocyte_raw"
INFILE="/exports/home/jolivieri/equivalence_testing_output/scripts/output/sample_tabular_data/HLCA_normal_1000_classical-monocyte_non-classical-monocyte_1000_100_raw.csv"

OUTPATH=/exports/home/jolivieri/equivalence_testing_output/scripts/output/perform_test_tabular_data/

python perform_test_tabular_data.py \
	--savename $SAVENAME \
	--outpath $OUTPATH \
	--infile $INFILE \
	--delta 1 \
        --log_scale 

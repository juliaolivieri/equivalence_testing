#!/usr/bin/env bash


#DATANAME="GSE134296"
DATANAME="Snf2_WT_6"
#INFILE=/Users/jolivie1/Desktop/Research/equivalence_testing/notebooks/output/Snf2_filtering/.csv
INFILE="/Users/jolivie1/Desktop/Research/equivalence_testing/notebooks/output/Snf2_filtering/${DATANAME}.csv"
#INFILE="~/Desktop/Research/Datasets/GSE206932_merged.counts.bulk.csv"
#INFILE="~/Desktop/Research/DESeq2/GSE206932_experimental_condition_normalized_counts.csv"
#META="~/Desktop/Research/Datasets/GSE206932_meta.csv"
META="/Users/jolivie1/Desktop/Research/equivalence_testing/notebooks/output/Snf2_filtering/meta_${DATANAME}.csv"
#META=/Users/jolivie1/Desktop/Research/equivalence_testing/notebooks/output/Snf2_filtering/meta_Snf2_WT_40.csv
METACOLUMN="condition"

OUTPATH_DESEQ="/Users/jolivie1/Desktop/Research/equivalence_testing_output/scripts/output/deseq_normalization/${DATANAME}"



date
Rscript deseq_normalization.R \
	${INFILE} \
	${META} \
	${METACOLUMN} \
	${OUTPATH_DESEQ}

date
echo "normalized"

DELTA=1
SAVENAME="${DATANAME}_${METACOLUMN}_${DELTA}"
INFILE_EQUIV="${OUTPATH_DESEQ}_deseq2_normalized_counts.csv"

OUTPATH_EQUIV=~/Desktop/Research/equivalence_testing_output/scripts/output/perform_test_tabular_data/


python perform_test_tabular_data.py \
	--savename $SAVENAME \
	--outpath $OUTPATH_EQUIV \
	--infile $INFILE_EQUIV \
	--delta $DELTA \
        --log_scale \
	--meta ${META} \
	--condition ${METACOLUMN}

date
echo "complete"

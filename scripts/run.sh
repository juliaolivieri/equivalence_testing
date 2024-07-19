#!/usr/bin/env bash


DATANAME="GSE206932"
INFILE="~/Desktop/Research/Datasets/GSE206932_merged.counts.bulk.csv"
#INFILE="~/Desktop/Research/DESeq2/GSE206932_experimental_condition_normalized_counts.csv"
META="~/Desktop/Research/Datasets/GSE206932_meta.csv"
METACOLUMN="age"

OUTPATH_DESEQ="~/Desktop/Research/equivalence_testing_output/scripts/output/deseq_normalization/${DATANAME}"

DELTA=1
SAVENAME="${DATANAME}_${METACOLUMN}_${DELTA}"
INFILE="${OUTPATH_DESEQ}_deseq2_normalized_counts.csv"

OUTPATH_EQUIV=~/Desktop/Research/equivalence_testing_output/scripts/output/perform_test_tabular_data/

date
Rscript deseq_normalization.R \
	${INFILE} \
	${META} \
	${METACOLUMN} \
	${OUTPATH_DESEQ}

date
echo "normalized"

python perform_test_tabular_data.py \
	--savename $SAVENAME \
	--outpath $OUTPATH_EQUIV \
	--infile $INFILE \
	--delta $DELTA \
        --log_scale \
	--meta ${META} \
	--condition ${METACOLUMN}

date
echo "complete"

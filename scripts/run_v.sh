#!/usr/bin/env bash


SAVENAME="GSE206932"
INFILE="~/Desktop/Research/Datasets/GSE206932_merged.counts.bulk.csv"
META="~/Desktop/Research/Datasets/GSE206932_meta.csv"
METACOLUMN="age"

OUTPATH_DESEQ="~/Desktop/Research/equivalence_testing_output/scripts/output/deseq_normalization/${SAVENAME}"

DELTA=1
INFILE="${OUTPATH_DESEQ}_deseq2_normalized_counts.csv"

OUTPATH_EQUIV=~/Desktop/Research/equivalence_testing_output/scripts/output/equiv_test_vectorized/

date
python equiv_test_vectorized.py \
	--savename $SAVENAME \
	--outpath $OUTPATH_EQUIV \
	--infile $INFILE \
	--delta $DELTA \
	--meta ${META} \
	--condition ${METACOLUMN}

date
echo "complete"

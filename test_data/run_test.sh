#!/usr/bin/env bash


SAVENAME="test"
INFILE="data.csv"
META="meta.csv"
METACOLUMN="condition"

OUTPATH_DESEQ="${SAVENAME}"

DELTA=1
INFILE_EQUIV="${OUTPATH_DESEQ}_deseq2_normalized_counts.csv"

OUTPATH_EQUIV="./"

date
a="Rscript ../scripts/deseq_normalization.R \
	${INFILE} \
	${META} \
	${METACOLUMN} \
	${OUTPATH_DESEQ}"
echo $a
eval $a

date
echo "normalized"

b="python ../scripts/equiv_test_vectorized.py \
	--savename $SAVENAME \
	--outpath $OUTPATH_EQUIV \
	--infile $INFILE_EQUIV \
	--delta $DELTA \
	--meta ${META} \
	--condition ${METACOLUMN}"

echo $b
eval $b

date
echo "complete"

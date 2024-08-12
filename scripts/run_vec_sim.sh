#!/usr/bin/env bash


SAVENAME="simulation"
#INFILE="/Users/jolivie1/Desktop/Research/equivalence_testing_output/notebooks/output/simulation_vectorize/data.csv"
META="/Users/jolivie1/Desktop/Research/equivalence_testing_output/notebooks/output/simulation_vectorize/meta.csv"
METACOLUMN="condition"

#OUTPATH_DESEQ="~/Desktop/Research/equivalence_testing_output/scripts/output/deseq_normalization/${SAVENAME}"

DELTA=1
INFILE_EQUIV="/Users/jolivie1/Desktop/Research/equivalence_testing_output/notebooks/output/simulation_vectorize/data.csv"

OUTPATH_EQUIV="~/Desktop/Research/equivalence_testing_output/scripts/output/equiv_test_vectorized/"

date
a="Rscript deseq_normalization.R \
	${INFILE} \
	${META} \
	${METACOLUMN} \
	${OUTPATH_DESEQ}"
#echo $a
#eval $a

date
#echo "normalized"

b="python equiv_test_vectorized.py \
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

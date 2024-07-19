#!/usr/bin/env bash


DATANAME="GSE206932"
INFILE="~/Desktop/Research/Datasets/GSE206932_merged.counts.bulk.csv"
#INFILE="~/Desktop/Research/DESeq2/GSE206932_experimental_condition_normalized_counts.csv"
META="~/Desktop/Research/Datasets/GSE206932_meta.csv"
METACOLUMN="experimental_condition"

OUTPATH="~/Desktop/Research/equivalence_testing_output/scripts/output/deseq_normalization/${DATANAME}"

Rscript deseq_normalization.R \
	${INFILE} \
	${META} \
	${METACOLUMN} \
	${OUTPATH}


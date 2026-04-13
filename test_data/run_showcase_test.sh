#!/usr/bin/env bash

SAVENAME="showcase"
INFILE="showcase_data.csv"
META="showcase_meta.csv"
METACOLUMN="condition"

OUTPATH_DESEQ="${SAVENAME}"
DELTA=1
INFILE_EQUIV="${OUTPATH_DESEQ}_deseq2_normalized_counts.csv"
OUTPATH_EQUIV="./"

Rscript ../scripts/deseq_normalization.R \
  "${INFILE}" \
  "${META}" \
  "${METACOLUMN}" \
  "${OUTPATH_DESEQ}"

python3 ../scripts/equiv_test_vectorized.py \
  --savename "${SAVENAME}" \
  --outpath "${OUTPATH_EQUIV}" \
  --infile "${INFILE_EQUIV}" \
  --delta "${DELTA}" \
  --meta "${META}" \
  --condition "${METACOLUMN}" \
  --deseq2_results "${SAVENAME}_${METACOLUMN}_deseq2_results.csv"

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.stats.multitest as multitest
import tqdm
from equiv_test_funcs import *

args = get_args()

df = pd.read_csv(args.infile, index_col = 0).T

meta = pd.read_csv(args.meta, index_col = 0)

# only include samples that are in the metadata
df = df[df.index.isin(meta.index)]



df["condition"] = df.index.map({k : v for k, v in zip(meta.index, meta[args.condition])})

print("SUBSET DF")
print(df.head)

if args.normalize:
  # normalize each row (sample) by the number of counts
  df = normalize(df)

if args.log_scale:
  df = log2_scale(df)

plot = False

# Remove columns that are all zero
nnz = (~(df == 0).all())
df = df[nnz[nnz].index]

# testing out function you can use in Python code
# perform_tests_for_df(df, args.delta)

out_df = loop_over_genes(df, args.delta, plot)
out_df = process_out_df(out_df, args.delta)

plot_results(out_df, args.outpath, args.savename, args.delta)

print(args.savename)
print("delta: {}".format(args.delta))
print(pd.crosstab(out_df['sig_equiv'], out_df['sig_diff']))

summarize_df_significance(out_df, args.delta, verbose = True)

out_df.sort_values("diff_pval_adj", inplace=True)
out_df.to_csv("{}{}_results.csv".format(args.outpath, args.savename),index = False)


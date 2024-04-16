import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.stats.multitest as multitest
import tqdm
from equiv_test_funcs import *

args = get_args()

df = pd.read_csv(args.infile, index_col = 0)

plot = False

out_df = loop_over_genes(df, args.delta, plot)
out_df = process_out_df(out_df)

plot_results(out_df, args.outpath, args.savename, args.delta)

print(args.savename)
print("delta: {}".format(args.delta))
print(pd.crosstab(out_df['sig_equiv'], out_df['sig_diff']))

out_df.sort_values("diff_pval_adj", inplace=True)
out_df.to_csv("{}{}_results.csv".format(args.outpath, args.savename),index = False)

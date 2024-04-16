import argparse
import matplotlib.pyplot as plt
import statsmodels.stats.multitest as multitest
import numpy as np
import pandas as pd
import scipy.stats as stats
import tqdm

def get_args():
  parser = argparse.ArgumentParser(description="perform equivalence test")
  parser.add_argument('--savename', type=str, help = 'dataname to process')
  parser.add_argument('--outpath', type=str, help = 'path to save files', default="")
  parser.add_argument('--infile', type=str, help = 'path to file to analyze')
  parser.add_argument('--delta', type=float, help = 'delta value to use', default=0.05)
  args = parser.parse_args()
  return args


def perform_t_tests(group1_vals, group2_vals, delta):
    mean1 = np.mean(group1_vals)
    mean2 = np.mean(group2_vals)
    
    # perform standard t-test (testing for difference)
    # consider switching to Welch's t-test (current assumes equal variances)
    diff_pval = stats.ttest_ind(a=group1_vals, b = group2_vals).pvalue

    # get t statistics for equivalence test
    # (X1 - X2 - delta)/sqrt(s1^2/n1 + s2^2/n2)
    t1 = (np.mean(group1_vals) - np.mean(group2_vals) - delta)/np.sqrt(np.var(group1_vals)/len(group1_vals)+np.var(group2_vals)/len(group2_vals))

    # (X1 - X2 + delta)/sqrt(s1^2/n1 + s2^2/n2)
    t2 = (np.mean(group1_vals) - np.mean(group2_vals) + delta)/np.sqrt(np.var(group1_vals)/len(group1_vals)+np.var(group2_vals)/len(group2_vals))

    # find degrees of freedom (n1 + n2 - 1)
    df = len(group1_vals) + len(group2_vals)

    # find p value for t1 (less than)
    p1 = stats.t.cdf(t1, df)

    # find p value for t2 (greater than)
    # sf is survival function: 1 - cdf
    p2 = stats.t.sf(t2, df)

    equiv_pval = np.max([p1, p2])
#     print("diff pvalue: {}\nequiv pvalue: {}".format(diff_pval, equiv_pval))
    return diff_pval, equiv_pval

def plot_groups(group1_vals, group2_vals):
    mean1 = np.mean(group1_vals)
    mean2 = np.mean(group2_vals)
    plt.hist(group1_vals,alpha=0.4, color=u'#1f77b4')
    plt.axvline(x=mean1, color=u'#1f77b4',linestyle="--")
    plt.hist(group2_vals,alpha = 0.4, color = u'#ff7f0e')
    plt.axvline(x=mean2, color=u'#ff7f0e',linestyle="--")
    # plt.axvline(x=np.mean(group2_vals), color="black")
    plt.title("difference in means: {:.3f}".format(mean1 - mean2))
    plt.close()

def loop_over_genes(df, delta, plot):
  grouped_df = df.groupby("cell_type")
  
  out = {"gene" : [], "nnz_group1" : [], "nnz_group2" : [], "avg_group1" : [], 
         "avg_group2" : [], "diff_pval" : [], "equiv_pval" : []}
  
  for gene in tqdm.tqdm([x for x in df.columns if x != "cell_type"]):
  
      samples = grouped_df[gene].apply(list)
      diff_pval, equiv_pval = perform_t_tests(samples[0], samples[1], delta)
      if plot and ((diff_pval < 0.05) or (equiv_pval > 0.05)):
          print(gene)
          plt.hist(df[gene])
          plt.show()
          plot_groups(samples[0], samples[1])
      out["gene"].append(gene)
      out["nnz_group1"].append(np.count_nonzero(samples[0]))
      out["nnz_group2"].append(np.count_nonzero(samples[1]))
      out["avg_group1"].append(np.mean(samples[0]))
      out["avg_group2"].append(np.mean(samples[1]))
      out["diff_pval"].append(diff_pval)
      out["equiv_pval"].append(equiv_pval)
  
  return pd.DataFrame(out)

def process_out_df(out_df):
  out_df.dropna(inplace=True)
  for pval in ["diff", "equiv"]:
      out_df["{}_pval_adj".format(pval)] = multitest.multipletests(out_df["{}_pval".format(pval)], method="fdr_bh")[1]
      out_df["sig_{}".format(pval)] = False
      out_df.loc[out_df["{}_pval_adj".format(pval)] < 0.05, "sig_{}".format(pval)] = True
      
  out_df["sig_both"] = out_df["sig_diff"] & out_df["sig_equiv"]
  return out_df

def plot_results(out_df, outpath, dataname, delta):
  plt.hist(out_df["diff_pval_adj"])
  plt.title("{}\ndiff_pval_adj".format(dataname))
  plt.savefig("{}{}_{}_diff_pval.png".format(outpath, dataname, delta))
  plt.close()
  
  plt.hist(out_df["equiv_pval_adj"])
  plt.title("{}\nequiv_pval_adj".format(dataname))
  plt.savefig("{}{}_{}_equiv_pval.png".format(outpath, dataname, delta))
  plt.close()

  sub = out_df[(out_df["sig_diff"] == False) & (out_df["sig_equiv"] == False)]
  plt.plot(sub["diff_pval"], sub["equiv_pval"],marker="o",linestyle="", label = "neither sig")
  
  sub = out_df[(out_df["sig_diff"] == True) & (out_df["sig_equiv"] == False)]
  plt.plot(sub["diff_pval"], sub["equiv_pval"],marker="o",linestyle="", label = "sig diff")
  
  sub = out_df[(out_df["sig_diff"] == False) & (out_df["sig_equiv"] == True)]
  plt.plot(sub["diff_pval"], sub["equiv_pval"],marker="o",linestyle="", label = "sig equiv")
  
  sub = out_df[(out_df["sig_diff"] == True) & (out_df["sig_equiv"] == True)]
  plt.plot(sub["diff_pval"], sub["equiv_pval"],marker="o",linestyle="", label = "sig both")
  
  plt.xlabel("diff_pval")
  plt.ylabel("equiv_pval")
  plt.title("{}\ndelta: {}".format(dataname,delta))
  plt.legend()
  plt.savefig("{}{}_{}_comparison.png".format(outpath, dataname, delta))
  plt.close()



import argparse
import math
import matplotlib.pyplot as plt
import statsmodels.stats.multitest as multitest
from statsmodels.stats.proportion import proportion_confint
import numpy as np
import pandas as pd
import scipy.stats as stats
import tqdm

def get_args():
  parser = argparse.ArgumentParser(description="perform equivalence test")
  parser.add_argument('--savename', type=str, help = 'dataname to process')
  parser.add_argument('--outpath', type=str, help = 'path to save files', default="")
  parser.add_argument('--infile', type=str, help = 'path to file to analyze')
  parser.add_argument('--meta', type=str, help = 'metadata to analyze')
  parser.add_argument('--condition', type=str, help = 'column of metadata to use for differential analysis')
  parser.add_argument('--delta', type=float, help = 'delta value to use', default=0.05)
  parser.add_argument('--log_scale', action="store_true", help = 'include if you want your data to be log2 scaled')
  parser.add_argument('--normalize', action="store_true", help = 'normalize each sample gene-wise')

  args = parser.parse_args()
  return args

def normalize(df, scale = 100):
  df.iloc[:,:-1] = 100*df.iloc[:,:-1].div(df.iloc[:,:-1].sum(axis=1), axis=0)
  return df


def log2_scale(df):
  # add 1 to every entry (before log scaling)
  df.iloc[:,:-1] = np.log2(df.iloc[:,:-1] + 1)
  return df

def perform_t_tests(group1_vals, group2_vals, delta):
    #try:
    mean1 = np.mean(group1_vals)
    mean2 = np.mean(group2_vals)
    #except:
    #  print("ERROR\n{}\n{}".format(group1_vals, group2_vals))
    # perform standard t-test (testing for difference)
    # consider switching to Welch's t-test (current assumes equal variances)
    diff_pval = stats.ttest_ind(a=group1_vals, b = group2_vals).pvalue
    
    if math.isnan(diff_pval):
      print("is nan")
      print(group1_vals)
      print(group2_vals)
      print(mean1, mean2)

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

    if math.isnan(equiv_pval):
      print(t1)
      print(t2)
      print("denom",np.sqrt(np.var(group1_vals)/len(group1_vals)+np.var(group2_vals)/len(group2_vals)))
      print("is nan equiv")
      print(group1_vals)
      print(group2_vals)
      print(mean1, mean2)
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

def loop_over_genes(df, delta, plot = False):
  grouped_df = df.groupby("condition")
  
  out = {"gene" : [], "nnz_group1" : [], "nnz_group2" : [], "avg_group1" : [], 
         "avg_group2" : [], "diff_pval" : [], "equiv_pval" : []}
  
  for gene in [x for x in df.columns if x != "condition"]:
  
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

def process_out_df(out_df, delta):


  # need to deal with NAs somehow
  print(out_df)
  out_df.dropna(inplace=True)

  out_df["eff_size"] = (out_df["avg_group1"] - out_df["avg_group2"]).abs()
  for pval in ["diff", "equiv"]:
      try:
        out_df["{}_pval_adj".format(pval)] = multitest.multipletests(out_df["{}_pval".format(pval)], method="fdr_bh")[1]
      except:
        print(out_df)

      out_df["sig_{}".format(pval)] = False
      out_df.loc[out_df["{}_pval_adj".format(pval)] < 0.05, "sig_{}".format(pval)] = True

  # include effect size filter on significance
  out_df.loc[out_df["eff_size"] <= delta, "sig_diff"] = False    
  out_df.loc[out_df["eff_size"] > delta, "sig_equiv"] = False    
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

def summarize_df_significance(out_df, delta, verbose = True):
  num_genes = out_df.shape[0]
  sig_equiv = out_df["sig_equiv"].sum()
  sig_diff = out_df["sig_diff"].sum()
  sig_neither = (~(out_df["sig_diff"] | out_df["sig_equiv"])).sum()
  
  if verbose:
    print("\nnum_genes: {}\n# equivalent (eff <= {}): {} ({:.2f}%)\n# different (eff > {}): {} ({:.2f}%)\n# neither: {} ({:.2f}%)\n".format(num_genes, delta, sig_equiv, 100*sig_equiv/num_genes, delta, sig_diff, 100*sig_diff/num_genes, sig_neither, 100*sig_neither/num_genes))

  return num_genes, sig_equiv, sig_diff, sig_neither



def perform_tests_for_df(df, delta, verbose = True):
  out_df = loop_over_genes(df, delta)
  out_df = process_out_df(out_df, delta)
  num_genes, sig_equiv, sig_diff, sig_neither = summarize_df_significance(out_df, delta, verbose)

  return num_genes, sig_equiv, sig_diff, sig_neither
  
############################################################
############################################################
##########                                        ##########
##########            For Simulation              ##########
##########                                        ##########
############################################################
############################################################

def ttests(df, delta, norm, log_scale):
    
    # run t tests from script (should be what's done in perform tests)

    if norm:
      # normalize each row (sample) by the number of counts
        df = normalize(df)

    if log_scale:
        df = log2_scale(df)

    plot = False

    # testing out function you can use in Python code
    # perform_tests_for_df(df, args.delta)

    out_df = loop_over_genes(df, delta, plot)

    out_df = process_out_df(out_df, delta)

    return out_df["sig_diff"][0], out_df["sig_equiv"][0]

def simulate(numReads, num_trials, geneA_fracs, num_samples, delta, norm = True, log_scale = True, seed = 123):
    np.random.seed(seed)

    out = {"numReads" : [], "frac_sig_diff" : [], "frac_sig_diff_lower" : [], "frac_sig_diff_upper" : [], 
           "frac_sig_equiv" : [], "frac_sig_equiv_lower" : [], "frac_sig_equiv_upper" : [], 
           "frac_sig_inc" : [], "frac_sig_inc_lower" : [], "frac_sig_inc_upper" : []}

    for num in numReads:
        print("num reads: {}".format(num))
        # collect p vals for each
        num_sig_diff = 0
        num_sig_equiv = 0
        num_inconclusive = 0

        for i in range(num_trials):
            data = {"sample" : [], "A" : [], "B" : [], "cell_type" : []}
            count = 0
            for j in range(len(geneA_fracs)):
                for k in range(num_samples):
                    count += 1

                    # add 1 so it's never zero
                    n = np.random.poisson(num) + 1
                    A = np.random.binomial(n, geneA_fracs[j])
                    data["sample"].append("sample" + str(count))
                    data["A"].append(A)
                    data["B"].append(n - A)
                    data["cell_type"].append(j)
            df = pd.DataFrame(data).set_index("sample")

            # get results
            diff_sig, equiv_sig = ttests(df, norm, log_scale, delta)

            # save output
            num_sig_diff += diff_sig
            num_sig_equiv += equiv_sig
            num_inconclusive += 1 - diff_sig - equiv_sig

        # get frac significant
        out["numReads"].append(num)

    #    num_sig_diff = sum(1 for p in diff_pvals if p < alpha)
        lower_ci, upper_ci = proportion_confint(num_sig_diff, num_trials)
        out["frac_sig_diff"].append(num_sig_diff/num_trials)
        out["frac_sig_diff_lower"].append(lower_ci)
        out["frac_sig_diff_upper"].append(upper_ci)

    #    num_sig_equiv = sum(1 for p in equiv_pvals if p < alpha)
        lower_ci, upper_ci = proportion_confint(num_sig_equiv, num_trials)
        out["frac_sig_equiv"].append(num_sig_equiv/num_trials)
        out["frac_sig_equiv_lower"].append(lower_ci)
        out["frac_sig_equiv_upper"].append(upper_ci)


    #    num_sig_equiv = sum(1 for p in equiv_pvals if p < alpha)
        lower_ci, upper_ci = proportion_confint(num_inconclusive, num_trials)
        out["frac_sig_inc"].append(num_inconclusive/num_trials)
        out["frac_sig_inc_lower"].append(lower_ci)
        out["frac_sig_inc_upper"].append(upper_ci)
    #print(out)
    out = pd.DataFrame(out)
    display(out)
    
    max_cols = out[["frac_sig_diff", "frac_sig_equiv", "frac_sig_inc"]].idxmax(axis=1)
    winning_col = max_cols.iloc[-1]

    winning_numReads = out.loc[len(max_cols[~(max_cols == winning_col)]), "numReads"]
    
    return out, winning_col, winning_numReads

def plot_sim(out, geneA_fracs, delta, outpath):
    plt.errorbar(out["numReads"], out["frac_sig_diff"], yerr = [out["frac_sig_diff"] - out["frac_sig_diff_lower"],
            out["frac_sig_diff_upper"] - out["frac_sig_diff"] ], marker = "o", label = "diff")
    plt.errorbar(out["numReads"], out["frac_sig_equiv"], yerr = [out["frac_sig_equiv"] - out["frac_sig_equiv_lower"],
                out["frac_sig_equiv_upper"] - out["frac_sig_equiv"] ], marker = "o", label = "equiv")
    plt.errorbar(out["numReads"], out["frac_sig_inc"], yerr = [out["frac_sig_inc"] - out["frac_sig_inc_lower"],
                out["frac_sig_inc_upper"] - out["frac_sig_inc"] ], marker = "o", label = "inc")
    plt.xscale("log")
    plt.axhline(y=0.05, linestyle="--", color="k", label = "$y = 0.05$")
    plt.axhline(y=0.0, color="k")

    plt.legend()
    plt.xlabel("number of reads")
    plt.ylabel("fraction called in category")
    plt.title("gene A pop 1: {}\ngene A pop 2: {}\nDelta: {}".format(geneA_fracs[0], geneA_fracs[1], delta))
    plt.savefig("{}pop1_{}_pop2_{}_delt_{}.png".format(outpath, *geneA_fracs, delta),bbox_inches='tight')
    plt.show()

import argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as stats
import statsmodels.stats.multitest as multitest

def get_args():
  parser = argparse.ArgumentParser(description="perform equivalence test")
  parser.add_argument('--savename', type=str, help = 'dataname to process')
  parser.add_argument('--outpath', type=str, help = 'path to save files', default="")
  parser.add_argument('--infile', type=str, help = 'path to file to analyze')
  parser.add_argument('--meta', type=str, help = 'metadata to analyze')
  parser.add_argument('--condition', type=str, help = 'column of metadata to use for differential analysis')
  parser.add_argument('--delta', type=float, help = 'delta value to use', default=1)
  parser.add_argument('--deseq2_results', type=str, help = 'include path to deseq2 results file if you want plots', default="none")

  args = parser.parse_args()
  return args

def load_data(infile, meta, condition):

    df = pd.read_csv(infile, index_col = 0)
    meta = pd.read_csv(meta, index_col = 0)

    # only include samples that are in the metadata
    df = df[meta.index]

    # find which samples are in each group
    groups = []
    for name, group in meta.groupby(condition):
        groups.append(list(group.index))

    # Remove genes (rows) that have <= 1 nonzero entry
    nnz_plus = ((df != 0).sum(axis=1) > 1)
    df = df.loc[nnz_plus[nnz_plus].index,]

    # log 2 scale all entries (after adding pseudocount of 1)
    df = np.log2(df + 1)
    return df, groups

def calc_basic_stats(df, groups):
    calc_df = pd.concat([df[groups[0]].mean(axis=1), # mean group 1
            df[groups[1]].mean(axis=1), # mean group 2
            df[groups[0]].var(axis=1), # sample variance group 1
            df[groups[1]].var(axis=1), # sample variance group 2
            (df[groups[0]] != 0).sum(axis=1), # nnz group 1
            (df[groups[1]] != 0).sum(axis=1) # nnz group 1
            ],axis=1).rename({0 : "muf", 1 : "mug", 2 : "sf2", 3 : "sg2", 4 : "nnz_group1", 5 : "nnz_group2"},axis=1)
    calc_df["nf"] = len(groups[0])
    calc_df["ng"] = len(groups[1])
    return calc_df


def degrees_freedom_welch(df):
    # see methods for formula
    # deviation from simple nf + ng - 2 because we're not assuming equal variances
    df["nu"] = (df["sf2"]/df["nf"] + df["sg2"]/df["ng"])**2/(df["sf2"]**2/(df["nf"]**2*(df["nf"] - 1)) + df["sg2"]**2/(df["nf"]**2*(df["ng"] - 1)))
    #df["nu"] = df["nf"]  + df["ng"]  - 2

def calc_t(df, delta):
    df["t{}".format(delta)] = ((df["muf"] - df["mug"]) + delta)/(np.sqrt((df["sf2"]/df["nf"]) + (df["sg2"]/df["ng"])))

def calc_diff_pval(df):
    df["cdf_t0"] = stats.t.cdf(df["t0"],df["nu"])
    df["diff_pval"] = 2*pd.concat([df["cdf_t0"], 1 - df["cdf_t0"]],axis=1).min(axis=1)

def calc_equiv_pval(df, delta):
    df["cdf_t1"] = stats.t.cdf(df["t{}".format(-delta)], df["nu"])
    df["cdf_t2"] = 1 - stats.t.cdf(df["t{}".format(delta)], df["nu"])
    df["equiv_pval"] = pd.concat([df["cdf_t1"], df["cdf_t2"]],axis=1).max(axis=1)

def adj_pvals(df,col):
    df[col + "_adj"] = None
    
    # don't calculate adjusted p values if all of the p values are NA
    if df[col].dropna().shape[0] > 0:
        df.loc[~df[col].isna(),col + "_adj"] = multitest.multipletests(df[col].dropna(), method="fdr_bh")[1]

def format_outdf(calc_df, delta):
    savecols = ["avg_group1_log2", "avg_group2_log2", "nnz_group1", "nnz_group2", 
            "eff_size", "avg_group1", "avg_group2", "fold_change", 
            "diff_pval", "equiv_pval", "diff_pval_adj", "equiv_pval_adj",
            "sig_diff", "sig_equiv", "category"]
    calc_df.rename({"muf" : "avg_group1_log2", "mug" : "avg_group2_log2"},axis=1,inplace=True)
    calc_df["eff_size"] =  calc_df["avg_group2_log2"] - calc_df["avg_group1_log2"]
    for i in range(1,3):
        calc_df["avg_group{}".format(i)] = 2**calc_df["avg_group{}_log2".format(i)]
    calc_df["fold_change"] = 2**calc_df["eff_size"]
    for v in ["diff", "equiv"]:
        calc_df["sig_" + v] = False
        calc_df.loc[calc_df["{}_pval_adj".format(v)] < 0.05, "sig_{}".format(v)] = True
    # include effect size filter on significance    
    calc_df.loc[(-delta <= calc_df["eff_size"]) & (calc_df["eff_size"] <= delta), "sig_diff"] = False    
    calc_df.loc[(calc_df["eff_size"] > delta) | (calc_df["eff_size"] < -delta), "sig_equiv"] = False 
    calc_df["category"] = "inconclusive"
    calc_df.loc[calc_df["sig_diff"],"category"] = "different"
    calc_df.loc[calc_df["sig_equiv"],"category"] = "equivalent"
    return calc_df[savecols]

def plot_deseq2(calc_df, deseq2_results, savepath):
    if deseq2_results == "none":
        return
    deseq2 = pd.read_csv(deseq2_results,index_col = 0)
    merged = pd.concat([calc_df, deseq2], join = "inner",axis=1)
    
    # sometimes interpreted as object; causes problems with correlation
    merged["diff_pval_adj"] = pd.to_numeric(merged["diff_pval_adj"], errors='coerce')

    plt.plot(merged["padj"],merged["diff_pval_adj"], marker = "o", linestyle="", alpha = 0.1)
    plt.xlabel("DESeq2 adjusted p value")
    plt.ylabel("adjusted differential expression p value")
    plt.title("Pearson correlation {:.3f}".format(merged["padj"].corr(merged["diff_pval_adj"])))
    plt.savefig("{}_deseq2_pvals.png".format(savepath), bbox_inches="tight", dpi=300)
    plt.close()

    col1 = "log2FoldChange"
    col2 = "eff_size"
    plt.plot(merged[col1],merged[col2], marker = "o", linestyle="", alpha = 0.1)
    plt.title("Pearson correlation {:.3f}".format(merged[col1].corr(merged[col2])))

    plt.xlabel("DESeq2 log2 fold change")
    plt.ylabel("current log2 fold change")
    plt.savefig("{}_deseq2_foldchange.png".format(savepath), bbox_inches="tight", dpi=300)
    plt.close()

def round_interval(interval):
    return pd.Interval(round(interval.left, 1), round(interval.right, 1), closed=interval.closed)

def plot_results(calc_df, savepath, delta):
    k = 10
    calc_df["avg"] = calc_df[["avg_group1", "avg_group2"]].mean(axis=1)
    calc_df["depth_quantile"] = pd.qcut(calc_df['avg'], k, labels=False)
    calc_df["depth_quantile_name"] = pd.qcut(calc_df['avg'], k)

    # Apply the rounding function to each interval in the column
    calc_df["depth_quantile_name_rounded"] = calc_df["depth_quantile_name"].apply(round_interval)

    # find number of genes in each category at each sequencing depth
    out_dict = {"quant_name" : [], "quant" : [], "num_genes" : [], "num_diff" : [], "num_equiv" : [], "num_incon" : []}
    for quant, quantdf in calc_df.groupby("depth_quantile"):
        out_dict["quant_name"].append(quantdf["depth_quantile_name_rounded"].astype(str)[0])
        out_dict["quant"].append(quant)
        out_dict["num_genes"].append(quantdf.shape[0])
        out_dict["num_diff"].append(quantdf[quantdf["category"] == "different"].shape[0])
        out_dict["num_equiv"].append(quantdf[quantdf["category"] == "equivalent"].shape[0])
        out_dict["num_incon"].append(quantdf[quantdf["category"] == "inconclusive"].shape[0])
    out = pd.DataFrame(out_dict)
    out["frac_diff"] = out["num_diff"]/out["num_genes"]
    out["frac_equiv"] = out["num_equiv"]/out["num_genes"]
    out["frac_incon"] = out["num_incon"]/out["num_genes"]

    # plot each category (with Wald binomial error bars, 95% confidence interval)
    for cat in ["diff", "equiv", "incon"]:
        plt.errorbar(out["quant"], out["frac_" + cat], yerr=1.96*np.sqrt(out["frac_" + cat]*(1 - out["frac_" + cat])/out["num_genes"]), fmt="o", label = cat)
    plt.ylabel("fraction of genes in category")
    plt.xlabel("log2 read depth range")
    plt.xticks(out["quant"], out["quant_name"], rotation='vertical')
    plt.legend()
    plt.title("Fraction from each category by read depth")
    plt.savefig("{}_readdepth.png".format(savepath), bbox_inches="tight", dpi=300)
    plt.close()

    calc_df["diff_pval_adj"] = pd.to_numeric(calc_df["diff_pval_adj"], errors='coerce')
    plt.axhline(y=1.301,color="lightgray", label = "0.05", linestyle="--")
    plt.axvline(x = delta, color = "lightgray", linestyle="-")
    plt.axvline(x = -delta, color = "lightgray", linestyle="-", label = "+/- delta")
    for cat, group in calc_df.groupby("category"):
        plt.plot(group["eff_size"], -np.log10(group["diff_pval_adj"]), marker = "o", linestyle="",alpha=0.1, label = cat)

    plt.legend()

    plt.xlabel("log2 fold change")
    plt.ylabel("-log10 difference p value")
    plt.title("volcano plot")
    plt.savefig("{}_volcano.png".format(savepath), bbox_inches="tight", dpi=300)
    plt.close()

def perform_full_analysis(infile, meta, condition, delta):

    # load data
    df, groups = load_data(infile, meta, condition)

    # calculate mean, variance, nnz, etc
    calc_df = calc_basic_stats(df, groups)

    # Calculate degrees of freedom
    degrees_freedom_welch(calc_df)

    # calculate 3 necessary t values
    calc_t(calc_df, 0)
    calc_t(calc_df, delta)
    calc_t(calc_df, -delta)

    # calculate p values
    calc_diff_pval(calc_df)
    calc_equiv_pval(calc_df, delta)

    # adjust p values
    for v in ["diff", "equiv"]:
        adj_pvals(calc_df,v + "_pval")

    calc_df.sort_values("diff_pval", inplace=True)

    calc_df = format_outdf(calc_df, delta)
    return calc_df
from equiv_test_funcs_vectorized import *

def main():
    args = get_args()
    savepath = "{}{}_{}_{}".format(args.outpath, args.savename, args.condition, args.delta)

    # moved to separate script to allow simulation to run
    calc_df = perform_full_analysis(
        args.infile,
        args.meta,
        args.condition,
        args.delta,
        low_count_threshold=args.low_count_threshold,
        min_nonzero_per_group=args.min_nonzero_per_group,
    )

    print(pd.crosstab(calc_df['sig_equiv'], calc_df['sig_diff']))
    plot_deseq2(calc_df, args.deseq2_results, savepath)
    plot_results(calc_df, savepath, args.delta)
    calc_df.to_csv("{}_results.csv".format(savepath))
    print("Number of genes: {}".format(calc_df.shape[0]))
    print("Number differentially expressed genes: {}".format((calc_df['category'] == 'different').sum()))
    print("Number equivalently expressed genes: {}".format((calc_df['category'] == 'equivalent').sum()))
    print("Number of inconclusive genes: {}".format((calc_df['category'] == 'inconclusive').sum()))
    print("Number of low-information genes: {}".format((calc_df['category'] == 'low_information').sum()))
main()

from equiv_test_funcs_vectorized import *

def main():
    args = get_args()

    # moved to separate script to allow simulation to run
    calc_df = perform_full_analysis(args.infile, args.meta, args.condition, args.delta)

    calc_df.to_csv("{}{}_{}_{}_results.csv".format(args.outpath, args.savename, args.condition, args.delta))
main()
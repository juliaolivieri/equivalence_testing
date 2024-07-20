import pandas as pd

identifier = "GSE117891"

data_files = pd.read_csv("../data/data_files.csv")

row = data_files[data_files["id"] == identifier]

delta = 1

savename = "{}_{}_{}".format(identifier, row["meta_column"].iloc[0], delta)

outpath_equiv = "~/Desktop/Research/equivalence_testing_output/scripts/output/perform_test_tabular_data/"

outpath_deseq = "~/Desktop/Research/equivalence_testing_output/scripts/output/deseq_normalization/{}".format(identifier)

f = open("run_auto.sh","w")

f.write("#!/usr/bin/env bash\n\ndate\n\n")

f.write("Rscript deseq_normalization.R {} {} {} {}".format(row["data_path"].iloc[0], row["twocat_meta_path"].iloc[0], row["meta_column"].iloc[0],outpath_deseq))

f.write("\n\ndate\necho 'normalized'\n\n")

f.write("python perform_test_tabular_data.py --savename {} --outpath {} --infile {} --delta {} --log_scale --meta {} --condition {}".format(
    savename, outpath_equiv, row["data_path"].iloc[0], delta, row["twocat_meta_path"].iloc[0],row["meta_column"].iloc[0] 
))


f.write("\n\ndate\necho 'complete'\n\n")
f.close()
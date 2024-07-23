import pandas as pd

identifier = "GSE206932"

data_files = pd.read_csv("../data/data_files.csv")

row = data_files[data_files["id"] == identifier]

delta = 1

epsilon = 0.1

savename = "{}_{}_{}_{}".format(identifier, row["meta_column"].iloc[0], delta, epsilon)

outpath_equiv = "~/Desktop/Research/equivalence_testing_output/scripts/output/perform_test_tabular_data/"

outpath_deseq = "~/Desktop/Research/equivalence_testing_output/scripts/output/deseq_normalization/{}".format(identifier)

f = open("run_auto.sh","w")

f.write("#!/usr/bin/env bash\n\ndate\n\n")

f.write("Rscript deseq_normalization.R {} {} {} {}".format(row["data_path"].iloc[0], row["twocat_meta_path"].iloc[0], row["meta_column"].iloc[0],outpath_deseq))

f.write("\n\ndate\necho 'normalized'\n\n")

f.write("python perform_test_tabular_data.py --savename {} --outpath {} --infile {}_deseq2_normalized_counts.csv --delta {} --epsilon {} --log_scale --meta {} --condition {}".format(
    savename, outpath_equiv, outpath_deseq, delta, epsilon, row["twocat_meta_path"].iloc[0],row["meta_column"].iloc[0] 
))


f.write("\n\ndate\necho 'complete'\n\n")
f.close()
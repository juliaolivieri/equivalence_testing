import pandas as pd
import numpy as np

identifier = "onco"

DESEQ2 = True

data_files = pd.read_csv("../data/data_files.csv")

row = data_files[data_files["id"] == identifier]

delta = 1
#delta = round(np.log2(3/2), 3)
print("delta", delta)

outpath_equiv = "output/equiv_test_vectorized/"

outpath_deseq = "output/deseq_normalization/{}".format(identifier)

f = open("run_auto.sh","w")

f.write("#!/usr/bin/env bash\n\ndate\n\n")

if DESEQ2:
    f.write("Rscript deseq_normalization.R {} {} {} {}".format(row["data_path"].iloc[0], row["twocat_meta_path"].iloc[0], row["meta_column"].iloc[0],outpath_deseq))

    f.write("\n\ndate\necho 'normalized'\n\n")

f.write("python equiv_test_vectorized.py --savename {} --outpath {} --infile {}_deseq2_normalized_counts.csv  --meta {} --condition {} --delta {} --deseq2_results {}_{}_deseq2_results.csv".format(
    identifier, outpath_equiv, outpath_deseq, row["twocat_meta_path"].iloc[0],row["meta_column"].iloc[0], delta,
    outpath_deseq, row["meta_column"].iloc[0] 
        ))


f.write("\n\ndate\necho 'complete'\n\n")
f.close()

import anndata
import pandas as pd

h5ad_file = "/exports/home/jolivieri/data/single_cell/HLCA/7a3f08f9-5d07-4ddd-a8fe-5967dd34f35f.h5ad"
outpath = "/exports/home/jolivieri/equivalence_testing_output/scripts/output/write_HLCA_csv/"

data = anndata.read_h5ad(h5ad_file)

df = pd.DataFrame.sparse.from_spmatrix(data.raw.X, columns = data.var["feature_name"], index = data.obs.index)

df.to_csv("{}HLCA_raw.csv".format(outpath))

# save metadata for genes
data.raw.var.to_csv("{}HLCA_genes.csv".format(outpath))

# save metadata for cells
data.obs.to_csv("{}HLCA_meta.csv".format(outpath))


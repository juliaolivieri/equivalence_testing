import anndata
import argparse
import numpy as np
import pandas as pd

def get_args():
  parser = argparse.ArgumentParser(description="subset HLCA data")
  parser.add_argument('--cell_types', type=str, nargs = '*', help = 'cell types to include')
  parser.add_argument('--gene_cutoff', type=int, help = 'sum of gene column must be greater than this value to be included')
  parser.add_argument('--raw', action="store_true", help = 'dictates whether to use raw counts')
  args = parser.parse_args()
  return args

def main():
  args = get_args()
  outpath = "/exports/home/jolivieri/equivalence_testing_output/scripts/output/HLCA_data_processing/"
  dataname = "HLCA_normal"
#  cell_types = ["CD4-positive, alpha-beta T cell","CD8-positive, alpha-beta T cell"]
#  gene_cutoff = 100
  print("args.cell_types",args.cell_types)
  print("args.gene_cutoff",args.gene_cutoff)

  type_suff = "_{}".format(args.gene_cutoff)
  for x in args.cell_types:
      type_suff += "_" + x.replace(" ", "-").replace(",","")

  # read in normal lung data
  h5ad_file = "/exports/home/jolivieri/data/single_cell/HLCA/7a3f08f9-5d07-4ddd-a8fe-5967dd34f35f.h5ad"
  data = anndata.read_h5ad(h5ad_file)
  




  # subset to only tissue with the largest number of cells
  print(data.obs["tissue"].value_counts())
#  data = data[data.obs.tissue == "lung parenchyma"]


  # subset to only two chosen cell types
  print(data.obs.cell_type.value_counts().head())
#  data = data[data.obs.cell_type.isin(args.cell_types)]

  cell_condition = (data.obs.cell_type.isin(args.cell_types)) & (data.obs.tissue == "lung parenchyma")  


  # subset to only genes with at least gene_cutoff count across these cells
  # doesn't need to be raw b/c just looking for zeros
  data.var["gene_count"] = data.X.sum(axis=0).tolist()[0]
#  data = data[:,data.var.gene_count > args.gene_cutoff]

  gene_condition = data.var.gene_count > 1000
  cell_indices = np.where(cell_condition)[0]
  gene_indices = np.where(gene_condition)[0]

  # shape of current df
  print("data.X.shape",data.X.shape)

  # save metadata for genes
  data.var.to_csv("{}{}{}_genes.csv".format(outpath, dataname, type_suff))

  # save metadata for cells
  data.obs.to_csv("{}{}{}_meta.csv".format(outpath, dataname, type_suff))

  # create pandas dataframe with just the selected data and save it

  if args.raw:
    df = pd.DataFrame.sparse.from_spmatrix(data.raw.X[cell_indices,:][:,gene_indices], columns = data.var["feature_name"][gene_indices], index = data.obs.iloc[cell_indices].index)    

    df.to_csv("{}{}{}_raw_data.csv".format(outpath, dataname, type_suff))

  else:
    df = pd.DataFrame.sparse.from_spmatrix(data.X[cell_indices,:][:,gene_indices], columns = data.var["feature_name"][gene_indices], index = data.obs.iloc[cell_indices].index)    

    df.to_csv("{}{}{}_data.csv".format(outpath, dataname, type_suff))

main()

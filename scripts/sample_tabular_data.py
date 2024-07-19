import argparse
import pandas as pd
import random

def get_args():
  parser = argparse.ArgumentParser(description="subset tabular data")
  parser.add_argument('--dataname', type=str, help = 'dataname to process')
  parser.add_argument('--num_cells', type=int, help = 'number of cells to include', default = 1000)
  parser.add_argument('--num_genes', type=int, help = 'number of genes to include', default = 100)
  parser.add_argument('--all_genes', action="store_true", help = "run for all genes")
  parser.add_argument('--all_cells', action="store_true", help = "run for all cells")
  parser.add_argument('--raw',action="store_true", help = "indicates that we're using the raw data")


  args = parser.parse_args()
  return args


def main():
  args = get_args()
  outpath = "../../../equivalence_testing_output/scripts/output/sample_tabular_data/"
  inpath = "../../../equivalence_testing_output/scripts/output/HLCA_data_processing/"

  cell_suff = args.num_cells
  gene_suff = args.num_genes
  if args.all_cells:
    cell_suff = "all"
  if args.all_genes:
    gene_suff = "all"

  if args.raw:
    df = pd.read_csv("{}{}_raw_data.csv".format(inpath, args.dataname), index_col=0)
  else:
    df = pd.read_csv("{}{}_data.csv".format(inpath, args.dataname), index_col=0)
  meta = pd.read_csv("{}{}_meta.csv".format(inpath, args.dataname), index_col = 0)
  df["cell_type"] = meta["cell_type"]
  df.to_csv("{}{}.csv".format(outpath, args.dataname))


  random.seed(123)
  
  # sample the given number of cells and genes
  # include cell type column
  if args.all_cells:
    sub_df = df
  else:
    sub_df = df.sample(axis=0, n=args.num_cells)
  cell_type = sub_df["cell_type"]
  if not args.all_genes:
    sub_df = sub_df.drop("cell_type", axis=1).sample(axis=1,n=args.num_genes)
  sub_df["cell_type"] = cell_type
  
  if args.raw:

    sub_df.to_csv("{}{}_{}_{}_raw.csv".format(outpath, args.dataname, cell_suff, gene_suff))
  else:

    sub_df.to_csv("{}{}_{}_{}.csv".format(outpath, args.dataname, cell_suff, gene_suff))

main()


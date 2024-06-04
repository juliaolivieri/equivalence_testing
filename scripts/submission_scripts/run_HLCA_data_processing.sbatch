#!/bin/bash
#SBATCH --time=24:30:00                     # how much time to run
#SBATCH --job-name=HLCA_data_processing                 # name the job jupyter host
#SBATCH --output=/exports/home/jolivieri/equivalence_testing_output/scripts/job_output/HLCA_processing.%j.out
#SBATCH --error=/exports/home/jolivieri/equivalence_testing_output/scripts/job_output/HLCA_processing.%j.err

source ~/.bashrc
conda deactivate
conda activate jup_env

#CT1="'classical monocyte'"
#CT2="'non-classical monocyte'"

CT1="'alveolar macrophage'"
CT2="'elicited macrophage'"

#CT1="'type I pneumocyte'"
#CT2="'type II pneumocyte'"

#CT1="'capillary endothelial cell'"
#CT2="'mast cell'"

#CT1="'CD4-positive, alpha-beta T cell'"
#CT2="'CD8-positive, alpha-beta T cell'"

#CT1="'ciliated columnar cell of tracheobronchial tree'"
#CT2="'CD8-positive, alpha-beta T cell'"

GENE_CUTOFF=1000

a="python -u ../HLCA_data_processing.py --cell_types ${CT1} ${CT2} --gene_cutoff ${GENE_CUTOFF} --raw"

date
echo $a
eval $a
date

NUM_CELLS=1000
NUM_GENES=100
DATANAME="HLCA_normal_1000_alveolar-macrophage_elicited-macrophage"
#DATANAME="HLCA_normal_${GENE_CUTOFF}_ciliated-columnar-cell-of-tracheobronchial-tree_CD8-positive-alpha-beta-T-cell"
#DATANAME="HLCA_normal_1000_classical-monocyte_non-classical-monocyte"
#DATANAME="HLCA_normal_1000_ciliated-columnar-cell-of-tracheobronchial-tree_CD8-positive-alpha-beta-T-cell"
#DATANAME="HLCA_normal_1000_capillary-endothelial-cell_mast-cell"

a="python -u ../sample_tabular_data.py --num_cells ${NUM_CELLS} --num_genes ${NUM_GENES} --dataname ${DATANAME} --raw"

date
echo $a
eval $a
date


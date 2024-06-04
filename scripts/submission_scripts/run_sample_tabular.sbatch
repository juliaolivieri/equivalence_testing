#!/bin/bash
#SBATCH --time=2:30:00                     # how much time to run
#SBATCH --job-name=sample_tabular_data                 # name the job jupyter host
#SBATCH --output=job_output/sample_tabular_data.%j.out
#SBATCH --error=job_output/sample_tabular_data.%j.err

source ~/.bashrc
conda deactivate
conda activate jup_env

NUM_CELLS=1000
NUM_GENES=1000
#DATANAME="HLCA_normal_1000_alveolar-macrophage_elicited-macrophage"
#DATANAME="HLCA_normal_1000_capillary-endothelial-cell_mast-cell"
DATANAME="HLCA_normal_1000_classical-monocyte_non-classical-monocyte"

a="python -u ../sample_tabular_data.py --num_cells ${NUM_CELLS} --num_genes ${NUM_GENES} --dataname ${DATANAME} --raw"

date
echo $a
eval $a
date


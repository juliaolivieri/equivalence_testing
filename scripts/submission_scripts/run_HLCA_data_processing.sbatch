#!/bin/bash
#SBATCH --time=8:30:00                     # how much time to run
#SBATCH --job-name=HLCA_data_processing                 # name the job jupyter host
#SBATCH --output=job_output/HLCA_processing.%j.out
#SBATCH --error=job_output/HLCA_processing.%j.err

source ~/.bashrc
conda deactivate
conda activate jup_env

CT1="'classical monocyte'"
CT2="'non-classical monocyte'"
GENE_CUTOFF=100

a="python -u ../HLCA_data_processing.py --cell_types ${CT1} ${CT2} --gene_cutoff ${GENE_CUTOFF}"

date
echo $a
eval $a
date


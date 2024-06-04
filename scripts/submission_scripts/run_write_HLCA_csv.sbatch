#!/bin/bash
#SBATCH --time=8:30:00                     # how much time to run
#SBATCH --job-name=write_csv                 # name the job jupyter host
#SBATCH --output=/exports/home/jolivieri/equivalence_testing_output/scripts/job_output/write_csv.%j.out
#SBATCH --error=/exports/home/jolivieri/equivalence_testing_output/scripts/job_output/write_csv.%j.err

source ~/.bashrc
conda deactivate
conda activate jup_env

a="python ../write_HLCA_csv.py"

date
echo $a
eval $a
date

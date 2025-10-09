#!/bin/bash

# Note: below script is used to gather computational model data points into a single file, and is designed to run on a Slurm-based cluster
# This is an example only. modifications may be required depending on your cluster specifications and conda installation.

#SBATCH -c 1                     # request 1 core
#SBATCH -t 0-01:00                  # 1 hour run time
#SBATCH -p short                  # use short partition
#SBATCH --mem 16G               # use 32Gb
#SBATCH --array=1-1             # Run array for indexes 1-n
#SBATCH -o slurm_output/slurm-%A_%a.out

module load conda/miniforge3/24.11.3-0
conda activate /path/to/.conda/envs/your_env


# Example commands to gather files for several parameterizations
# For example purposes only. Requires modification according to run names, file names, etc.

python gather_grid.py  --dir=grid_4param/ --name=grid_PL_final --write=pickle --final --B_dist

python gather_grid.py  --dir=grid_3M_interp/ --name=grid_3M_PL_interp_final --write=pickle --final --B_dist

python gather_grid.py  --dir=grid_4param_xaxis/ --name=grid_4param_xaxis --final --write=pickle --sep='tE|tC|LE|LC' --check_complete='grid_group_jobs_xaxis.pickle'

python gather_grid.py  --dir=grid_sparse_constantspeedup/ --name=grid_PL_constant_final --write=pickle --final --col_list='[8,10,12,14]' --check_complete='jobs_constantspeedup.pickle'
python gather_grid.py  --dir=grid_sparse_constantspeedup/ --name=grid_PL_constant --write=pickle --col_list='[8,10,12,14]'

# Stochastic runs
python gather_grid.py  --dir=grid_top95c_pro_stoch/1/ --name=grid_top95c_stoch_1 --write=pickle
python gather_grid.py  --dir=grid_top95c_pro_stoch/2/ --name=grid_top95c_stoch_2 --write=pickle
python gather_grid.py  --dir=grid_top95c_pro_stoch/3/ --name=grid_top95c_stoch_3 --write=pickle

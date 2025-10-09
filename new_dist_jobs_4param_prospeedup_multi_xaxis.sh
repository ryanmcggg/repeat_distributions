#!/bin/bash
#SBATCH -c 20                     # request 20 cores (max for short)
#SBATCH -t 0-06:00                  # 6 hours run time (12hr max for short)
#SBATCH -p short                  # use short partition
#SBATCH --mem 1000M               # use 2000Mb
#SBATCH --array=0-19             # Run array for indexes 1-999 (Slurm doesn't accept more than 10000 jobs per batch)
#SBATCH -o slurm_output/slurm-%A_%a.out

module load conda/miniforge3/24.11.3-0 || module load conda/miniforge3 || sleep 5 && module --ignore-cache load conda/miniforge3
conda activate /home/rjm44/.conda/envs/rep_env_3

# Run with '${mult}' etc. as the variable names
python3 computational_model_script_4params.py --jobgroup=$SLURM_ARRAY_TASK_ID  --jobfile='grid_group_jobs_xaxis_v2.pickle' --starting_counts='subonly' --rates_function='powerlaw_lambda' --speedup=3 --A_bins=200 --B_bins=200 --dir=grid_4param_xaxis_v3/

#python3 computational_model_script_4params.py --jobgroup=-1  --jobfile='jobs_incomplete.pickle' --starting_counts='subonly' --rates_function='powerlaw_x' --speedup=3 --A_bins=200 --B_bins=200 --dir=grid_4param_xaxis_v3/

python3 computational_model_script_4params.py --jobgroup=$SLURM_ARRAY_TASK_ID  --jobfile='jobs_incomplete_split.pickle' --starting_counts='subonly' --rates_function='powerlaw_lambda' --speedup=3 --A_bins=200 --B_bins=200 --dir=grid_4param_xaxis_v3/

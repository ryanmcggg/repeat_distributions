#!/bin/bash
#SBATCH -c 20                     # request 20 cores (max for short)
#SBATCH -t 0-12:00                  # 5 days run time (12hr max for short, 5 days max for medium)
#SBATCH -p short                  # use short or medium partition
#SBATCH --mem 1000M               # use 2000Mb
#SBATCH --array=26-47             # Run array for indexes 1-999 (Slurm doesn't accept more than 10000 jobs per batch)
#SBATCH -o slurm_output/slurm-%A_%a.out

module load conda/miniforge3/24.11.3-0 || module load conda/miniforge3 || sleep 5 && module --ignore-cache load conda/miniforge3
conda activate /home/rjm44/.conda/envs/rep_env_3

# Run with '${mult}' etc. as the variable names
#python computational_model_script_4params.py --jobfile='grid_group_jobs_sparse.pickle' --jobgroup=$SLURM_ARRAY_TASK_ID  --starting_counts='uniform' --speedup=1 --ceiling=0.1 --rounds=9 --A_bins=200 --B_bins=200  --dir=grid_sparse_constantspeedup_uniform/ --constantspeedup

#python computational_model_script_4params.py --jobfile='grid_group_jobs_top95c.pickle' --jobgroup=$SLURM_ARRAY_TASK_ID  --starting_counts='subonly' --speedup=1 --ceiling=0.1 --rounds=9 --A_bins=200 --B_bins=200  --dir=grid_constantspeedup_top95c/ --constantspeedup

# use --array=11-14
#python computational_model_script_4params.py --jobfile='grid_group_jobs_top95c.pickle' --jobgroup=$SLURM_ARRAY_TASK_ID  --starting_counts='uniform' --speedup=1 --ceiling=0.1 --rounds=9 --A_bins=200 --B_bins=200  --dir=grid_constantspeedup_top95c/ --constantspeedup

# use array=0-47
#python computational_model_script_4params.py --jobfile='3Mgrid_jobs_constantspeedup_short.pickle' --jobgroup=$SLURM_ARRAY_TASK_ID  --starting_counts='subonly' --speedup=1 --ceiling=0.1 --rounds=9 --A_bins=200 --B_bins=200  --dir=grid_sparse_constantspeedup/ --constantspeedup
python computational_model_script_4params.py --jobfile='3Mgrid_jobs_constantspeedup_short.pickle' --jobgroup=$SLURM_ARRAY_TASK_ID  --starting_counts='uniform' --speedup=1 --ceiling=0.1  --rounds=9 --A_bins=200 --B_bins=200  --dir=grid_sparse_constantspeedup_uniform/ --constantspeedup

#!/bin/bash
#SBATCH -c 20                     # request 20 cores (max for short)
#SBATCH -t 0-09:00                  # 6 hours run time (12hr max for short)
#SBATCH -p short                  # use short partition
#SBATCH --mem 800M               # use 2000Mb
#SBATCH --array=0-134             # Run array for indexes 1-999 (Slurm doesn't accept more than 10000 jobs per batch)
#SBATCH -o slurm_output/slurm-%A_%a.out

module load conda/miniforge3/24.11.3-0 || module load conda/miniforge3 || sleep 5 && module --ignore-cache load conda/miniforge3
conda activate /home/rjm44/.conda/envs/rep_env_3

# Run with '${mult}' etc. as the variable names
#python computational_model_script_4params.py --jobgroup=$SLURM_ARRAY_TASK_ID  --starting_counts='uniform' --boundary_count=200000 --speedup=5 --rounds=1 --A_bins=200 --B_bins=200 --interp --dir=grid_4param_interp_uniform/

#python computational_model_script_4params.py --jobgroup=$SLURM_ARRAY_TASK_ID  --starting_counts='subonly' --jobfile='fullgrid_jobs_prospeedup.pickle' --speedup=3 --A_bins=200 --B_bins=200 --dir=grid_4param_v6/

# use array=0-200
#python computational_model_script_4params.py --jobgroup=$SLURM_ARRAY_TASK_ID  --starting_counts='subonly' --jobfile='3Mgrid_jobs_prospeedup.pickle' --speedup=3 --A_bins=200 --B_bins=200 --interp  --dir=grid_3M_interp/

# use array = 0-999
#python3 computational_model_script_4params.py --jobgroup=$SLURM_ARRAY_TASK_ID  --starting_counts='subonly' --rates_function='log' --jobfile='grid_group_jobs_log.pickle' --speedup=3 --A_bins=200 --B_bins=200 --dir=grid_4param_log_v4/
python3 computational_model_script_4params.py --jobgroup=$SLURM_ARRAY_TASK_ID  --starting_counts='subonly' --rates_function='log' --jobfile='jobs_incomplete_split4.pickle' --speedup=3 --A_bins=200 --B_bins=200 --dir=log_annex/

# use array=0-14
#python computational_model_script_4params.py --jobgroup=$SLURM_ARRAY_TASK_ID  --starting_counts='subonly' --jobfile='grid_group_jobs_top95c2.pickle' --speedup=3 --A_bins=200 --B_bins=200 --dir=grid_top95c_pro_stoch/10/ --stochastics=10


# use array=0-200
#python computational_model_script_4params.py --jobgroup=$SLURM_ARRAY_TASK_ID  --starting_counts='subonly' --jobfile='3Mgrid_jobs_prospeedup.pickle' --speedup=3 --A_bins=200 --B_bins=200 --stochastics=1  --dir=grid_3M_stoch/

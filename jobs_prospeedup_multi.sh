#!/bin/bash

# Note: below script is used to run computational model data points in parallel, and is designed to run on a Slurm-based cluster
# This is an example only. Modifications may be required depending on your cluster specifications and conda installation.

#SBATCH -c 20                    # request 20 cores
#SBATCH -t 0-12:00               # 12 hours run time
#SBATCH -p short                 # use short partition
#SBATCH --mem 800M               # RAM usage
#SBATCH --array=0-134            # Run array for indexes 0-n
#SBATCH -o slurm_output/slurm-%A_%a.out

module load conda/miniforge3/24.11.3-0
conda activate /path/to/.conda/envs/your_env

# Run with '${mult}' etc. as the variable names

# Jobs run with progressive time rescaling

# use array = 0-1002
python computational_model_script_4params.py --jobgroup=$SLURM_ARRAY_TASK_ID  --starting_counts='subonly' --jobfile='fullgrid_jobs_prospeedup.pickle' --speedup=3 --A_bins=200 --B_bins=200 --dir=grid_4param/

# use array=0-200
#python computational_model_script_4params.py --jobgroup=$SLURM_ARRAY_TASK_ID  --starting_counts='subonly' --jobfile='3Mgrid_jobs_prospeedup.pickle' --speedup=3 --A_bins=200 --B_bins=200 --interp  --dir=grid_3M_interp/

# use array = 0-999
python3 computational_model_script_4params.py --jobgroup=$SLURM_ARRAY_TASK_ID  --starting_counts='subonly' --rates_function='log' --jobfile='grid_group_jobs_log.pickle' --speedup=3 --A_bins=200 --B_bins=200 --dir=grid_4param_log/

python3 computational_model_script_4params.py --jobgroup=$SLURM_ARRAY_TASK_ID  --jobfile='grid_group_jobs_xaxis.pickle' --starting_counts='subonly' --rates_function='powerlaw_lambda' --speedup=3 --A_bins=200 --B_bins=200 --dir=grid_4param_xaxis/


# Fast-running jobs for constant time rescaling
# use array=0-47
python computational_model_script_4params.py --jobfile='3Mgrid_jobs_constantspeedup_short.pickle' --jobgroup=$SLURM_ARRAY_TASK_ID  --starting_counts='subonly' --speedup=1 --ceiling=0.1 --rounds=9 --A_bins=200 --B_bins=200  --dir=grid_sparse_constantspeedup/ --constantspeedup
python computational_model_script_4params.py --jobfile='3Mgrid_jobs_constantspeedup_short.pickle' --jobgroup=$SLURM_ARRAY_TASK_ID  --starting_counts='uniform' --speedup=1 --ceiling=0.1  --rounds=9 --A_bins=200 --B_bins=200  --dir=grid_sparse_constantspeedup_uniform/ --constantspeedup

#!/bin/bash

# Note: below script is used to run computational model data points in parallel, and is designed to run on a Slurm-based cluster
# This is an example only. Modifications may be required depending on your cluster specifications and conda installation.
# Below job examples require significant computational time (due to constant time rescaling).

#SBATCH -c 20                    # request 20 cores
#SBATCH -t 3-00:00               # 3 days run time (12hr max for short, 5 days max for medium)
#SBATCH -p medium                # use short or medium partition
#SBATCH --mem 800M               # RAM usage
#SBATCH --array=0-14             # Run array for indexes 1-n (set according to job list file)
#SBATCH -o slurm_output/slurm-%A_%a.out

module load conda/miniforge3/24.11.3-0
conda activate /path/to/.conda/envs/your_env


# use --array=0-102
#python computational_model_script_4params.py --jobfile='3Mgrid_jobs_constantspeedup_medium.pickle' --jobgroup=$SLURM_ARRAY_TASK_ID  --starting_counts='subonly' --speedup=1 --ceiling=0.1 --rounds=9 --A_bins=200 --B_bins=200  --dir=grid_sparse_constantspeedup/ --constantspeedup
#python computational_model_script_4params.py --jobfile='3Mgrid_jobs_constantspeedup_medium.pickle' --jobgroup=$SLURM_ARRAY_TASK_ID  --starting_counts='uniform' --speedup=1 --ceiling=0.1 --rounds=9 --A_bins=200 --B_bins=200  --dir=grid_sparse_constantspeedup_uniform/ --constantspeedup

# use --array=0-14
python computational_model_script_4params.py --jobfile='95hdr_all.pickle' --jobgroup=$SLURM_ARRAY_TASK_ID  --starting_counts='subonly' --speedup=1 --ceiling=0.1 --rounds=9 --A_bins=200 --B_bins=200  --recnum=1001 --dir=grid_constantspeedup_top95c/subonly_rerun/ --constantspeedup

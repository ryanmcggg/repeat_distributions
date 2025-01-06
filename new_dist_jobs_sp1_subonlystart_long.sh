#!/bin/bash
#SBATCH -c 1                     # request one core
#SBATCH -t 20-00:00                # 3:00 hours run time
#SBATCH -p long                 # use short partition
#SBATCH --mem 1000M               # use 100Mb
#SBATCH --array=1-4249           # Run array for indexes 1-3108
#SBATCH -o slurm_output/slurm-%A_%a.out

source activate repeat_dist_env

# required parameters: '--mult' '--exp' '--con'


# Specify the path to the config file
config=grid_A_L9_subonlystart_long/jobs_sp1.csv


# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
mult=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
exp=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
con=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)

# Run with '${mult}' etc. as the variable names
python simulation_script.py --exp=${exp} --con=${con} --mult=${mult} --subonly_start --speedup=1 --rounds=8 --A_bins=100 --dir=grid_A_L9_subonlystart_long/

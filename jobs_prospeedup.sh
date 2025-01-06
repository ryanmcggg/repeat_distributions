#!/bin/bash
#SBATCH -c 1                     # request one core
#SBATCH -t 0-08:00                # 0:06 hours run time
#SBATCH -p short                 # use short partition
#SBATCH --mem 100M               # use 100Mb
#SBATCH --array=1-4056           # Run array for indexes 1-18207
#SBATCH -o slurm_output/slurm-%A_%a.out

source activate repeat_dist_env

# required parameters: '--mult' '--exp' '--con'


# Specify the path to the config file
#config=grid_A_L9_subonlystart_prospeedup/jobs.csv
#config=grid_A_L9_subonlystart_prospeedup/jobs2.csv
#config=grid_A_L9_subonlystart_prospeedup/jobs_incomplete.csv
config=grid_A_L9_subonlystart_prospeedup/jobs_short.csv

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
mult=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
exp=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
con=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)


# Run with '${mult}' etc. as the variable names
#python simulation_script_speedup.py --exp=${exp} --con=${con} --mult=${mult}  --subonly_start --speedup=5 --rounds=1 --A_bins=200 --B_bins=200 --dir=grid_A_L9_subonlystart_prospeedup_200/

# grid with stochastics
python simulation_script_speedup.py --exp=${exp} --con=${con} --mult=${mult}  --stochastics=1 --subonly_start --speedup=5 --rounds=1 --A_bins=200 --B_bins=200 --dir=grid_A_L9_subonlystart_prospeedup_200_stoch/

# sample with stochastics
#python simulation_script_speedup.py --exp=1.5 --con=1.8 --mult=2  --stochastics=$SLURM_ARRAY_TASK_ID --subonly_start --speedup=5 --rounds=1 --A_bins=200 --B_bins=200 --dir=grid_A_L9_subonlystart_prospeedup_200/stoch/


# poisson for single condition
#python simulation_script_speedup.py --exp=1.5 --con=1.8 --mult=2 --subonly_start --speedup=5 --rounds=1 --boot=$SLURM_ARRAY_TASK_ID --A_bins=200 --B_bins=200 --dir=grid_A_L9_subonlystart_prospeedup_200/poisson/



# testing
#python simulation_script_speedup.py --exp=1.5 --con=1.8 --mult=2 --subonly_start --speedup=5 --rounds=1 --dir=grid_A_L9_subonlystart_prospeedup/

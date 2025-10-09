#!/bin/bash
#SBATCH -c 1                     # request 2 cores (max 20 for short)
#SBATCH -t 0-01:00                  # 1 hour run time (12hr max for short)
#SBATCH -p short                  # use short partition
#SBATCH --mem 16G               # use 32Gb
#SBATCH --array=1-1             # Run array for indexes 1-n (Slurm doesn't accept more than 10000 jobs per batch)
#SBATCH -o slurm_output/slurm-%A_%a.out

module load conda/miniforge3/24.11.3-0
conda activate /home/rjm44/.conda/envs/rep_env_3

#python gather_grid.py  --dir=grid_4param_interp_uniform/ --name=grid_4param_interp_uniform --write=pickle
#python gather_grid.py  --dir=grid_4param_interp_uniform/ --name=grid_4param_interp_uniform_final --write=pickle --final --B_dist
#python gather_grid.py  --dir=grid_4param_interp/ --name=grid_4param_interp_final --write=pickle --final --B_dist

#python gather_grid.py  --dir=grid_4param_log_v2/ --name=grid_4param_log_final --write=pickle --final --B_dist

#python gather_grid.py  --dir=grid_4param_xaxis/ --name=grid_4param_xaxis_final --write=pickle --final --B_dist --sep='tE|tC|LE|LC' --check_complete


#python gather_grid.py  --dir=grid_4param_xaxis/ --name=grid_4param_xaxis --write=pickle --sep='tE|tC|LE|LC' --check_complete



#python gather_grid.py  --dir=grid_4param_v5/ --name=grid_4param_PL_final --write=pickle --final --B_dist
#python gather_grid.py  --dir=grid_4param_log_v3/ --name=grid_4param_log_final --write=pickle --final --B_dist

#python gather_grid.py  --dir=grid_4param_v5/ --name=grid_4param_PL --write=pickle
#python gather_grid.py  --dir=grid_4param_log_v3/ --name=grid_4param_log --write=pickle


#python gather_grid.py  --dir=grid_sparse_constantspeedup_uniform/ --name=grid_4param_PL_constant_uniform_final --write=pickle --final
#python gather_grid.py  --dir=grid_sparse_constantspeedup_uniform/ --name=grid_4param_PL_constant_uniform --write=pickle

#python gather_grid.py  --dir=grid_4param_xaxis_v2/ --name=grid_4param_xaxis_v2 --final --write=pickle --sep='tE|tC|LE|LC' --check_complete='grid_group_jobs_xaxis_v2.pickle'

#python gather_grid.py  --dir=grid_4param_v6/ --name=grid_newpro_PL_final --write=pickle --final --B_dist
#python gather_grid.py  --dir=grid_4param_v6/ --name=grid_newpro_PL --write=pickle


#python gather_grid.py  --dir=grid_top95c_pro_stoch/ --name=grid_top95c_stoch_2 --write=pickle
#python gather_grid.py  --dir=grid_top95c_pro_stoch/3/ --name=grid_top95c_stoch_3 --write=pickle
#python gather_grid.py  --dir=grid_top95c_pro_stoch/4/ --name=grid_top95c_stoch_4 --write=pickle
#python gather_grid.py  --dir=grid_top95c_pro_stoch/5/ --name=grid_top95c_stoch_5 --write=pickle
#python gather_grid.py  --dir=grid_top95c_pro_stoch/6/ --name=grid_top95c_stoch_6 --write=pickle
#python gather_grid.py  --dir=grid_top95c_pro_stoch/7/ --name=grid_top95c_stoch_7 --write=pickle
#python gather_grid.py  --dir=grid_top95c_pro_stoch/8/ --name=grid_top95c_stoch_8 --write=pickle
#python gather_grid.py  --dir=grid_top95c_pro_stoch/9/ --name=grid_top95c_stoch_9 --write=pickle
#python gather_grid.py  --dir=grid_top95c_pro_stoch/10/ --name=grid_top95c_stoch_10 --write=pickle

#python gather_grid.py  --dir=grid_3M_interp/ --name=grid_3M_PL_subonly_interp_final --write=pickle --final --B_dist

#python gather_grid.py  --dir=grid_sparse_constantspeedup/ --name=grid_4param_PL_constant_subonly_final --write=pickle --final --col_list='[8,10,12,14]'
#python gather_grid.py  --dir=grid_sparse_constantspeedup/ --name=grid_4param_PL_constant_subonly --write=pickle --col_list='[8,10,12,14]'


#python gather_grid.py  --dir=grid_3M_stoch/ --name=grid_3M_PL_subonly_stoch_final --write=pickle --final --B_dist

#python gather_grid.py  --dir=grid_sparse_constantspeedup_uniform/ --name=test --write=pickle --final --col_list='[8,10,12,14]' --check_complete='3Mgrid_jobs_constantspeedup_short.pickle'
#python gather_grid.py  --dir=grid_sparse_constantspeedup_uniform/ --name=grid_4param_PL_constant_uniform_final --write=pickle --final --col_list='[8,10,12,14]'
#python gather_grid.py  --dir=grid_sparse_constantspeedup_uniform/ --name=grid_4param_PL_constant_uniform --write=pickle --col_list='[8,10,12,14]'


#python gather_grid.py  --dir=grid_4param_xaxis_v3/ --name=grid_4param_xaxis_v3 --final --write=pickle --sep='tE|tC|LE|LC' --check_complete='grid_group_jobs_xaxis_v2.pickle'
#python gather_grid.py  --dir=grid_4param_xaxis_v3/ --name=grid_4param_xaxis_v3 --final --write=pickle --sep='tE|tC|LE|LC'

#python gather_grid.py  --dir=grid_4param_log_v4/ --name=grid_4param_log_v4_final --final --write=pickle --check_complete='fullgrid_jobs_prospeedup.pickle'


#python gather_grid.py  --dir=grid_constantspeedup_top95c/subonly_rerun/ --name=top95constant_subonlyrerun --write=pickle --col_list='[8,10,12,14]' --check_complete='95hdr_all.pickle'
python gather_grid.py  --dir=grid_constantspeedup_top95c/subonly_rerun/ --name=top95constant_subonlyrerun_final --final --write=pickle --col_list='[8,10,12,14]' --check_complete='95hdr_all.pickle'

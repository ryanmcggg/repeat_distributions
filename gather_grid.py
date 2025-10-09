#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import os
import pickle

def gather_files(folder, run_name, write='out', final_only = False, check_complete = None, B_dist = False, separator = 'tE|tC|iE|iC', col_list = [6,8,10,12]):
    grid_files = os.listdir(folder)
    grid_files_A = sorted([file for file in grid_files if 'Adist' in file])
    filename_start = pd.Series(grid_files_A).str.split(separator)[0][0]
    filename_end = '_' + pd.Series(grid_files_A).str.split(separator)[0][4].split('_')[1]
    parameter_info = pd.Series(grid_files_A).str.split(separator + '|_', expand=True)[col_list].astype(float)
    parameter_info.index = grid_files_A
    parameter_info.columns = ['exp_p', 'con_p', 'exp_i', 'con_i']
    parameter_info = parameter_info[['exp_i', 'con_i', 'exp_p', 'con_p']]
    parameter_info['exp_i'] = parameter_info['exp_i'].astype(int); parameter_info['con_i'] = parameter_info['con_i'].astype(int)
    parameter_info = parameter_info.sort_index().drop_duplicates(keep = 'first') # if any duplicates exist, keep one with lowest speedup
    #parameters_inference_cols = parameter_info.groupby(['exp_i', 'con_i', 'exp_p', 'con_p'])['con_p'].count().index
    parameters_inference_cols = pd.Series([(a,b,c,d) for a,b,c,d in zip(parameter_info['exp_i'], parameter_info['con_i'], parameter_info['exp_p'], parameter_info['con_p'])], index = parameter_info.index)
#    parameter_info.to_pickle('grid_present_' + str(run_name) + '.pickle')

    if check_complete is not None:
        # check for missing files
        jobs = pd.read_pickle(check_complete)
        jobs_incomplete = jobs.loc[jobs.isin(parameters_inference_cols) == False].reset_index(drop=True)
        if len(jobs_incomplete) > 0:
            print('missing ' + str(len(jobs_incomplete)) + ' files')
            jobs_incomplete.to_pickle('jobs_incomplete.pickle')
            return jobs_incomplete
        else:
            print('grid is complete')

    if write != False:
        inference_grid = dict()
        inference_grid_Bdist = dict()
        for file in parameter_info.index:
            try:
                if final_only == False:
                    inference_grid[file] = pd.read_pickle(folder + file)
                else:
                    inference_grid[file] = pd.read_pickle(folder + file).T.tail(1).T
            except pickle.UnpicklingError:
                os.remove(folder + file)
                print('deleted: ' + str(file))
            if B_dist == True:
                inference_grid_Bdist[file] = pd.read_pickle(folder + file.replace('Adist', 'Bdist')).T.tail(1).T
        inference_grid = pd.concat(inference_grid, axis=1).sort_index(axis=1)
        cols = parameter_info.reindex(inference_grid.columns.get_level_values(0))
        cols['round'] = inference_grid.columns.get_level_values(1)
        cols = pd.MultiIndex.from_frame(cols)
        inference_grid.columns = cols
        inference_grid.index +=1
        if write == 'pickle':
            inference_grid.to_pickle('completed_grids/' + run_name + '.pickle')
        if write == 'csv':
            inference_grid.to_csv('completed_grids/' + run_name + '.csv', sep = '\t')

        if B_dist == True:
            inference_grid_Bdist = pd.concat(inference_grid_Bdist, axis=1).sort_index(axis=1)
            inference_grid_Bdist.columns = cols
            inference_grid_Bdist.index +=1
            if write == 'pickle':
                inference_grid_Bdist.to_pickle('completed_grids/' + run_name + '_Bdist.pickle')
            if write == 'csv':
                inference_grid_Bdist.to_csv('completed_grids/' + run_name + '_Bdist.csv', sep = '\t')

        if write == 'out':
            return inference_grid, inference_grid_Bdist


import argparse
parser = argparse.ArgumentParser(description='repeat distribution simulation - gather data')
parser.add_argument('--dir', action="store", dest='dir', default = 'simulations/grid_output/', type=str)
parser.add_argument('--name', action="store", dest='name', default = 'grid_name', type=str)
parser.add_argument('--write', action="store", dest='write', default = 'out', type=str)
parser.add_argument('--final', default=False, action="store_true")
parser.add_argument('--check_complete', default=None, action="store", type = str)
parser.add_argument('--B_dist', default=False, action="store_true")
parser.add_argument('--sep', action="store", dest='sep', default = 'tE|tC|iE|iC', type=str)
parser.add_argument('--col_list', action="store", dest='col_list', default = '[6,8,10,12]', type=str)

args = parser.parse_args()
args.col_list = eval(args.col_list)
gather_files(args.dir, args.name, args.write, args.final, args.check_complete, args.B_dist, args.sep, args.col_list)

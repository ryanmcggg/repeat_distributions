#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import pickle
import os
import multiprocessing

import warnings
warnings.filterwarnings('ignore', '(overflow|invalid)')

# #### Load data
CHM13_counts = pd.read_pickle('repeat_distributions/CHM13_counts.pickle')
CHM13_counts_B = pd.read_pickle('repeat_distributions/CHM13_counts_Bdist.pickle')
random_counts = pd.read_pickle('repeat_distributions/random_counts.pickle')
random_counts_B = pd.read_pickle('repeat_distributions/random_counts_Bdist.pickle')

denovo_exp_rate = pd.read_pickle('denovo/denovo_exp_rate.pickle')
denovo_con_rate = pd.read_pickle('denovo/denovo_con_rate.pickle')
denovo_nonexp_rate = pd.read_pickle('denovo/denovo_nonexp_rate.pickle')

denovo_exp_rate_poisson = pd.read_pickle('denovo/denovo_exp_rate_poisson.pickle')
denovo_con_rate_poisson = pd.read_pickle('denovo/denovo_con_rate_poisson.pickle')
denovo_nonexp_rate_poisson = pd.read_pickle('denovo/denovo_nonexp_rate_poisson.pickle')

denovo_substitution_context_rate = pd.read_pickle('denovo/denovo_mut_freq_triplets.pickle')
denovo_substitution_context_rate_poisson = pd.read_pickle('denovo/denovo_mut_freq_triplets_poisson.pickle')

rates_mu_nu = pd.read_pickle('denovo/denovo_mut_freq_mu_nu_remake.pickle')

intercept_list = dict()
intercept_list['A'] = [denovo_con_rate['A'][8]] + [denovo_exp_rate['A'][8] * (denovo_exp_rate['A'][8]/ denovo_con_rate['A'][8])**x for x in range(7)]
intercept_list['AGC'] = [denovo_con_rate['AGC'][4]] + [denovo_exp_rate['AGC'][4] * (denovo_exp_rate['AGC'][4]/ denovo_con_rate['AGC'][4])**x for x in range(7)]

subonly_counts = pd.read_pickle('repeat_distributions/subonly_counts_remake.pickle')
counts_uniform = pd.read_pickle('repeat_distributions/subonly_counts_uniformrange.pickle')


# evolve function
def mut_evolve_dist_AB(A_count_input, B_count_input, starting_conditions, boot = None, input_nuc = 'A', mut = True, mutonly = False, speedup_multiplier = 1, output_components = False, stochastics = None, reflective = True, boundary_count = 1000, first_bin = 0):
    exp_rate_A_AA, con_rate_A_AA, nonexp_rate_A_AB, B_indel_rates = starting_conditions
    A_count_output = A_count_input.copy(); B_count_output = B_count_input.copy()
    A_bins = len(A_count_input); B_bins = len(B_count_input)
    boundary_flag = False
    A_length_array = np.array(range(1,A_bins+3))
    A_length_array_bases = np.array(range(1,A_bins+3)) * len(input_nuc) ### including motif length
    B_length_array = np.array(range(1,B_bins+3))
    B_length_array_bases = np.array(range(1,B_bins+3)) * len(input_nuc) ### including motif length
    A_count_input = np.insert(A_count_input, A_bins, [0,0])
    B_count_input = np.insert(B_count_input, B_bins, [0,0])
    #A_count_input = A_count_input.astype('int64')
    #B_count_input = B_count_input.astype('int64')

    if boot is None:
        denovo_sub = denovo_substitution_context_rate.loc[input_nuc]
    else:
        denovo_sub = denovo_substitution_context_rate_poisson.loc[boot].loc[input_nuc]

    # distribution info
    total_A_bases = (A_count_input[:A_bins] * A_length_array_bases[:A_bins]).sum()
    total_B_bases = (B_count_input[:B_bins] * B_length_array_bases[:B_bins]).sum()
    # try to quit before encountering overflow errors
    if (total_A_bases > 1e300) | (total_B_bases > 1e300) | (total_A_bases < 1e-300) | (total_B_bases < 1e-300):
        return A_count_output, B_count_output, True, boundary_flag

    A_nonflank_base_portion = (A_count_input[2:A_bins+2] * A_length_array_bases[:A_bins]).sum() / total_A_bases
    A_flank_base_portion = (A_count_input[1:A_bins] * 2 * len(input_nuc)).sum() / total_A_bases ### including motif length
    B_L1_base_portion = ((B_count_input[0] * len(input_nuc)) / (B_count_input[:B_bins]* B_length_array_bases[:B_bins]).sum()) ### including motif length
    B_nonflank_base_portion = (B_count_input[2:B_bins+2] * B_length_array_bases[:B_bins]).sum() / total_B_bases  ### include portion of triplets 1nt away???
    B_flank_base_portion = (B_count_input[1:B_bins] * 2 * len(input_nuc)).sum() / total_B_bases ### including motif length

    total_A_change_in = np.array([0.0]*A_bins); total_B_change_in = np.array([0.0]*B_bins)
    total_A_change_out = np.array([0.0]*A_bins); total_B_change_out = np.array([0.0]*B_bins)

    if mut == True:
        # A>B which adds to the A count locally. add these to A
        A_mut_in_local_A_B = 2 * len(input_nuc) * denovo_sub['Acontraction'] * A_count_input[1:]
        A_mut_out_local_A_B = -2 * len(input_nuc) * denovo_sub['Acontraction'] * A_count_input
        A_mut_out_local_A_B[0] = -1 * len(input_nuc) * A_count_input[0] * denovo_sub['A10']
        #A_mut_in_local_A_B = A_mut_out_local_A_B[1:]

        # total number of A>B fission events
        A_mut_out_fission = np.insert((-denovo_sub['Afission'] * A_count_input[2:] * A_length_array_bases[:A_bins]), 0, [0, 0]) # used to subtract from A_count, starting from L=3 (with 0 for L=1,2)
        # each fission creates 2 As. add these to A
        A_mut_in_fission = ((2/A_length_array[:A_bins]) * -A_mut_out_fission[2:A_bins+2])[::-1].cumsum()[::-1]
 
        # B>A which adds to the A count locally (which must come from B_L>1)
        # A from B>A leaving the -1 bin
        A_len_freq = (A_count_input / A_count_input.sum())[:A_bins]
        A_mut_out_local_B_A = -denovo_sub['Aexpansion'] * B_flank_base_portion * total_B_bases * A_len_freq
        # B>A creating A_L=1 from B_L>2
        B_A_into_L1 = total_B_bases * B_nonflank_base_portion * denovo_sub['A01']
        A_mut_in_local_B_A = np.insert(-A_mut_out_local_B_A, 0, B_A_into_L1)
        
        # fusion process for A
        A_len_freq = (A_count_input / A_count_input.sum())[:A_bins]
        A_fusion_freq_in = np.bincount((np.add.outer(A_length_array[:A_bins], A_length_array[:A_bins])+1).ravel(), weights = np.outer(A_len_freq, A_len_freq).ravel())[1:]
        A_mut_in_fusion_A_B = A_fusion_freq_in * denovo_sub['Afusion'] * B_L1_base_portion * total_B_bases
        A_mut_out_fusion_A_B = (-2) *A_len_freq * denovo_sub['Afusion'] * B_L1_base_portion * total_B_bases
        
        
        # total B>A
        # B>A which adds to the B count locally. add these to B
        B_mut_in_local_B_A = 2 * len(input_nuc) * denovo_sub['Aexpansion'] * B_count_input[1:]
        B_mut_out_local_B_A = -2 * len(input_nuc) * denovo_sub['Aexpansion'] * B_count_input
        B_mut_out_local_B_A[0] = -1 * B_L1_base_portion * total_B_bases * denovo_sub['Afusion']

        # total number of B>A fission events
        B_mut_out_fission = np.insert((-denovo_sub['A01'] * B_count_input[2:] * B_length_array_bases[:B_bins]), 0, [0, 0]) # used to subtract from B_count, starting from L=3 (with 0 for L=1,2)
        # each fission creates 2 Bs. add these to B
        B_mut_in_fission = ((2/B_length_array[:B_bins]) * -B_mut_out_fission[2:B_bins+2])[::-1].cumsum()[::-1]

        # A>B which adds to the B count locally (which must come from A_L>1)
        # B from A>B leaving the -1 bin
        B_len_freq = (B_count_input / B_count_input.sum())[:B_bins]
        B_mut_out_local_A_B = -denovo_sub['Acontraction'] * A_flank_base_portion * total_A_bases * B_len_freq
        # A>B creating B_L=1 from A_L>2
        A_B_into_L1 = total_A_bases * A_nonflank_base_portion * denovo_sub['Afission']
        B_mut_in_local_A_B = np.insert(-B_mut_out_local_A_B, 0, A_B_into_L1)
        
        # fusion process for B
        B_len_freq = (B_count_input / B_count_input.sum())[:B_bins]
        B_fusion_freq_in = np.bincount((np.add.outer(B_length_array[:B_bins], B_length_array[:B_bins])+1).ravel(), weights = np.outer(B_len_freq, B_len_freq).ravel())[1:]
        B_mut_in_fusion_B_A = B_fusion_freq_in * denovo_sub['A10'] * A_count_input[0] * len(input_nuc)
        B_mut_out_fusion_B_A = (-2) * B_len_freq * denovo_sub['A10'] * A_count_input[0] * len(input_nuc)

        # update counts for next round (with absorbing boundary)
        total_A_change_in += A_mut_in_local_A_B[:A_bins] + A_mut_in_local_B_A[:A_bins] + A_mut_in_fission[:A_bins] + A_mut_in_fusion_A_B[:A_bins]
        total_B_change_in += B_mut_in_local_B_A[:B_bins] + B_mut_in_local_A_B[:B_bins] + B_mut_in_fission[:B_bins] + B_mut_in_fusion_B_A[:B_bins]
        total_A_change_out += A_mut_out_local_A_B[:A_bins] + A_mut_out_local_B_A[:A_bins] + A_mut_out_fission[:A_bins] + A_mut_out_fusion_A_B[:A_bins]
        total_B_change_out += B_mut_out_local_B_A[:B_bins] + B_mut_out_local_A_B[:B_bins] + B_mut_out_fission[:B_bins] + B_mut_out_fusion_B_A[:B_bins]

        # apply reflecting boundary
        if reflective == True:
            total_A_change_in[A_bins-1] += A_mut_in_local_A_B[A_bins:].sum() + A_mut_in_local_B_A[A_bins:].sum() + A_mut_in_fission[A_bins:].sum() + A_mut_in_fusion_A_B[A_bins:].sum()
            total_B_change_in[B_bins-1] += B_mut_in_local_B_A[B_bins:].sum() + B_mut_in_local_A_B[B_bins:].sum() + B_mut_in_fission[B_bins:].sum() + B_mut_in_fusion_B_A[B_bins:].sum()
            total_A_change_out[A_bins-1] += A_mut_out_local_A_B[A_bins:].sum() + A_mut_out_local_B_A[A_bins:].sum() + A_mut_out_fission[A_bins:].sum() + A_mut_out_fusion_A_B[A_bins:].sum()
            total_B_change_out[B_bins-1] += B_mut_out_local_B_A[B_bins:].sum() + B_mut_out_local_A_B[B_bins:].sum() + B_mut_out_fission[B_bins:].sum() + B_mut_out_fusion_B_A[B_bins:].sum()
           
    if mutonly == False:
        # A expansions in and out
        A_exp_out = A_count_input[:A_bins] * -exp_rate_A_AA[:A_bins]
        A_exp_in = np.insert(-A_exp_out, 0, B_indel_rates[0]*total_B_bases)

        # A contractions in and out
        A_con_out = A_count_input[:A_bins+1] * -con_rate_A_AA[:A_bins+2]
        A_con_in = -A_con_out[1:]

        # A fusions from B1>B0 deletions
        if (mut != True):
            A_len_freq = (A_count_input / A_count_input.sum())[:A_bins]
            A_fusion_freq_in = np.bincount((np.add.outer(A_length_array[:A_bins], A_length_array[:A_bins])+1).ravel(), weights = np.outer(A_len_freq, A_len_freq).ravel())[1:]
        A_mut_in_fusion_Bdel = A_fusion_freq_in[1:A_bins+1] * B_indel_rates[1] * B_L1_base_portion * total_B_bases
        A_mut_out_fusion_Bdel = (-2) *A_len_freq * B_indel_rates[1] * B_L1_base_portion * total_B_bases

        # A fission events from insertions
        A_nonexp_out_fission = -A_count_input * nonexp_rate_A_AB # used to calculate fission_in, starting with L=2 going to 2x L=1
        # each fission creates 2 As. add these to A
        A_nonexp_in_fission = ((2/A_length_array[:A_bins]) * -A_nonexp_out_fission[1:A_bins+1])[::-1].cumsum()[::-1]

        # B expansions in and out
        B_exp_out = B_count_input[:B_bins] * -B_indel_rates[2] * B_length_array[:B_bins] # B>BB rates are flat, per base
        B_exp_in = np.insert(-B_exp_out, 0, A_nonexp_out_fission.sum())

        # B contractions in and out
        B_con_out = B_count_input[:B_bins+1] * -B_indel_rates[1] * B_length_array[:B_bins+1] # B>_ rates are flat, per base
        B_con_in = -B_con_out[1:]

        # B fusions from A1>A0 deletions
        if (mut != True):
            B_len_freq = (B_count_input / B_count_input.sum())[:B_bins]
            B_fusion_freq_in = np.bincount((np.add.outer(B_length_array[:B_bins], B_length_array[:B_bins])+1).ravel(), weights = np.outer(B_len_freq, B_len_freq).ravel())[1:]
        B_mut_in_fusion_Adel = B_fusion_freq_in[1:B_bins+1] * A_count_input[0] * con_rate_A_AA[0]
        B_mut_out_fusion_Adel = 2 *B_len_freq * A_count_input[0] * -con_rate_A_AA[0]
                    
        # B fission events from insertions
        B_nonexp_out_fission = -B_count_input * B_indel_rates[0] * B_length_array # used to calculate fission_in, starting with L=2 going to 2x L=1
        # each fission creates 2 Bs. add these to B
        B_nonexp_in_fission = ((2/B_length_array[:B_bins]) * -B_nonexp_out_fission[1:B_bins+1])[::-1].cumsum()[::-1]

       # update counts for next round (with absorbing boundary)
        total_A_change_in += A_exp_in[:A_bins] + A_con_in[:A_bins] + A_mut_in_fusion_Bdel[:A_bins] + A_nonexp_in_fission[:A_bins]
        total_B_change_in += B_exp_in[:B_bins] + B_con_in[:B_bins] + B_mut_in_fusion_Adel[:B_bins] + B_nonexp_in_fission[:B_bins]
        total_A_change_out += A_exp_out[:A_bins] + A_con_out[:A_bins] + A_mut_out_fusion_Bdel[:A_bins] + A_nonexp_out_fission[:A_bins]
        total_B_change_out += B_exp_out[:B_bins] + B_con_out[:B_bins] + B_mut_out_fusion_Adel[:B_bins] + B_nonexp_out_fission[:B_bins]

        # apply reflecting boundary
        if reflective == True:
            total_A_change_in[A_bins-1] += A_exp_in[A_bins:].sum() + A_con_in[A_bins:].sum() + A_mut_in_fusion_Bdel[A_bins:].sum() + A_nonexp_in_fission[A_bins:].sum()
            total_B_change_in[B_bins-1] += B_exp_in[B_bins:].sum() + B_con_in[B_bins:].sum() + B_mut_in_fusion_Adel[B_bins:].sum() + B_nonexp_in_fission[B_bins:].sum()
            total_A_change_out[A_bins-1] += A_exp_out[A_bins:].sum() + A_con_out[A_bins:].sum() + A_mut_out_fusion_Bdel[A_bins:].sum() + A_nonexp_out_fission[A_bins:].sum()
            total_B_change_out[B_bins-1] += B_exp_out[B_bins:].sum() + B_con_out[B_bins:].sum() + B_mut_out_fusion_Adel[B_bins:].sum() + B_nonexp_out_fission[B_bins:].sum()

    # flag to stop the simulation if more repeats are removed from a bin than exist in that bin (excluding the last 10 noisy bins)
    flag = ((np.abs(total_A_change_out[:A_bins-10]) * speedup_multiplier > A_count_output[:A_bins-10]).sum()) > 0

    # apply speedup
    total_A_change_in *= speedup_multiplier; total_A_change_out *= speedup_multiplier
    total_B_change_in *= speedup_multiplier; total_B_change_out *= speedup_multiplier
    
    if stochastics is not None:
        # the sum of poisson random variables is poisson-distributed. not necessary to run n poisson samples
        total_A_change_in = np.random.poisson(total_A_change_in.clip(0))
        total_A_change_out = -1 * np.random.poisson(np.abs(total_A_change_out.clip(max=0)))
        total_B_change_in = np.random.poisson(total_B_change_in.clip(0))
        total_B_change_out = -1 * np.random.poisson(np.abs(total_B_change_out.clip(max=0)))   
    
    total_A_change = total_A_change_in + total_A_change_out
    total_B_change = total_B_change_in + total_B_change_out

    if first_bin != 0:
        total_A_change[:first_bin] = 0; total_B_change[:first_bin] = 0
        
    # update counts for next round
    A_count_output = A_count_output[:A_bins] + total_A_change[:A_bins]
    B_count_output = B_count_output[:B_bins] + total_B_change[:B_bins]

    # remove negative values
    A_count_output[A_count_output <0] = 0            
    B_count_output[B_count_output <0] = 0

    boundary_flag = ((A_count_output[A_bins-1] > boundary_count) | (B_count_output[B_bins-1] > boundary_count))
    
    if output_components == True:
        if mutonly == False:
            return  A_mut_in_local_A_B[:A_bins], A_mut_out_local_A_B[:A_bins], A_mut_in_local_B_A[:A_bins], A_mut_out_local_B_A[:A_bins], A_mut_in_fission[:A_bins], A_mut_out_fission[:A_bins], A_mut_in_fusion_A_B[:A_bins], A_mut_out_fusion_A_B[:A_bins], A_exp_in[:A_bins], A_exp_out[:A_bins], A_con_in[:A_bins], A_con_out[:A_bins], A_mut_in_fusion_Bdel[:A_bins], A_mut_out_fusion_Bdel[:A_bins], A_nonexp_in_fission[:A_bins], A_nonexp_out_fission[:A_bins]
        else:
            return  A_mut_in_local_A_B[:A_bins], A_mut_out_local_A_B[:A_bins], A_mut_in_local_B_A[:A_bins], A_mut_out_local_B_A[:A_bins], A_mut_in_fission[:A_bins], A_mut_out_fission[:A_bins], A_mut_in_fusion_A_B[:A_bins], A_mut_out_fusion_A_B[:A_bins]
    else:
        return A_count_output, B_count_output, flag, boundary_flag


def pin_power_law(power, pin_rate, pin_len=9, start_len=1, end_len=200):
    denom = (pin_len**power) / pin_rate
    return pd.Series([i**power for i in range(start_len, end_len+1)], index = list(range(start_len,end_len+1))) / denom

def intercept_then_powerlaw(exp_power, con_power, exp_int, con_int, pin_len = 9, A_bins = 200, boot = None, motif = 'A', interp = False, nonexp_factor = 0.01, intercept_group = 'motif'):
    if intercept_group == 'motif':
        intercept_group = intercept_list[motif]
    if boot is None:
        bootname = ''
        denovo_exp_rate_current = denovo_exp_rate[motif].copy()
        denovo_con_rate_current = denovo_con_rate[motif].copy()
        denovo_nonexp_rate_current = denovo_nonexp_rate[motif].copy()
    else:
        bootname = '_boot'+str(boot)
        denovo_exp_rate_current = denovo_exp_rate_poisson[motif][boot].copy()
        denovo_con_rate_current = denovo_con_rate_poisson[motif][boot].copy()
        denovo_nonexp_rate_current = denovo_nonexp_rate_poisson[motif][boot].copy()
    exp = pd.concat([denovo_exp_rate_current.reindex(range(pin_len)), pin_power_law(exp_power, intercept_group[exp_int], start_len=pin_len, end_len=A_bins+3)])
    con = pd.concat([denovo_con_rate_current.reindex(range(pin_len)), pin_power_law(con_power, intercept_group[con_int], start_len=pin_len, end_len=A_bins+3)])
    nonexp = pd.concat([denovo_nonexp_rate_current.reindex(range(pin_len)), pin_power_law(exp_power, intercept_group[exp_int] * nonexp_factor, start_len=pin_len, end_len=A_bins+3)])
    nonexpname = '_nex' + str(nonexp_factor)
    nonexp.loc[1] = 0

    if interp == True:
        interpname = '_interp'
        exp.loc[8:13] = np.nan
        con.loc[8:13] = np.nan
        nonexp.loc[8:13] = np.nan
        exp = exp.interpolate(method = 'quadratic')
        con = con.interpolate(method = 'quadratic')
        nonexp = nonexp.interpolate(method = 'quadratic')
    else:
        interpname = ''

    name = '_interceptPL_tE' + str(exp_power) + '_tC' + str(con_power) +'_iE' + str(exp_int) + '_iC' + str(con_int) + nonexpname + interpname + bootname
    return name, exp, con, nonexp


#### alternate rate functional forms
def power_law_rate_at_L(power, L, rate = 1e-8, bins = 200):
    return pd.Series([np.exp(np.log(rate) - power * np.log(L) + np.log(length)*power) for length in range(1,bins)], index = range(1,bins))

def power_law_rate_at_L_v2(power, Lam, mu = 4.117e-09, bins = 200):
    return pd.Series([((2*mu)/Lam) * (length/Lam)**power for length in range(1,bins)], index = range(1,bins))

def log_curve_from_L9(power, intercept, motif = 'A', bins = 200):
    return pd.Series([intercept_list[motif][intercept] * (np.log(1+L-8) / np.log(1+9 -8))**power for L in range(9,bins+1)], index = range(9,bins+1))


# setup function
def setup_evolve(exp_power=0, con_power=0, exp_int=0, con_int=0, boot = None, stochastics = None, interp = False, pin_len = 9, nonexp_factor = 0.01, A_bins = 200, B_bins = 200, input_nuc = 'A', mutonly = False, exp_zero = False, con_zero = False, nonexp_zero = False, starting_counts = 'random', ceiling = None, rates_function = 'powerlaw'):
# set up counts
    A_length_array = np.array(range(1,A_bins+1))
    B_length_array = np.array(range(1,A_bins+1))
    if starting_counts == 'random':
        A_count_input = np.nan_to_num(random_counts[len(input_nuc)][input_nuc].reindex(range(1,A_bins+1)).values)
        B_count_input = np.nan_to_num(random_counts_B[input_nuc].reindex(range(1,B_bins+1)).values)
    if starting_counts == 'subonly':
        A_count_input = np.nan_to_num(subonly_counts['A'][input_nuc].reindex(range(1,A_bins+1)).values)
        B_count_input = np.nan_to_num(subonly_counts['B'][input_nuc].reindex(range(1,B_bins+1)).values)
    if starting_counts == 'CHM13':
        A_count_input = np.nan_to_num(CHM13_counts[len(input_nuc)][input_nuc].reindex(range(1,A_bins+1)).values)
        B_count_input = np.nan_to_num(CHM13_counts_B[input_nuc].reindex(range(1,B_bins+1)).values)
    if starting_counts == 'uniform':
        A_count_input = np.nan_to_num(counts_uniform['A'][input_nuc].reindex(range(1,A_bins+1)).values)
        B_count_input = np.nan_to_num(counts_uniform['B'][input_nuc].reindex(range(1,B_bins+1)).values)
    if starting_counts not in ['random', 'subonly', 'CHM13', 'uniform']:
        A_count_input = np.nan_to_num(starting_counts[0].reindex(range(1,A_bins+1)).values)
        B_count_input = np.nan_to_num(starting_counts[1].reindex(range(1,B_bins+1)).values)
# set up rates    
    if mutonly == False:
        if rates_function == 'powerlaw':
            name, exp_rate, con_rate, nonexp_rate = intercept_then_powerlaw(exp_power = exp_power, con_power = con_power, exp_int = exp_int, con_int=con_int, pin_len=pin_len, A_bins = A_bins, boot = boot, motif = input_nuc, interp = interp, nonexp_factor = nonexp_factor)
            B_indel_rate = np.array([exp_rate[0], con_rate[0], nonexp_rate[0]])        
        if rates_function == 'powerlaw_x':
            name = '_PL1e-8atL_tE' + str(exp_power) + '_tC' + str(con_power) +'_LE' + str(exp_int) + '_LC' + str(con_int)
            exp_rate = power_law_rate_at_L(exp_power, exp_int, rate = rates_mu_nu['B>A'], bins = A_bins+3).reindex(range(A_bins+3))
            con_rate = power_law_rate_at_L(con_power, con_int, rate = rates_mu_nu['A>B'], bins = A_bins+3).reindex(range(A_bins+3))
            nonexp_rate = exp_rate * nonexp_factor
            B_indel_rate = np.array([denovo_exp_rate[input_nuc][0], denovo_con_rate[input_nuc][0], denovo_nonexp_rate[input_nuc][0]])
        if rates_function == 'powerlaw_lambda':
            name = '_PLmuatL_tE' + str(exp_power) + '_tC' + str(con_power) +'_LE' + str(exp_int) + '_LC' + str(con_int)
            exp_rate = power_law_rate_at_L_v2(exp_power, exp_int, mu = rates_mu_nu['B>A'], bins = A_bins+3).reindex(range(A_bins+3))
            con_rate = power_law_rate_at_L_v2(con_power, con_int, mu = rates_mu_nu['B>A'], bins = A_bins+3).reindex(range(A_bins+3))
            nonexp_rate = exp_rate * nonexp_factor
            B_indel_rate = np.array([denovo_exp_rate[input_nuc][0], denovo_con_rate[input_nuc][0], denovo_nonexp_rate[input_nuc][0]])
        if rates_function == 'log':
            name = '_logrates_tE' + str(exp_power) + '_tC' + str(con_power) +'_iE' + str(exp_int) + '_iC' + str(con_int)
            exp_rate = pd.concat([denovo_exp_rate['A'].loc[:8], log_curve_from_L9(exp_power, exp_int, input_nuc, A_bins+3)])
            con_rate = pd.concat([denovo_con_rate['A'].loc[:8], log_curve_from_L9(con_power, con_int, input_nuc, A_bins+3)])
            nonexp_rate = pd.concat([denovo_nonexp_rate['A'].loc[:8], log_curve_from_L9(exp_power, exp_int, input_nuc, A_bins+3) * nonexp_factor])
            B_indel_rate = np.array([denovo_exp_rate[input_nuc][0], denovo_con_rate[input_nuc][0], denovo_nonexp_rate[input_nuc][0]])
        if rates_function == 'custom':
            exp_rate, con_rate, name = custom_rates
            exp_rate = exp_rate.reindex(range(A_bins+3))
            con_rate = con_rate.reindex(range(A_bins+3))
            nonexp_rate = exp_rate * nonexp_factor
            B_indel_rate = np.array([denovo_exp_rate[input_nuc][0], denovo_con_rate[input_nuc][0], denovo_nonexp_rate[input_nuc][0]])
        # change rates from per unit to per STR
        exp_rate = exp_rate.values[1:A_bins+1] * A_length_array
        con_rate = con_rate.values[1:A_bins+2] * np.array(range(1,A_bins+2))
        nonexp_rate = nonexp_rate.values[1:A_bins+3] * np.array(range(1,A_bins+3))
        if ceiling != None:
            ceiling_loc = []
            if (exp_rate > ceiling).sum() > 0:
                ceiling_loc.append(pd.Series(exp_rate > ceiling).idxmax())
            if (con_rate > ceiling).sum() > 0:
                ceiling_loc.append(pd.Series(con_rate > ceiling).idxmax())
            if len(ceiling_loc) > 0:
                ceiling_loc = min(np.array(ceiling_loc))
                print( '\r' + 'rate ceiling reached at L=' + str(ceiling_loc), end = ' ')
                exp_rate[ceiling_loc:] = ceiling
                con_rate[ceiling_loc:] = ceiling
                nonexp_rate[ceiling_loc:] = ceiling
            name = name + '_ceiling_' + str(ceiling)
        if exp_zero == True:
            exp_rate *= 0
        if con_zero == True:
            con_rate *= 0
        if nonexp_zero == True:
            nonexp_rate *= 0
        if starting_counts == 'random':
            name = name + '_randomstart'        
        if starting_counts == 'subonly':
            name = name + '_subonlystart'
        if starting_counts == 'CHM13':
            name = name + '_CHM13start'
        if starting_counts == 'uniform':
            name = name + '_uniformstart'
        if stochastics is None:
            name = name
        else:
            name = name + '_stochastics_' + str(stochastics)
        return A_count_input, B_count_input, exp_rate, con_rate, nonexp_rate, B_indel_rate, name
    else:
        if boot is None:
            name = 'mutonly'
        else:
            name = 'mutonly_boot' + str(boot)
        if starting_counts == 'random':
            name = name + '_randomstart'        
        if starting_counts == 'subonly':
            name = name + '_subonlystart'
        if starting_counts == 'CHM13':
            name = name + '_CHM13start'
        if starting_counts == 'uniform':
            name = name + '_uniformstart'
        if stochastics is None:
            name = name
        else:
            name = name + '_stochastics_' + str(stochastics)
        return A_count_input, B_count_input, None, None, None, None, name


def run_simulation_constantspeedup(exp_int, con_int, exp_power, con_power, min_speedup = 1, rounds = 9, ceiling = 0.1, interp = False, nonexp_factor = 0.01, A_bins = 200, B_bins = 200, input_nuc = 'A', mutonly = False, stochastics = None, boot = None, boundary_count = 1000, overwrite = False, starting_counts = 'random', reflective = True, sim_dir = 'simulations/testing/', rates_function = 'powerlaw', first_bin = 0, write = True, startfrom = None, recnum = 11, pin_len = 9):
    starting_conditions = setup_evolve(exp_power, con_power, exp_int, con_int, nonexp_factor = nonexp_factor, A_bins = A_bins, B_bins = B_bins, input_nuc = input_nuc, mutonly = mutonly, starting_counts = starting_counts, boot = boot, ceiling = None, interp = interp, rates_function = rates_function, pin_len = pin_len)
    # set speedup to maximum allowable or reduce A_bins to avoid rate ceiling
    if mutonly != True:
        rates_sum = starting_conditions[2][:A_bins] + starting_conditions[3][:A_bins] + starting_conditions[4][:A_bins]
    else:
        rates_sum = denovo_substitution_context_rate.loc[input_nuc].sum() * 100
    speedup = -2 - round(np.log10(rates_sum.max()) - 0.5)
    if speedup < min_speedup:
        speedup = min_speedup; A_bins = np.array(list(np.where(rates_sum * 10**speedup >= ceiling)[0])).min() -1
        starting_conditions = setup_evolve(exp_power, con_power, exp_int, con_int, nonexp_factor = nonexp_factor, A_bins = A_bins, B_bins = B_bins, input_nuc = input_nuc, mutonly = mutonly, starting_counts = starting_counts, boot = boot, ceiling = None, interp = interp, rates_function = rates_function, pin_len = pin_len)
    rounds = rounds - speedup # default is 10^9 generations = 10^(iterations+speedup)

    if overwrite == False:
        if 'Adist_'+input_nuc+'_bins'+str(A_bins)+'_sp1e'+str(speedup)+'_rounds1e'+str(rounds)+'_'+starting_conditions[6]+'.pickle' in os.listdir(sim_dir):
            print('already done: ' + starting_conditions[6])
            return None
        else:
            pass
    else:
        pass

    if startfrom is None:
        A_counts_timeseries = dict(); B_counts_timeseries = dict()
        A_counts_timeseries[0] = starting_conditions[0][:A_bins]; B_counts_timeseries[0] = starting_conditions[1]
        A_counts_current = A_counts_timeseries[0]; B_counts_current = B_counts_timeseries[0]; flag = False; boundary_flag = False
        startrep = 0
        print('\r' + '         ' + starting_conditions[6], end = '     ')

    else:
        A_counts_timeseries, B_counts_timeseries = startfrom
        startrep = list(A_counts_timeseries.keys())[-1]
        A_counts_current = pd.Series(A_counts_timeseries[startrep]).reindex(range(A_bins)).fillna(0).values; B_counts_current = B_counts_timeseries[startrep]; flag = False; boundary_flag = False
        
    for rep in range(1, 1 + 10**rounds):
        if (flag == False):# & (max(A_counts_current[~np.isnan(A_counts_current)]) < 1e12):
            A_counts_current, B_counts_current, flag, boundary_flag = mut_evolve_dist_AB(A_counts_current, B_counts_current, starting_conditions[2:6], boot=boot, input_nuc = input_nuc, mutonly=mutonly, speedup_multiplier=10**speedup, stochastics = stochastics, reflective = reflective, boundary_count = boundary_count, first_bin = first_bin)
            # save 10 evenly-spaced (in log scale) timepoints per log10 rounds
            if rep in np.linspace(0,10**rounds,recnum).astype(int):
                print('\r' + str(rep), end = '   ')
                A_counts_timeseries[startrep + rep * (10**speedup)], B_counts_timeseries[startrep + rep * (10**speedup)] = A_counts_current, B_counts_current
        else:
            print('\r' + 'ending due to numerical error at round '+str(rep))
            break
    if write == True:
        A_counts_timeseries = pd.DataFrame(A_counts_timeseries)
        B_counts_timeseries = pd.DataFrame(B_counts_timeseries)
        A_counts_timeseries.to_pickle(sim_dir + 'Adist_'+input_nuc+'_bins'+str(A_bins)+'_sp1e'+str(speedup)+'_rounds1e'+str(rounds)+'_'+starting_conditions[6]+'.pickle')
        B_counts_timeseries.to_pickle(sim_dir + 'Bdist_'+input_nuc+'_bins'+str(A_bins)+'_sp1e'+str(speedup)+'_rounds1e'+str(rounds)+'_'+starting_conditions[6]+'.pickle')
    return A_counts_timeseries, B_counts_timeseries



def run_simulation_prospeedup(exp_int, con_int, exp_power, con_power, min_speedup = 4, interp = False, nonexp_factor = 0.01, A_bins = 200, B_bins = 200, input_nuc = 'A', mutonly = False, stochastics = None, boot = None, boundary_count = 1000, overwrite = False, starting_counts = 'subonly', reflective = True, sim_dir = 'simulations/testing/', rates_function = 'powerlaw', first_bin = 0, pin_len= 9):
    starting_conditions = setup_evolve(exp_power, con_power, exp_int, con_int, nonexp_factor = nonexp_factor, A_bins = A_bins, B_bins = B_bins, input_nuc = input_nuc, mutonly = mutonly, starting_counts = starting_counts, ceiling = None, interp = interp, boot = boot, rates_function = rates_function, pin_len = pin_len)
    if overwrite == False:
        if 'Adist_'+input_nuc+'_prospeedup_'+starting_conditions[6]+'.pickle' in os.listdir(sim_dir):
            print('already done: ' + starting_conditions[6])
            return None
        else:
            pass
    else:
        pass

    # set speedup to minimum, then increase to maximum allowable or reduce A_bins to avoid rate ceiling
    exp_rate = starting_conditions[2].copy(); con_rate = starting_conditions[3].copy()
    max_speedup=0
    exp_rate *= 10**max_speedup
    con_rate *= 10**max_speedup
    while exp_rate.max() + con_rate.max() < 0.01:
        max_speedup += 1
        exp_rate *= 10
        con_rate *= 10
    
    A_counts_timeseries, B_counts_timeseries = run_simulation_constantspeedup(exp_int=exp_int, con_int=con_int, exp_power=exp_power, con_power=con_power, min_speedup = min_speedup, rounds = 9, ceiling = 0.1, interp = interp, nonexp_factor = nonexp_factor, A_bins = A_bins, B_bins = B_bins, input_nuc = input_nuc, mutonly = mutonly, stochastics = stochastics, boot = boot, boundary_count = boundary_count, overwrite = None, starting_counts = starting_counts, reflective = reflective, sim_dir = None, rates_function = rates_function, first_bin = first_bin, write = False, recnum = 5, pin_len = pin_len)
    if max_speedup < min_speedup:
        for speedup in list(range(max_speedup, min_speedup))[::-1]:
            A_counts_timeseries, B_counts_timeseries = run_simulation_constantspeedup(exp_int=exp_int, con_int=con_int, exp_power=exp_power, con_power=con_power, min_speedup = speedup, rounds = 9 - min_speedup +speedup, ceiling = 0.1, interp = interp, nonexp_factor = nonexp_factor, A_bins = A_bins, B_bins = B_bins, input_nuc = input_nuc, mutonly = mutonly, stochastics = stochastics, boot = boot, boundary_count = boundary_count, overwrite = None, starting_counts = starting_counts, reflective = reflective, sim_dir = None, rates_function = rates_function, first_bin = first_bin, write = False, startfrom = (A_counts_timeseries, B_counts_timeseries), recnum = 5, pin_len = pin_len)
            
    A_counts_timeseries = pd.DataFrame.from_dict(A_counts_timeseries, orient = 'index').T
    B_counts_timeseries = pd.DataFrame.from_dict(B_counts_timeseries, orient = 'index').T
    A_counts_timeseries.to_pickle(sim_dir + 'Adist_'+input_nuc+'_prospeedup_'+starting_conditions[6]+'.pickle')
    B_counts_timeseries.to_pickle(sim_dir + 'Bdist_'+input_nuc+'_prospeedup_'+starting_conditions[6]+'.pickle')
    return A_counts_timeseries, B_counts_timeseries




import argparse
parser = argparse.ArgumentParser(description='repeat distribution simulation')

parser.add_argument('--dir', action="store", dest='dir', default = 'simulations/grid_output/', type=str)
parser.add_argument('--exp_p', action="store", dest='exp_p', type=float)
parser.add_argument('--con_p', action="store", dest='con_p', type=float)
parser.add_argument('--exp_i', action="store", dest='exp_i', type=int)
parser.add_argument('--con_i', action="store", dest='con_i', type=int)
parser.add_argument('--motif', action="store", dest='motif', default = 'A', type=str)
parser.add_argument('--speedup', action="store", dest='speedup', default = 5, type=int)
parser.add_argument('--rounds', action="store", dest='rounds', default = 0, type=int)
parser.add_argument('--A_bins', action="store", dest='A_bins', default = 200, type=int)
parser.add_argument('--B_bins', action="store", dest='B_bins', default = 200, type=int)
parser.add_argument('--boot', action="store", dest='boot', type = int)
parser.add_argument('--boundary_count', action="store", dest='boundary_count', default = 1000, type=int)
parser.add_argument('--mutonly', default=False, action="store_true")
parser.add_argument('--starting_counts', action="store", dest='starting_counts', default = 'random', type=str)
parser.add_argument('--overwrite', default=False, action="store_true")
parser.add_argument('--stochastics', action="store", dest='stochastics', type = int)
parser.add_argument('--reflective', default=True, action="store_false")
parser.add_argument('--ceiling', action="store", dest='ceiling', default=None, type=float)
parser.add_argument('--constantspeedup', default=False, action="store_true")
parser.add_argument('--interp', default=False, action="store_true")
parser.add_argument('--jobgroup', action="store", dest='jobgroup', type=int)
parser.add_argument('--rates_function', action="store", dest='rates_function', default = 'powerlaw', type=str)
parser.add_argument('--jobfile', action="store", dest='jobfile', default = 'grid_group_jobs.pickle', type=str)
parser.add_argument('--firstbin', action="store", dest='firstbin', default = 0, type=int)
parser.add_argument('--pin_len', action="store", dest='pin_len', default = 9, type=int)
parser.add_argument('--recnum', action="store", dest='recnum', default = 11, type=int)


args = parser.parse_args()
finished = os.listdir(args.dir)

if args.constantspeedup == False:
    if args.jobgroup is not None:
        if args.jobgroup == -1:
            grid_group = pd.read_pickle(args.jobfile)
        else:
            grid_group = pd.read_pickle(args.jobfile).loc[args.jobgroup]
        if __name__ == '__main__':
            pool = multiprocessing.Pool()
            for exp_i, con_i, exp_p, con_p in grid_group:
                pool.apply_async(run_simulation_prospeedup, args=(exp_i, con_i, exp_p, con_p, args.speedup, args.interp, 0.01, args.A_bins, args.B_bins, args.motif, args.mutonly, args.stochastics, args.boot, args.boundary_count, args.overwrite, args.starting_counts, args.reflective, args.dir, args.rates_function, args.firstbin, args.pin_len))
            pool.close()
            pool.join()
    else:
        run_simulation_prospeedup(exp_power=args.exp_p, con_power=args.con_p, exp_int=args.exp_i, con_int=args.con_i, boot = args.boot, input_nuc = args.motif, A_bins = args.A_bins, B_bins = args.B_bins, mutonly = args.mutonly, starting_counts = args.starting_counts, min_speedup = args.speedup, sim_dir = args.dir, overwrite = args.overwrite, stochastics = args.stochastics, reflective = args.reflective, rates_function = args.rates_function, first_bin = args.firstbin)
else:
    if args.jobgroup is not None:
        if args.jobgroup == -1:
            grid_group = pd.read_pickle(args.jobfile)
        else:
            grid_group = pd.read_pickle(args.jobfile).loc[args.jobgroup]
        if __name__ == '__main__':
            pool = multiprocessing.Pool()
            for exp_i, con_i, exp_p, con_p in grid_group:
                pool.apply_async(run_simulation_constantspeedup, args=(exp_i, con_i, exp_p, con_p, args.speedup, args.rounds, args.ceiling, args.interp, 0.01, args.A_bins, args.B_bins, args.motif, args.mutonly, args.stochastics, args.boot, args.boundary_count, args.overwrite, args.starting_counts, args.reflective, args.dir, args.rates_function, args.firstbin, True, None, args.recnum, args.pin_len))
            pool.close()
            pool.join()
    else:
        run_simulation_constantspeedup(exp_power=args.exp_p, con_power=args.con_p, exp_int=args.exp_i, con_int=args.con_i, boot = args.boot, input_nuc = args.motif, A_bins = args.A_bins, B_bins = args.B_bins, mutonly = args.mutonly, starting_counts = args.starting_counts, min_speedup = args.speedup, rounds = args.rounds, ceiling = args.ceiling, sim_dir = args.dir, overwrite = args.overwrite, stochastics = args.stochastics, reflective = args.reflective, rates_function = args.rates_function, first_bin = args.firstbin)

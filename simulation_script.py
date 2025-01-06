#!/usr/bin/env python
import pandas as pd
import numpy as np
import pickle
import os


# load data
CHM13_counts = pd.read_pickle('repeat_distributions/CHM13_counts.pickle')
random_counts = pd.read_pickle('repeat_distributions/random_counts.pickle')
subonly_counts = pd.read_pickle('repeat_distributions/subonly_counts.pickle')

denovo_exp_rate = pd.read_pickle('denovo/denovo_exp_rate.pickle')
denovo_con_rate = pd.read_pickle('denovo/denovo_con_rate.pickle')
denovo_nonexp_rate = pd.read_pickle('denovo/denovo_nonexp_rate.pickle')

#denovo_mut_freq_AB = pd.read_pickle('denovo/denovo_mut_freq_AB.pickle')
#denovo_mut_freq_AB_poisson = pd.read_pickle('denovo/denovo_mut_freq_AB_poisson.pickle')

denovo_exp_rate_poisson = pd.read_pickle('denovo/denovo_exp_rate_poisson.pickle')
denovo_con_rate_poisson = pd.read_pickle('denovo/denovo_con_rate_poisson.pickle')
denovo_nonexp_rate_poisson = pd.read_pickle('denovo/denovo_nonexp_rate_poisson.pickle')

decode_exp_rate_poisson = pd.read_pickle('decode/decode_expansion_rates_poisson.pickle')
decode_con_rate_poisson = pd.read_pickle('decode/decode_contraction_rates_poisson.pickle')

denovo_substitution_context_rate = pd.read_pickle('denovo/denovo_mut_freq_triplets.pickle')
denovo_substitution_context_rate_poisson = pd.read_pickle('denovo/denovo_mut_freq_triplets_poisson.pickle')


# evolve function
def mut_evolve_dist_AB(A_count_input, B_count_input, starting_conditions, boot = None, input_nuc = 'A', mut = True, mutonly = False, speedup_multiplier = 1, output_components = False, stochastics = None, reflective = True):
    exp_rate_A_AA, con_rate_A_AA, nonexp_rate_A_AB, B_indel_rates = starting_conditions
    A_count_output = A_count_input.copy(); B_count_output = B_count_input.copy()
    A_bins = len(A_count_input)
    B_bins = len(B_count_input)
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
    total_B_bases = (B_count_input[:B_bins] * B_length_array_bases[:B_bins]).sum()
    B_L1_base_portion = ((B_count_input[0] * len(input_nuc)) / (B_count_input[:B_bins]* B_length_array_bases[:B_bins]).sum()) ### including motif length
    B_nonflank_base_portion = (B_count_input[2:B_bins+2] * B_length_array_bases[:B_bins]).sum() / total_B_bases  ### include portion of triplets 1nt away???
    B_flank_base_portion = (B_count_input[1:B_bins] * 2 * len(input_nuc)).sum() / total_B_bases ### including motif length
    
    total_A_bases = (A_count_input[:A_bins] * A_length_array_bases[:A_bins]).sum()
    A_nonflank_base_portion = (A_count_input[2:A_bins+2] * A_length_array_bases[:A_bins]).sum() / total_A_bases
    A_flank_base_portion = (A_count_input[1:A_bins] * 2 * len(input_nuc)).sum() / total_A_bases ### including motif length
    
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
        A_mut_in_fission =  np.array([np.sum(((2/A_length_array[:A_bins]) * -A_mut_out_fission[2:A_bins+2])[L-1:]) for L in A_length_array[:A_bins]]) ### use length_array_bases???
 
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
        B_mut_in_fission =  np.array([np.sum(((2/B_length_array[:B_bins]) * -B_mut_out_fission[2:B_bins+2])[L-1:]) for L in B_length_array[:B_bins]]) ### use length_array_bases???

        # A>B which adds to the B count locally (which must come from A_L>1)
        # B from A>B leaving the -1 bin
        B_len_freq = (B_count_input / B_count_input.sum())[:B_bins]
        B_mut_out_local_A_B = -denovo_sub['Acontraction'] * A_flank_base_portion * total_A_bases * B_len_freq
        # A>B creating B_L=1 from A_L>2
        A_B_into_L1 = total_A_bases * A_nonflank_base_portion * denovo_sub['Afission']
#        A_B_into_L1 = -A_mut_out_fission.sum()
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
        A_nonexp_fissions_out = -A_count_input * nonexp_rate_A_AB # used to calculate fission_in, starting with L=2 going to 2x L=1
        # each fission creates 2 As. add these to A
        A_nonexp_in_fission =  np.array([np.sum(((2/A_length_array[:A_bins]) * -A_nonexp_fissions_out[1:A_bins+1])[L-1:]) for L in A_length_array[:A_bins]])


        # B expansions in and out
        B_exp_out = B_count_input[:B_bins] * -B_indel_rates[2] * B_length_array[:B_bins] # B>BB rates are flat, per base
        B_exp_in = np.insert(-B_exp_out, 0, A_nonexp_fissions_out.sum())

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
        B_nonexp_fissions_out = -B_count_input * B_indel_rates[0] * B_length_array # used to calculate fission_in, starting with L=2 going to 2x L=1
        # each fission creates 2 Bs. add these to B
        B_nonexp_in_fission =  np.array([np.sum(((2/B_length_array[:B_bins]) * -B_nonexp_fissions_out[1:B_bins+1])[L-1:]) for L in B_length_array[:B_bins]])
            
       # update counts for next round (with absorbing boundary)
        total_A_change_in += A_exp_in[:A_bins] + A_con_in[:A_bins] + A_mut_in_fusion_Bdel[:A_bins] + A_nonexp_in_fission[:A_bins]
        total_B_change_in += B_exp_in[:B_bins] + B_con_in[:B_bins] + B_mut_in_fusion_Adel[:B_bins] + B_nonexp_in_fission[:B_bins]
        total_A_change_out += A_exp_out[:A_bins] + A_con_out[:A_bins] + A_mut_out_fusion_Bdel[:A_bins] + A_nonexp_fissions_out[:A_bins]
        total_B_change_out += B_exp_out[:B_bins] + B_con_out[:B_bins] + B_mut_out_fusion_Adel[:B_bins] + B_nonexp_fissions_out[:B_bins]

        # apply reflecting boundary
        if reflective == True:
            total_A_change_in[A_bins-1] += A_exp_in[A_bins:].sum() + A_con_in[A_bins:].sum() + A_mut_in_fusion_Bdel[A_bins:].sum() + A_nonexp_in_fission[A_bins:].sum()
            total_B_change_in[B_bins-1] += B_exp_in[B_bins:].sum() + B_con_in[B_bins:].sum() + B_mut_in_fusion_Adel[B_bins:].sum() + B_nonexp_in_fission[B_bins:].sum()
            total_A_change_out[A_bins-1] += A_exp_out[A_bins:].sum() + A_con_out[A_bins:].sum() + A_mut_out_fusion_Bdel[A_bins:].sum() + A_nonexp_fissions_out[A_bins:].sum()
            total_B_change_out[B_bins-1] += B_exp_out[B_bins:].sum() + B_con_out[B_bins:].sum() + B_mut_out_fusion_Adel[B_bins:].sum() + B_nonexp_fissions_out[B_bins:].sum()

    
    # flag to stop the simulation if more repeats are removed from a bin than exist in that bin (excluding the last 10 noisy bins)
    flag = ((np.abs(total_A_change_out[:A_bins-10]) * speedup_multiplier > A_count_output[:A_bins-10]).sum()) > 0

    # apply speedup
    total_A_change_in *= speedup_multiplier; total_A_change_out *= speedup_multiplier
    total_B_change_in *= speedup_multiplier; total_B_change_out *= speedup_multiplier
    
    if stochastics is not None:
        # the sum of poisson random variables is poisson-distributed. not necessary to run n poisson samples
        total_A_change_in = np.random.poisson(total_A_change_in)
        total_A_change_out = -1 * np.random.poisson(np.abs(total_A_change_out))
        total_B_change_in = np.random.poisson(total_B_change_in)
        total_B_change_out = -1 * np.random.poisson(np.abs(total_B_change_out))   
    
    total_A_change = total_A_change_in + total_A_change_out
    total_B_change = total_B_change_in + total_B_change_out
    
    # update counts for next round
    A_count_output = A_count_output[:A_bins] + total_A_change[:A_bins]
    B_count_output = B_count_output[:B_bins] + total_B_change[:B_bins]

    # remove negative values
    A_count_output[A_count_output <0] = 0            
    B_count_output[B_count_output <0] = 0
    
    if output_components == True:
        if mutonly == False:
            return  A_mut_in_local_A_B[:A_bins], A_mut_out_local_A_B[:A_bins], A_mut_in_local_B_A[:A_bins], A_mut_out_local_B_A[:A_bins], A_mut_in_fission[:A_bins], A_mut_out_fission[:A_bins], A_mut_in_fusion_A_B[:A_bins], A_mut_out_fusion_A_B[:A_bins], A_exp_in[:A_bins], A_exp_out[:A_bins], A_con_in[:A_bins], A_con_out[:A_bins], A_mut_in_fusion_Bdel[:A_bins], A_mut_out_fusion_Bdel[:A_bins], A_nonexp_in_fission[:A_bins], A_nonexp_fissions_out[:A_bins]
        else:
            return  A_mut_in_local_A_B[:A_bins], A_mut_out_local_A_B[:A_bins], A_mut_in_local_B_A[:A_bins], A_mut_out_local_B_A[:A_bins], A_mut_in_fission[:A_bins], A_mut_out_fission[:A_bins], A_mut_in_fusion_A_B[:A_bins], A_mut_out_fusion_A_B[:A_bins]
    else:
        return A_count_output, B_count_output, flag


def extend_power_law(power, start_rate, start_len, end_len=100):
    denom = (start_len**power) / start_rate
    return pd.Series([i**power for i in range(start_len+1, end_len+1)], index = list(range(start_len+1,end_len+1))) / denom

def multiply_then_powerlaw(exp_power, con_power, mult, A_bins = 100, boot = None, L_mult = 9, L_mult_nonexp = 9, motif = 'A', fill = False, nonexp_factor = False):
    if boot is None:
        bootname = ''
        denovo_exp_rate_current = denovo_exp_rate[motif]
        denovo_con_rate_current = denovo_con_rate[motif]
        denovo_nonexp_rate_current = denovo_nonexp_rate[motif]
    else:
        bootname = '_boot'+str(boot)
        denovo_exp_rate_current = denovo_exp_rate_poisson[motif][boot]
        denovo_con_rate_current = denovo_con_rate_poisson[motif][boot]
        denovo_nonexp_rate_current = denovo_nonexp_rate_poisson[motif][boot]
    if fill == True:
        fillname = '_fill'
        denovo_exp_rate_current = denovo_exp_rate_current.replace(0, np.nan).interpolate(method = 'from_derivatives')
        denovo_con_rate_current = denovo_con_rate_current.replace(0, np.nan).interpolate(method = 'from_derivatives')
        denovo_nonexp_rate_current = denovo_nonexp_rate_current.replace(0, np.nan).interpolate(method = 'from_derivatives')
    else:
        fillname = ''
    exp = pd.concat([denovo_exp_rate_current.reindex(range(L_mult)), pd.Series(denovo_exp_rate_current[L_mult-1] * mult, index = [L_mult]), extend_power_law(exp_power, denovo_exp_rate_current[L_mult-1]*mult, L_mult, A_bins+3)])
    con = pd.concat([denovo_con_rate_current.reindex(range(L_mult)), pd.Series(denovo_con_rate_current[L_mult-1] * mult, index = [L_mult]), extend_power_law(con_power, denovo_con_rate_current[L_mult-1]*mult, L_mult, A_bins+3)])
    if nonexp_factor == False:
        nonexpname = ''
        nonexp = pd.concat([denovo_nonexp_rate_current.reindex(range(L_mult_nonexp)), pd.Series(denovo_nonexp_rate_current[L_mult_nonexp-1] * mult, index = [L_mult_nonexp]), extend_power_law(exp_power, denovo_nonexp_rate_current[L_mult_nonexp-1]*mult, L_mult_nonexp, A_bins+3)])
    else:
        nonexpname = '_nonexp_x' + str(nonexp_factor)
        nonexp = exp * nonexp_factor
    nonexp.loc[1] = 0
    name = 'mult_L'+str(L_mult)+'_x'+str(mult)+'_extend_pl_' + str(exp_power) + '_' + str(con_power) + nonexpname + fillname + bootname
    return name, exp, con, nonexp


# setup function
def setup_evolve(exp_power=1, con_power=1, mult=1, boot = None, stochastics = None, L_mult = 9, L_mult_nonexp = 9, fill = False, nonexp_factor = False, A_bins = 100, B_bins = 100, input_nuc = 'A', mutonly = False, exp_zero = False, con_zero = False, nonexp_zero = False, different_input = False, random_start = False, subonly_start = False, ceiling = None):
# set up counts
    A_length_array = np.array(range(1,A_bins+1))
    B_length_array = np.array(range(1,A_bins+1))
    if different_input == False:
        if random_start == False:
            if subonly_start == False:
                A_count_input = np.nan_to_num(CHM13_counts['A'][input_nuc].reindex(range(1,A_bins+1)).values)
                B_count_input = np.nan_to_num(CHM13_counts['B'][input_nuc].reindex(range(1,B_bins+1)).values)
            if subonly_start == True:
                A_count_input = np.nan_to_num(subonly_counts['A'][input_nuc].reindex(range(1,A_bins+1)).values)
                B_count_input = np.nan_to_num(subonly_counts['B'][input_nuc].reindex(range(1,B_bins+1)).values)
        if random_start == True:
            A_count_input = np.nan_to_num(random_counts['A'][input_nuc].reindex(range(1,A_bins+1)).values)
            B_count_input = np.nan_to_num(random_counts['B'][input_nuc].reindex(range(1,B_bins+1)).values)
    else:
        A_count_input = np.nan_to_num(different_input[0].reindex(range(1,A_bins+1)).values)
        B_count_input = np.nan_to_num(different_input[1].reindex(range(1,B_bins+1)).values)
# set up rates    
    if mutonly == False:
        name, exp_rate, con_rate, nonexp_rate = multiply_then_powerlaw(exp_power = exp_power, con_power = con_power, A_bins = A_bins, mult=mult, boot = boot, L_mult = 9, L_mult_nonexp = 9, motif = input_nuc, fill = fill, nonexp_factor = nonexp_factor)
        B_indel_rate = np.array([exp_rate[0], con_rate[0], nonexp_rate[0]])        
        # change rates from per unit to per STR
        exp_rate = exp_rate.values[1:A_bins+1] * A_length_array
        con_rate = con_rate.values[1:A_bins+2] * np.array(range(1,A_bins+2))
        nonexp_rate = nonexp_rate.values[1:A_bins+3] * np.array(range(1,A_bins+3))
        if ceiling != None:
            exp_rate[exp_rate > ceiling] = ceiling
            con_rate[con_rate > ceiling] = ceiling
            nonexp_rate[nonexp_rate > ceiling] = ceiling
            name = name + '_ceiling_' + str(ceiling)
        if exp_zero == True:
            exp_rate *= 0
        if con_zero == True:
            con_rate *= 0
        if nonexp_zero == True:
            nonexp_rate *= 0
        if random_start == False:
            name = name + '_CHM13start'
        if random_start == True:
            name = name + '_randomstart'
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
        if random_start == False:
            name = name + '_CHM13start'
        if random_start == True:
            name = name + '_randomstart'
        if stochastics is None:
            name = name
        else:
            name = name + '_stochastics_' + str(stochastics)
        return A_count_input, B_count_input, None, None, None, None, name


def run_simulation(exp_power=1, con_power=1, mult=1, boot = None, L_mult = 9, L_mult_nonexp = 9, fill = False, nonexp_factor = False, A_bins = 100, B_bins = 100, input_nuc = 'A', mutonly = False, speedup = 3, rounds = 5, overwrite = False, stochastics = None, random_start = False, subonly_start = False, ceiling = False, reflective = True, sim_dir = 'grid_output/'):
    starting_conditions = setup_evolve(exp_power=exp_power, con_power=con_power, mult=mult, boot = boot, stochastics = stochastics, L_mult = L_mult, L_mult_nonexp = L_mult_nonexp, fill = fill, nonexp_factor = nonexp_factor, A_bins = A_bins, B_bins = B_bins, input_nuc = input_nuc, mutonly = mutonly, random_start = random_start, subonly_start = subonly_start, ceiling = ceiling)

    if overwrite == False:
        if 'Adist_bins'+str(A_bins)+'_sp1e'+str(speedup)+'_rounds1e'+str(rounds)+'_'+starting_conditions[6]+'.pickle' in finished:
            print('already done: ' + starting_conditions[6])
            return None
        else:
            pass
    else:
        pass        
    print('\r' + '         ' + starting_conditions[6], end = '     ')
    A_counts_timeseries = dict(); B_counts_timeseries = dict()
    A_counts_timeseries[0] = starting_conditions[0]; B_counts_timeseries[0] = starting_conditions[1]
    A_counts_current = A_counts_timeseries[0]; B_counts_current = B_counts_timeseries[0]; flag = False
    for rep in range(1, 1 + 10**rounds):
        if (flag == False):# & (max(A_counts_current[~np.isnan(A_counts_current)]) < 1e12):
            A_counts_current, B_counts_current, flag = mut_evolve_dist_AB(A_counts_current, B_counts_current, starting_conditions[2:6], boot=boot, input_nuc = input_nuc, mutonly=mutonly, speedup_multiplier=10**speedup, stochastics = stochastics, reflective = reflective)
            if rep%int(max(1, 1e6/(10**speedup))) == 0:
                print('\r' + str(rep), end = '   ')
                A_counts_timeseries[rep], B_counts_timeseries[rep] = A_counts_current, B_counts_current
        else:
            print('\r' + 'ending due to numerical error at round '+str(rep))
            break
    A_counts_timeseries = pd.DataFrame(A_counts_timeseries)
    B_counts_timeseries = pd.DataFrame(B_counts_timeseries)
    A_counts_timeseries.to_pickle(sim_dir + 'Adist_'+input_nuc+'_bins'+str(A_bins)+'_sp1e'+str(speedup)+'_rounds1e'+str(rounds)+'_'+starting_conditions[6]+'.pickle')
    B_counts_timeseries.to_pickle(sim_dir + 'Bdist_'+input_nuc+'_bins'+str(A_bins)+'_sp1e'+str(speedup)+'_rounds1e'+str(rounds)+'_'+starting_conditions[6]+'.pickle')
    return A_counts_timeseries, B_counts_timeseries


import argparse
parser = argparse.ArgumentParser(description='repeat distribution simulation')

parser.add_argument('--dir', action="store", dest='dir', default = 'simulations/grid_output/', type=str)
parser.add_argument('--mult', action="store", dest='mult', type=float)
parser.add_argument('--exp', action="store", dest='exp', type=float)
parser.add_argument('--con', action="store", dest='con', type=float)
parser.add_argument('--L_mult', action="store", dest='L_mult', default = 9, type=int)
parser.add_argument('--L_mult_nonexp', action="store", dest='L_mult_nonexp', default = 9, type=int)
parser.add_argument('--motif', action="store", dest='motif', default = 'A', type=str)
parser.add_argument('--speedup', action="store", dest='speedup', default = 3, type=int)
parser.add_argument('--rounds', action="store", dest='rounds', default = 5, type=int)
parser.add_argument('--A_bins', action="store", dest='A_bins', default = 100, type=int)
parser.add_argument('--B_bins', action="store", dest='B_bins', default = 200, type=int)
parser.add_argument('--boot', action="store", dest='boot', type = int)
parser.add_argument('--mutonly', default=False, action="store_true")
parser.add_argument('--random_start', default=False, action="store_true")
parser.add_argument('--subonly_start', default=False, action="store_true")
parser.add_argument('--overwrite', default=False, action="store_true")
parser.add_argument('--stochastics', action="store", dest='stochastics', type = int)
parser.add_argument('--ceiling', action="store", dest='ceiling', default=None, type=float)
parser.add_argument('--reflective', default=True, action="store_false")

args = parser.parse_args()


finished = os.listdir(args.dir)


run_simulation(mult=args.mult, exp_power=args.exp, con_power=args.con, boot = args.boot, L_mult = args.L_mult, L_mult_nonexp = args.L_mult_nonexp, input_nuc = args.motif, A_bins = args.A_bins, B_bins = args.B_bins, mutonly = args.mutonly, random_start = args.random_start, subonly_start = args.subonly_start, speedup = args.speedup, rounds = args.rounds, sim_dir = args.dir, overwrite = args.overwrite, stochastics = args.stochastics, ceiling = args.ceiling, reflective = args.reflective)

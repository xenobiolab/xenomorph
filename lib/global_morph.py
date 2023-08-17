########################################################################
########################################################################
"""
global_morph.py

Description: This script is used for performing global alternative 
hypthesis testing by method of MLE. Global level values are calculated by
taking the average kmer level across all observations of the same sequence. 
Unlike morph, global_morph does not calculate per-read statistics. 

To use global instead of per-read, add the -g flag to morph: 
xenomorph.py morph -g 

Title: Synthesis and Sequencing of 12-Letter Supernumerary DNA

By: H. Kawabe, C. Thomas, S. Hoshika, Myong-Jung Kim, Myong-Sang Kim, L. Miessner, J. M. Craig, 
J. Gundlach, A. Laszlo,  S. A. Benner, J. A. Marchand

Updated: 2/14/23
"""
########################################################################
########################################################################


import pandas as pd
import numpy as np
import itertools
import seaborn as sns
import os
import sys
from itertools import compress
from string import ascii_lowercase
import scipy.stats as stats
import matplotlib.pyplot as plt
import statistics
import warnings
warnings.filterwarnings("ignore")
import multiprocessing
import random
import math
from alive_progress import alive_bar
from xm_params import * 
from xm_tools import *


def kmer2level(kmers, kmer_model,n_iter):
    level=[]

    for i in range(0,len(kmers)):
        kmer_select = kmers[i]
                
        #Fetch mean from kmer model
        kmer_s = kmer_model[kmer_model['KXmer']==kmer_select]
        mu = float(kmer_s['mu'].iloc[0]) #Fetch mean
        
        #Fetch standard deivation from kmer model
        noise = float(kmer_s['sigma'].iloc[0])
        
        #Generate error
        st = np.random.normal(0,noise, n_iter)
        
        #Simulate an observation
        mus = mu + st
        
        #Store observed level in array
        level.append(mus)
        
    #Return simulated levels
    return(level)

#Converts a sequence seq into kmers of length k
def seq2kmer(seq, k):
    kmers = []
    for i in range(0,len(seq)-k+1):
        kmers.append(seq[i:i+k].upper())
    return kmers
    
def kmer2seq(kmers):
    seq=kmers[0]
    for i in range(1, len(kmers)):
        seq=seq+kmers[i][-1]
    return(seq)


#Samples levels from a raw_kmer_file. Can be used to simulate data instead of kmer2level
#Return row_wize set of level simulations
def sample2levels(kmers, raw_kmers, n_iter):


    level_sets =[] 
    for i in range(0,len(kmers)):
        kmer_select = kmers[i]
                
        #Fetch mean from kmer model
        kmer_s = raw_kmers[raw_kmers['kmer_xy']==kmer_select]
        lev= kmer_s.iloc[0]['mean_level'].replace('\n',' ').replace("'",'').replace(']','').replace('[','').split(' ')
        


        level=[]
        for i in range(0,len(lev)): 
            try:
                level.append(float(lev[i]))
            except: 
                ex = 0 
        level_sets.append(random.choices(level, k = n_iter))

    return level_sets


def readlevel2kmerlevel(read_level, kmer_mask): 
    z_pos = kmer_mask.find('x')
    read_level = read_level.replace('\n',' ').replace("'",'').replace(']','').replace('[','').split(' ')
    read_level = [i for i in read_level if i]
    read_level_kmers = read_level[z_pos:len(read_level)-len(kmer_mask[z_pos:-1])]
    return read_level_kmers

#Input a Kmer string, normalized level, and a kmer_distribution CSV file. Outputs probability of observing event (P e|x)
def kmerpdf(kmer_select,event,kmer_model):
    #Generate pdf of kmer levels with mean mu and standard deviation std
    kmer_s = kmer_model[kmer_model['KXmer']==kmer_select]
    mu =kmer_s['mu'].iloc[0]
    st = kmer_s['sigma'].iloc[0]
    return stats.norm(float(mu),float(st)).pdf(event)

 
 
#Input a path of kmers to take, and output the log probability of that path
def chart_path(kmer_path, levels, kmer_model, kmer_weights):
    p = 0

    for i in range(0,len(kmer_path)):
        kmer_select = kmer_path[i]
        event = float(levels[i])

        pj = (np.log10(kmerpdf(kmer_select,event,kmer_model)))*kmer_weights[i]
        p = p + pj

    return (p)


def chart_path_robust(kmer_path_can, kmer_path_alt, levels, kmer_model,kmer_weights ):


    orllr = 0
    for i in range(0,len(levels)):

        #Convert level to float, stored as string up until this point
        level = float(levels[i])   

        #Kmer from canonical distribution (this is the variable sequence)
        kmer_can = kmer_path_can[i]

        #Kmer from alternative distribution (this is the fixed, XNA-containing sequence)
        kmer_alt = kmer_path_alt[i] 

        #Likelihood that observed level belongs to canonical_kmer
        LL_can = np.log10(kmerpdf(kmer_can,level,kmer_model))

        #Likelihood that observed level belongs to alt_kmer
        LL_alt = np.log10(kmerpdf(kmer_alt,level,kmer_model))

        #Get LLR - Negative means the null hypothesis is more likely
        LLR = (LL_alt - LL_can)*kmer_weights[i]

        #Canonical mean
        cmean = kmer_model[kmer_model['KXmer']==kmer_path_can[i]]['mu'].values[0]

        #Alternative base mean 
        amean = kmer_model[kmer_model['KXmer']==kmer_path_alt[i]]['mu'].values[0]

        #Canonical standard deviation (should be fixed at set sigma)
        csig = (kmer_model[kmer_model['KXmer']==kmer_path_can[i]]['sigma'].values[0])

        #Alternative standard deviation (should be fixed at set sigma)
        asig = (kmer_model[kmer_model['KXmer']==kmer_path_alt[i]]['sigma'].values[0])

        #Take average of standard deviation, then square to get an average variance
        acvar = ((csig+asig)/2)**2 

        #Scale difference between means and observed level
        scaleDiff = float(level) - ((cmean+amean)/2)

        #Difference in means 
        meanDiff = np.abs(cmean - amean) 

        #Corrected outlier robust from Tombo documentation - Negative exponential included
        outrobLLR = (np.exp(-(scaleDiff **2) / (Sf * (acvar))) * LLR) / ((acvar) * (meanDiff ** Sp) * Sf2)

        #Sum outlier-robust log likelihoods over the Kmers
        orllr = orllr + outrobLLR

    #For alt-base = canonical-base, set to 100 (here, positive orllr = not likely)
    if math.isnan(orllr)==True: 
        orllr =50
    return orllr
    


def alt_type(bases, xna_base, alt_base_type): 
    if alt_base_type == 'all': 
        set_alt_base = bases
    elif alt_base_type == 'standard':
        set_alt_base = standard_bases
    elif alt_base_type =='pyrpur': 
    

        if xna_base_rc(xna_base,confounding_pairs) in 'AG': 
            set_alt_base = ['A','G']
        elif xna_base_rc(xna_base,confounding_pairs) in 'TC': 
            set_alt_base = ['T','C']
        else: 
            print("Xenomorph Status - [Error] XNA base not classified as purine or pyrimidine (setting alt_base = all). Check lib/xm_params.py confounding pairs.")
            set_alt_base = bases
    elif alt_base_type =='confounding': 
        set_alt_base = xna_base_rc(xna_base,confounding_pairs)
    else: 
        print("Xenomorph Status - [Error] Incorrect alt_base_type supplied (setting alt_base = all). Check lib/xm_params.py.")
        set_alt_base = bases


    #Add XNA base to the list
    set_alt_base = list(set_alt_base)
    set_alt_base.append(xna_base)
    set_alt_base = list(set(set_alt_base))

    return set_alt_base

    
#Strickly uses 7th base position. Given a sequence and list of levels, will calculate all alternative likelyhoods 
def gen_alt_all(sequence, kmer_levels, kmer_model, all_bases, xbase_pos, kmer_len):
    
    #Store np array of log likelihoods 
    lpp=[]
    lrr=[]

    #Setup alternative hypothesis testing bases based on specifications from alt_base_type parameter. 
    xna_base = sequence[xbase_pos]


    all_bases = alt_type(all_bases, xna_base, alt_base_type)



    for j in range(0,len(all_bases)):

        #Generate alternative sequence - Uses fixed 7mer index 
        seq_alt = sequence[0:xbase_pos]+all_bases[j]+sequence[xbase_pos+1:]


        #Convert sequence to kmer list, then calculate liklihoods 
        lp = chart_path(seq2kmer(seq_alt,kmer_len), kmer_levels, kmer_model, kmer_weights)

        #store likelihood output in lpp array 
        lpp.append(lp)


        #Calculate log likelihood robust ratio 
        lr = chart_path_robust(seq2kmer(seq_alt,kmer_len), seq2kmer(sequence,kmer_len), kmer_levels, kmer_model, kmer_weights)
        
        #Set XNA base some higher value since it should be giving nan LLR
        if all_bases[j]==xna_base:
            lr = 100

        #store likelihood output in lpp array 
        lrr.append(lr)


    x_sub = sequence[xbase_pos]


    #Determine most likley bsae 
    if likelihood_ratio == 'Normal': 
        #Calculate log-likelihood ratio that the predicted basecall is right 
        ppl_bc = list(compress(all_bases, lpp == np.amax(lpp, axis=0)))[0]
        x_sub = sequence[xbase_pos]
        x_lpp=lpp[all_bases.index(x_sub)]


        lpp.remove(x_lpp)
        a_lpp=np.max(lpp)
        llhr = x_lpp - a_lpp
        most_likely_base = ppl_bc 

    
    if likelihood_ratio == 'Outlier-Robust': 
        #If all values of LLR are >0 new basecall is decided
        if all(x >0 for x in lrr):
            lrr_bc = all_bases[lrr.index(100)]
        elif np.min(lrr)>= likelihood_threshold: 
            lrr_bc = all_bases[lrr.index(100)]
        else:
            lrr_bc = all_bases[lrr.index(np.min(lrr))]
        most_likely_base = lrr_bc
        llhr = np.min(lrr)

    return most_likely_base, llhr
    


        
 
############################
#Load kmer measurements
############################



##Load Kmer model distribution file
kmer_model_input_file = sys.argv[1]
kmer_model_input = pd.read_csv(kmer_model_input_file, sep=',')


##Load read levels raw 
read_level_summary = sys.argv[2]
read_level_summary= pd.read_csv(read_level_summary, sep=',')

##Filter out reads that do not pass quality settings in xm_params
print('Xenomorph Status - [Stats] Filtering reads with q score > '+str(qscore_filter))
print('Xenomorph Status - [Stats] Filtering reads with signal score < '+str(signal_filter))
read_level_summary=read_level_summary[read_level_summary['read_q-score']>qscore_filter]
read_level_summary=read_level_summary[read_level_summary['read_signal_match_score']<signal_filter]


#Output file name to save 
out_fn = sys.argv[3]


#Load all possible bases from model file 
all_bases = list(set(''.join(kmer_model_input['KXmer'].to_list())))



###################################

##############################################

#Model configuration

sig = kmer_model_input['Std level']

if isinstance(sigma, str):
    if sigma == 'Global-Median':
        sig = np.median(sig)
    if sigma == 'Global-Mean':
        sig = np.mean(sig)
else:
    sig = sigma
    
kmer_model_input['Std level'] = sig

kmer_model=kmer_model_input[['KXmer',mu,'Std level']]

kmer_model.rename(columns = {mu:'mu', 'Std level':'sigma'}, inplace = True)

kmer_list = kmer_model['KXmer']

#Check which bases are loaded into model set
model_bases_loaded = list(set(''.join(kmer_list)))

#Filter out level file sequences that contain bases model cannot decode 
all_xna = list(set(read_level_summary['read_xna']))
decodeable_bases= list(set(model_bases_loaded) & set(all_xna))
read_level_summary=read_level_summary[read_level_summary['read_xna'].isin(decodeable_bases)]


##############################################
#Generate kmer data sampling from real dataset
##############################################


#Generate every possible nnn-mer
xbase_pos = nmer_model.find('x') 

#Note - Generate this warning specifying method of kmer mask being used
print('Xenomorph Status - [Global Morph] Using default kmer mask nxnn')


ref_seqs = list(set(read_level_summary['reference_sequence']))


#Setup output dataframe
output_column_names=['reference_sequence','n_reads','mean_q-score','mean_signal_match_score','reference_locus','strand','xna_sequence','xna', 'basecall', 'llhr','xna_is_concensus', 'model_sigma', 'model_mean', 'model_file']
output_summary = pd.DataFrame(columns = output_column_names)

#numbers to calculate consensus count
consensus_count=0
total_count= {}
global_count = {}


ref_seqs.sort()
with alive_bar(len(ref_seqs), force_tty=True) as bar: 
    for i in range(0,len(ref_seqs)): 

        #Set of reads matching reference
        ref_i = read_level_summary[read_level_summary['reference_sequence']==ref_seqs[i]]

        #Reference sequence 
        ref_i_id = ref_seqs[i] 

        #Strands sequenced (+,-)
        ref_s = list(set(ref_i['read_xna_strand']))

        for j in range(0,len(ref_s)): 

            #Strand
            ref_strand=ref_s[j]

            #Reads for strand 
            ref_i_s = ref_i[ref_i['read_xna_strand']==ref_strand]

            #Number of reads averaged for measurement
            ref_i_n = len(ref_i_s)

            #Read reference sequence (= to locus if (+) strand, reverse complement of locus for (-) strand)
            sequence = ref_i_s.iloc[0]['read_xna_sequence']


            #Check that model can decode sequence 
            if len(list(set(sequence)-set(model_bases_loaded)))==0:

                #Extract read-kmer levels and take mean 
                ref_i_kmer_levels =[]
                for k in range(0,len(ref_i_s)):
                    ref_i_levels = ref_i_s.iloc[k]['read_levels']
                    ref_i_kmer_level_s = readlevel2kmerlevel(ref_i_levels,kmer_mask[0])
                    ref_i_kmer_levels.append([float(i) for i in ref_i_kmer_level_s])


                #Levels (averaged) 
                if 'mean' in (mu_global).lower(): 
                    levels = np.mean(ref_i_kmer_levels, axis = 0)
                if 'median' in (mu_global).lower(): 
                    levels = np.median(ref_i_kmer_levels, axis = 0)
        
                #Corresponding kmers 
                kmer = seq2kmer(sequence,len(kmer_mask[0]))

                #Calculate log liklihoods
                lpp = gen_alt_all(sequence,levels, kmer_model, all_bases,xbase_pos, len(kmer_mask[0]))

                #Basecall made
                ref_i_bc = lpp[0]

                #Basecall made
                ref_i_llhr = lpp[1]

                #XNA prior
                ref_i_x = sequence[xbase_pos]

                #Read reference locus
                ref_i_loc = ref_i_s.iloc[0]['read_reference_locus']

                #Sequence 
                ref_i_seq = sequence

                #Check if XNA is one of the conensus base calls
                ref_i_x_is_concensus = ref_i_bc == sequence[xbase_pos]

                #Average q-score
                ref_i_q = np.mean(ref_i_s['read_q-score'].to_list())

                #Average signal score
                ref_i_sc = np.mean(ref_i_s['read_signal_match_score'].to_list())

                #model sigma used 
                ref_i_sig = sig

                #model mean used
                ref_i_mean = mu

                #model used
                ref_i_model = kmer_model_input_file 

                ref_out = [ref_i_id, ref_i_n, ref_i_q, ref_i_sc,ref_i_loc, ref_strand, ref_i_seq, ref_i_x, ref_i_bc, ref_i_llhr, ref_i_x_is_concensus, ref_i_sig, ref_i_mean, ref_i_model]

                output_summary.loc[len(output_summary)]=ref_out

                #calculate consensus number 
                if ref_i_n >=concensus_stat_filter: 

                    if ref_i_x not in global_count:
                        global_count[ref_i_x]=0
                        
                    try: 
                        total_count[ref_i_x] = total_count[ref_i_x]+1
                    except: 
                        total_count[ref_i_x] = 1
                    consensus_count +=1

                    if ref_i_x_is_concensus ==True:
                        try: 
                            global_count[ref_i_x] = global_count[ref_i_x]+1
                        except: 
                            global_count[ref_i_x] = 1
                        consensus_count +=1
                #print(sequence[xbase_pos]+' : '+lpp[0])

            
        bar()




output_summary.to_csv(out_fn)
for xna_key in global_count: 
    print("########################################################")
    print("Xenomorph Status - [Global Morph] - Summary for ["+xna_key+"] global morph")
    print("Xenomorph Status - [Global Morph] - The number of correct consensus is (for n >"+str(concensus_stat_filter)+"): " + str(global_count[xna_key]) + "/" + str(total_count[xna_key])+' ('+str((global_count[xna_key]/total_count[xna_key])*100)+')')
    print("########################################################")




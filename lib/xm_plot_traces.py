########################################################################
########################################################################
"""
xm_plot_traces.py 

Description: This sample script uses signal-to-sequenced aligned traces 
and plots them against the canonical model and XNA model. Various options 
available for what to plot and how to plot it. Plotting options modified 
internally in script below. Randomly picks a read from a dataset to plot. 
Requires basecalled level file as input (python xenomorph.py morph).

Title: Synthesis and Sequencing of a 12-letter supernumerary DNA

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
import scipy.stats
from scipy.signal import argrelextrema
from sklearn.neighbors import KernelDensity
from sklearn.cluster import MeanShift


##############################################################
##plotting options and additional settings 
##############################################################
#Plot all raw traces 
plot_all_raw = False 

#Plot canonical model
plot_canonical_model = True

#Plot xna model
plot_xna_model =True 

#Save figure
save_figure = False
save_sub_figure = True
#Subtract canonical model
canonical_subtract = True

#Number of bases to show in place 
num_bases = 11#7 or 9 or 11

#Minimum number of reads
min_reads = 10

#This array denotes what base each XNA base should be plotted against 
confounding_pairs = ['BA','SA','PG','ZC','JC','VG','XA','KG']

#Required XNA in sequence plotted
required_xna = 'K'

#Kmer model to use
kmer_model='ATGCBSPZJVXK'

#Load kmer model based on all_bases
model_fn='active_kmer_model_sim.csv'

##Load read levels raw 
read_level_summary = ' ' 
read_level_summary= pd.read_csv(read_level_summary, sep=',')

#Kmer model mean value to use (Mean level, Median level, KDE Mean level)
mu = 'Mean level'

#KDE step size
kde_step_size =0.001


##############################################################


def gen_model(model_bases, model_fn): 

    print('Xenomorph Status - [Morph] Compiling model from input model base abbreviations.')

    #Parse model, get paths of each model file component
    mf = parse_model_files(model_bases,True)

    #Compile model from file paths into a master kmer model
    km = compile_model(mf, kmer_model) 

    #Save model 
    km.to_csv(model_fn, index=False)



def kmer2level(kmers, kmer_model,n_iter):
    level=[]

    for i in range(0,len(kmers)):
        kmer_select = kmers[i]
                
        #Fetch mean from kmer model
        kmer_s = kmer_model[kmer_model['KXmer']==kmer_select]
        mu = float(kmer_s['mu'].iloc[0]) #Fetch mean
        
       # print(kmer_model[kmer_model['KXmer']=='CKCC'])
        #Fetch standard deivation from kmer model
        noise = 0 #float(kmer_s['sigma'].iloc[0])
        
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
        orllr =100
    return orllr

def getKernelDensityEstimation(values, kde_step_size, kernel):
    x = np.arange(-5,5,kde_step_size)
    level=values
    iqr = scipy.stats.iqr(level)
    a_factor = 0.9*np.min([iqr/1.34,np.std(level)])
    bandwidth = a_factor*len(level)**(-1/5)
    if bandwidth == 0: 
        bandwidth = 0.2

    model = KernelDensity(kernel = kernel, bandwidth=bandwidth)
    model.fit(np.array(values).reshape((-1,1)))
    logprob= model.score_samples(np.array(x).reshape((-1,1)))

    clustering = MeanShift(bandwidth=bandwidth).fit(np.array(values).reshape((-1,1)))

    #Get maxima
    ma = argrelextrema(logprob, np.greater)[0]

    #Get min of max
    logma = list(logprob[ma])
    logma.sort()

    try:
        minin = list(logprob[ma]).index(logma[-rank_density])
    except: 
        minin = list(logprob[ma]).index(logma[-1])
    kernel_mean = x[ma][minin]

    return kernel_mean


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


def sub_seq(sequence, xbase_pos, sub_base):
    
    seq_alt = sequence[0:xbase_pos]+sub_base+sequence[xbase_pos+1:]


    return seq_alt

    
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
        else:
            lrr_bc = all_bases[lrr.index(np.min(lrr))]
        most_likely_base = lrr_bc
        llhr = np.min(lrr)

    return most_likely_base, llhr
    


        
 
############################
#Load kmer measurements
############################

#Generate model 
gen_model(kmer_model, model_fn)
kmer_model_input_file = (model_fn)
kmer_model_input = pd.read_csv(kmer_model_input_file, sep=',')
kmer_model_input_sim = pd.read_csv(kmer_model_input_file, sep=',')


##Sequence to plot 
sequence_to_plot = 'ATGPCAT'


#Output file name to save 
out_fn = '' #sys.argv[3]


#Load all possible bases from model file 
all_bases = list(set(''.join(kmer_model_input['KXmer'].to_list())))



##############################################
#Configure kmer model and kmer simulation model
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
##############################################
#Generate kmer data sampling from real dataset
##############################################


#Generate every possible nnn-mer
xbase_pos = nmer_model.find('x') 


base_call = [] 
llikelihood_ratio = [] 


#Random read
read_level_summary = read_level_summary[read_level_summary['read_xna_sequence'].str.contains(required_xna)]
sequence_list = list(set(read_level_summary['read_xna_sequence']))

n_count=0
while n_count <=min_reads: 
    rc = random.choice(sequence_list)
    read_level_summary_choice = read_level_summary[read_level_summary['read_xna_sequence'] == random.choice(sequence_list)]
    n_count = len(read_level_summary_choice)
read_level_summary = read_level_summary_choice 





level_per_read = []

print(read_level_summary)
with alive_bar(len(read_level_summary), force_tty=True) as bar: 
    for i in range(0,len(read_level_summary)): 




        #Setup sequence 
        sequence=read_level_summary.iloc[i]['read_xna_sequence']
        print(sequence)



        #Get read levels
        read_levels=read_level_summary.iloc[i]['read_levels']



        #break up sequence into kmer 
        kmer=[]; levels =[] 
        for k in range(0,len(kmer_mask)): 
            kmer = kmer+seq2kmer(sequence,len(kmer_mask[k]))
            levels = levels + readlevel2kmerlevel(read_levels, kmer_mask[k])





        level=[]
        for i in range(0,len(levels)): 
            try:
                level.append(float(levels[i]))
            except: 
                ex = 0 
        level.insert(0,level[0])


        #Append level 
        level_per_read.append(level)


        if plot_all_raw == True: 
            xx = np.arange(len(sequence)-2)
            sns.lineplot(xx,level, drawstyle='steps-pre', linestyle ="None", color='k')

        bar()


#X-values: X-range 
xx = np.arange(len(sequence)-2)


#Set up plotting - global formatting 
font_fam = 'Courier New'
plt.rcParams["font.family"] = font_fam
plt.rcParams.update({'font.size': 14})
plt.rcParams["figure.figsize"] = (2,5)
#sns.set_style({'grid.linestyle': '--'})


#Y-values-1: Take mean and std of the collected raw data 
if mu == 'Mean level':
    yy_global = np.mean(level_per_read,axis=0)
    yy_std = np.std(level_per_read, axis=0)
    yy_err = (yy_std[1:])
    #Plot mean+std error of raw data
    plt.figure('Zscore - '+sequence) 
    yy_shft = yy_global[1:]; yy_err = (yy_std[1:])
    ax = sns.lineplot(xx,yy_global, drawstyle='steps-pre',  color='orange')
    ax.errorbar(xx[1:]-0.5,yy_shft, yy_err, color ='orange',  linestyle = 'None', marker = 's', markersize = 5, linewidth = 2)




if mu == 'Median level': 
    yy_global = np.median(level_per_read,axis=0)
    yy_std = np.std(level_per_read, axis=0)
    yy_err = (yy_std[1:])
    #Plot mean+std error of raw data
    plt.figure('Zscore - '+sequence) 
    yy_shft = yy_global[1:]; yy_err = (yy_std[1:])
    ax = sns.lineplot(xx,yy_global, drawstyle='steps-pre',  color='orange')
    ax.errorbar(xx[1:]-0.5,yy_shft, yy_err, color ='orange',  linestyle = 'None', marker = 's', markersize = 5, linewidth = 2)





#Y-values-2: XNA model 
#print(kmer)
kmod = kmer2level(kmer, kmer_model,1)
kmod = list(np.concatenate(kmod).ravel())
kmod.insert(0,kmod[0])


#Y-values-3: Canonical model
xpos =int((len(sequence)-1)/2)
xna_base = sequence[xpos]
subkmer =  list(np.concatenate(kmer2level(seq2kmer(sub_seq(sequence,xpos,xna_base_rc(xna_base,confounding_pairs)),len(kmer_mask[0])), kmer_model,1)).ravel())
smod=subkmer #Standard model in regular line space
subkmer.insert(0,subkmer[0]) #For step plot, include the 0-point







#Optional - Calculate KDE estimates
if mu == 'KDE Mean level': 

    k_levels = []
    levelt = np.transpose(level_per_read)
    for ii in range(0,len(levelt)):
        kde_mean = getKernelDensityEstimation(levelt[ii], kde_step_size,'gaussian')
        k_levels.append(kde_mean)
    kk_global = k_levels
    kk_shft= kk_global[1:]
    yy_std = np.std(level_per_read, axis=0)
    yy_err = (yy_std[1:])
    plt.figure('KDE Mean - '+sequence) 
    ax = sns.lineplot(xx,kk_global, drawstyle='steps-pre',  color='orange', linewidth = 3)
    ax.errorbar(xx[1:]-0.5,kk_shft, yy_err, color ='orange',  linestyle = 'None', marker = 's',markeredgewidth = 3, markersize = 5, linewidth = 3)


#Plot XNA model 
if plot_xna_model ==True: 
    ax = sns.lineplot(xx,kmod, drawstyle='steps-pre',   color='r', alpha = 0.4, linestyle = '-', linewidth = 3)
    #ax = sns.lineplot(xx,kmod, drawstyle='steps-pre',   color='r', linestyle = '--', linewidth = 6)

    print(kmod) 
#Plot model from sustituted XNA 
if plot_canonical_model == True: 
    sns.lineplot(xx,subkmer, drawstyle='steps-pre',   color='k', alpha = 0.4, linestyle = '-', linewidth = 3)
    #sns.lineplot(xx,subkmer, drawstyle='steps-pre',   color='k', linestyle = '--', linewidth = 2)
ylim_low = -4
ylim_high = 4

#Formatting plot 
ax.set_ylim(ylim_low, ylim_high)
ax.set_yticks(np.arange(ylim_low,ylim_high+1))
ax.set_yticklabels(np.arange(ylim_low,ylim_high+1))
ax.tick_params(axis="y",labelsize = 12)
#plt.ylabel('Normalized signal')


ax.set_xticks(np.arange(len(sequence)-3)+0.5)
ax.set_xticks(np.arange(len(sequence)-3), minor = True)
ax.set_xticklabels(sequence[1:len(sequence)-2], va = 'top')
ax.tick_params(axis="x",direction="in", pad=-15, size = 0, labelsize = 12)
ax.set_xlim(0, len(xx)-1)
#ax.grid(which='minor', linewidth=1.5)

#Set x-axis range if desired bases is < sequence length 
if num_bases < len(xx): 
    n_base_shft = num_bases//2
    n_mid = len(xx)//2
    ax.set_xlim(n_mid-n_base_shft, n_mid+n_base_shft+1)

if save_figure == True: 
    fn = 'publication/figures/Figure 4X - Model Traces - v2/'+required_xna+'_'+sequence+'_'+str(len(level_per_read))+'.pdf'
    plt.savefig(fn)
#plt.show()





if canonical_subtract == True: 
    yy_global = yy_global - smod 


    #Plot mean+std error of raw data
    plt.figure('Canonical model subtracted - '+sequence, figsize = (3,3)) 

    yy_shft = yy_global[1:]; yy_err = (yy_std[1:])
    plt.axhline(y = 0, color = [0, 0, 0, 0.5], linestyle = '--', linewidth= 2)
    ax2 = sns.lineplot(xx,yy_global, drawstyle='steps-pre',  color="orange", linewidth = 2)
    ax2.errorbar(xx[1:]-0.5,yy_shft, yy_err, color ='orange',  linestyle = 'None', marker = 's', markersize = 5)


    #Formatting plot 
    ax2.set_yticks(np.arange(-5,6))
    ax2.set_yticklabels(np.arange(-5,6))
    ax2.tick_params(axis="y",labelsize = 12)
    plt.ylabel('<Z>-Zc')
    ax2.set_ylim(-2, 2)


    ax2.set_xticks(np.arange(len(sequence)-3)+0.5)
    ax2.set_xticks(np.arange(len(sequence)-3), minor = True)
    ax2.set_xticklabels(sequence[1:len(sequence)-2], va = 'top')
    ax2.tick_params(axis="x",direction="in", pad=-15, size = 0, labelsize = 12)
    ax2.set_xlim(0, len(xx)-1)




    if num_bases < len(xx): 
        n_base_shft = num_bases//2
        n_mid = len(xx)//2
        ax2.set_xlim(n_mid-n_base_shft, n_mid+n_base_shft+1)

   # ax.grid(which='minor', linewidth=1.5)
    plt.tight_layout()
    if save_sub_figure == True: 
        fn = 'publication/figures/Figure 2X - Deviations - v2/'+required_xna+'_'+sequence+'_'+str(len(level_per_read))+'.pdf'
        plt.savefig(fn)

plt.show()





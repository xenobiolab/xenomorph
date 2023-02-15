########################################################################
########################################################################
"""
xenosim.py

Description: Generates a toy set of levels + noise (std)  given a 
sequence of kmers as input and specified model. Various user supplied 
settings can be specified including what models to use for simulation, 
specification of signal noise (kmer-specific or global), sampled data or 
simulate data, etc. Script should be modified below and ran using: 
python xenosim.py. 

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
from string import ascii_lowercase
import scipy.stats as stats
import matplotlib.pyplot as plt
import statistics
import warnings
warnings.filterwarnings("ignore")
import multiprocessing
import random
from alive_progress import alive_bar

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


def gen_model(model_bases, model_fn): 

    print('Xenomorph Status - [Morph] Compiling model from input model base abbreviations.')

    #Parse model, get paths of each model file component
    mf = parse_model_files(model_bases,True)

    #Compile model from file paths into a master kmer model
    km = compile_model(mf, model_bases) 

    #Save model 
    km.to_csv(model_fn, index=False)

#Samples levels from a raw_kmer_file. Can be used to simulate data instead of kmer2level
#Return row_wize set of level simulations
#Sampling is optional
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



#Input a Kmer string, normalized level, and a kmer_distribution CSV file. Outputs probability of observing event (P e|x)
def kmerpdf(kmer_select,event,kmer_model):
    #Generate pdf of kmer levels with mean mu and standard deviation std
    kmer_s = kmer_model[kmer_model['KXmer']==kmer_select]
    mu =kmer_s['mu'].iloc[0]
    st = kmer_s['sigma'].iloc[0]
    return stats.norm(float(mu),float(st)).pdf(event)

 
 
#Input a path of kmers to take, and output the log probability of that path
def chart_path(kmer_path, levels, kmer_model):
    p = 0
    for i in range(0,len(kmer_path)):
        kmer_select = kmer_path[i]
        event = levels[i]
        pj = np.log10(kmerpdf(kmer_select,event,kmer_model))
        p = p + pj
    return (p)
    

    
#Strickly uses 7th base position. Given a sequence and list of levels, will calculate all alternative likelyhoods 
def gen_alt_all(sequence, kmer_levels, kmer_model, all_bases,n_iter, xbase_pos,kmer_size):
    
    #Store np array of log likelihoods 
    lpp=[]
    for j in range(0,len(all_bases)):

        #Generate alternative sequence - Uses fixed 7mer index 
        seq_alt = sequence[0:xbase_pos]+all_bases[j]+sequence[xbase_pos+1:8]

        #Convert sequence to kmer list, then calculate liklihoods 
        lp = chart_path(seq2kmer(seq_alt,kmer_size), kmer_levels, kmer_model)


        if all_bases[j] in xna_bases: 
            lp = lp+rocd

        #store likelihood output in lpp array 
        lpp.append(lp)
        


    #Generate bool table of best for each simulation
    ppl_max =  np.multiply(lpp == np.amax(lpp, axis=0),1)



    #Calculate %basecall for each alternative base - column follows all_bases order 
    if analysis_mode != 'global':
        ppl_max =  np.sum(ppl_max,axis=1)/n_iter


    #Return %accuracy for each base 
    return ppl_max
    


        
 
############################
#Load kmer measurements
############################

xna_bases = ['A','T','G','C','B','S','P','Z','X','K','J','V']
all_bases = ['A','T','G','C','B','S','P','Z','X','K','J','V']
model_bases = 'ATGCBSPZXKJV'


#Analysis mode
analysis_mode = 'per-read' #per-read or global

#Raw kmer file input (for sampling-based generation)
raw_kmer_file = ''

#Output file name to save 
out_fn = ''

#Nmer model to simulate (nnnxnnn, nxnn, nnxnn etc) 
nmer_model = 'nnnxnnn'

#Kmer models to use (e.g. [4, 5, 6] will use 4kmer, 5kmer, and 6kmer models in ML calculations 
kmer_sizes = [4] 

#Set number of iterations to run in simulation
n_iter=1000

#Kmer levels 
kmer_levels = 'Simulate' #Sample (from real data) or Simulate kmer levels from PDF

#Set mu when simulating reads as: KDE Mean level, Mean level, or Median level 
mu = 'KDE Mean level'

#Set sigma as 'Global' for global median (recommended),  'Kmer' for kmer-specific, or set as float for fixed
sigma = 0.4

#Set mu for basecalling model as: KDE Mean level, Mean level, or Median level 
mu_sim = 'KDE Mean level'

#Configure simulated levels model (if kmer_levels = 'Simulate'). Kmer sigma more realistic for simulated reads. 
sigma_sim = 0.4

#Load kmer model based on all_bases
model_fn='active_kmer_model_sim.csv'
gen_model(''.join(all_bases), model_fn)
kmer_model_input_file = (model_fn )
kmer_model_input = pd.read_csv(kmer_model_input_file, sep=',')
kmer_model_input_sim = pd.read_csv(kmer_model_input_file, sep=',')

#bases to allow as trimers before and after NNNXNNN
nnn_bases = ['A', 'T', 'G', 'C']


#KDE step size
kde_step_size = 0.01



#One solution is you have an XNA sequence then convert to regular, then sample reg kmers from that 

##############################################
#Configure kmer model and kmer simulation model
##############################################

#Model configuration

sig = kmer_model_input['Std level']

if isinstance(sigma, str):
    if sigma == 'Global':
        sig = np.median(sig)
else:
    sig = sigma
    
kmer_model_input['Std level'] = sig

kmer_model=kmer_model_input[['KXmer',mu,'Std level']]

kmer_model.rename(columns = {mu:'mu', 'Std level':'sigma'}, inplace = True)

kmer_list = kmer_model['KXmer']



if kmer_levels =='Simulate': 
    #Simulation model configuation
    sig_sim = kmer_model_input_sim['Std level']
    if isinstance(sigma_sim, str):
        if sigma_sim == 'Global':
            sig_sim = np.mean(kmer_model_input['Std level'])
    else:
        sig_sim = sigma_sim

    kmer_model_input_sim['Std level'] = sig_sim
        
    kmer_model_sim=kmer_model_input_sim[['KXmer',mu,'Std level']]

    kmer_model_sim.rename(columns = {mu:'mu', 'Std level':'sigma'}, inplace = True)


#Load raw kmer file - if sampling levels from raw kmers 
if kmer_levels == 'Sample': 
    print('[Sampling] - Loading Xenosim input model and kmer files.')
    #raw_kmer_file = 'xenomorph_pub/PZ_libv2_FLG001GC_kmers.csv'
    #raw_kmer_file = 'xenomorph_pub/ATGC_blunt_libv2_FLG001GC_kmers.csv'
    raw_kmer_samples = pd.read_csv(raw_kmer_file, sep=',', dtype={'kmer_xy': str, 'mean_level': object}).dropna(subset = ['mean_level'])




##############################################
#Generate kmer data sampling from real dataset
##############################################

#Create empty dataframe to store results 
simu_out = pd.DataFrame(columns = ['Sequence', 'NNN1', 'XNA', 'NNN2', 'Iterations'])
simu_out[all_bases]=[]


#Generate every possible nnn-mer
xbase_pos = nmer_model.find('x') 
n1 = nnn = list(set([''.join(i) for i in itertools.product(nnn_bases, repeat = xbase_pos)]))
n1.sort()
n2 = nnn = list(set([''.join(i) for i in itertools.product(nnn_bases, repeat = len(nmer_model)-xbase_pos-1)]))
n2.sort()
nnnxnnn = list(itertools.product(n1,xna_bases,n2))

tp =[]
fp =[]

print('Running basecall simulations on '+nmer_model)
with alive_bar(len(nnnxnnn), force_tty=True) as bar: 
    for i in range(0,len(nnnxnnn)): 
        sequence = nnnxnnn[i][0]+nnnxnnn[i][1]+nnnxnnn[i][2]

        #break up sequence into kmer 
        kmer=[]
        for k in range(0,len(kmer_sizes)): 
            kmer = kmer+seq2kmer(sequence,kmer_sizes[k])


       # print(kmer)
        #Generate kmer levels 
        if kmer_levels =='Sample':  
            levels = sample2levels(kmer, raw_kmer_samples, n_iter)

        if kmer_levels =='Simulate': 
            levels = kmer2level(kmer, kmer_model_sim, n_iter) 


        if analysis_mode == 'global': 
            k_levels = []
            for ii in range(0,len(levels)):
                kde_mean = getKernelDensityEstimation(levels[ii], kde_step_size,'gaussian')
                k_levels.append(kde_mean)
            #Calculate log liklihoods
            levels = k_levels

        lpp = gen_alt_all(sequence,levels, kmer_model, all_bases, n_iter,xbase_pos,kmer_sizes[0])



        #For saving into output 
        nnn1=nnnxnnn[i][0]
        x=nnnxnnn[i][1]
        nnn2=nnnxnnn[i][2]
        simu = np.concatenate([[sequence, nnn1, x, nnn2, n_iter],lpp])

        simu_out.loc[len(simu_out.index)]=simu
        bar()

#Save output 
simu_out.to_csv(out_fn)




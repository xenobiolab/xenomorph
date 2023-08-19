########################################################################
########################################################################
"""
parse_kmer.py

Description: This script parses extracted kmer output files to generate 
a kmer model file. To generate input file, run xm_extract_levels on the 
level output file generated from xenomorph.py preprocess. The kmer model
summary generated can be used directly for basecalling purposes. The location
of the model needs to be specified in models/config_model.csv

Title: Synthesis and Sequencing of 12-Letter Supernumerary DNA

By: H. Kawabe, C. Thomas, S. Hoshika, Myong-Jung Kim, Myong-Sang Kim, L. Miessner, J. M. Craig, 
J. Gundlach, A. Laszlo,  S. A. Benner, J. A. Marchand

Updated: 8/18/23
"""
########################################################################
########################################################################

import numpy as np
import pandas as pd
import itertools 
import seaborn as sns
import os
import sys
import scipy.stats
from string import ascii_lowercase
from alive_progress import alive_bar
from sklearn.cluster import MeanShift
from scipy.signal import argrelextrema
from sklearn.neighbors import KernelDensity


#from tombo import tombo_helper, tombo_stats, resquiggle


print('Xenomorph Status - [Model building] Calculating kmer model statistics')
#Read kmer level output file as an input
input_file = sys.argv[1]
output_file = sys.argv[2]
bandwidth = 0.2
step_size = 0.005
rank_density = 1 


#Get mean of maxima from KDE
def getKernelDensityEstimation(values, x, bandwidth, kernel):

    #Generate kernel density
    model = KernelDensity(kernel = kernel, bandwidth=bandwidth)
    model.fit(np.array(values).reshape((-1,1)))
    logprob= model.score_samples(np.array(x).reshape((-1,1)))


    #Cluster
    clustering = MeanShift(bandwidth=bandwidth).fit(np.array(values).reshape((-1,1)))
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




pd.set_option('display.max_colwidth', None) 
kmer_level = pd.read_csv(input_file, sep=',', dtype={'kmer_xy': str, 'mean_level': object})
kmer_level = kmer_level.drop(kmer_level.columns[0], axis=1) 
kmer_level = kmer_level.dropna(subset = ['mean_level'])


if len(sys.argv)==4: 
    kmer_level = kmer_level[kmer_level['kmer_xy'].str.contains('^[ATGC]+$')]



kmer_output=pd.DataFrame()


maout=[]
miout=[] 
lmout=[]
lsout=[]
kxout=[]
coout=[] 
mnout=[]
mkout=[] #KDE mean
nkout=[] #KDE n
skout=[] #KDE std

datalist =[] 



with alive_bar(len(kmer_level), force_tty=True) as bar: 
	for i in range(0,len(kmer_level)):

	    bar()
	    kmx = kmer_level.iloc[i]['kmer_xy']

	    lev = kmer_level.iloc[i]['mean_level'].replace('\n',' ').replace("'",'').replace(']','').replace('[','').split(' ')


	    level=[]
	    for i in range(0,len(lev)): 
		    try:
			    level.append(float(lev[i]))
		    except: 
			    ex = 0 
	    datalist.append(np.array(level))


	    #Make unique - It might be worth going back to make sure we are using unique signals in the original xombo_py.py script 
	    level = list(set(level))
	    
	    ma = np.max(level)
	    mi = np.min(level)
	    lm = np.mean(level)
	    ls = np.std(level) 
	    mn = np.median(level)
	    cov =len(level) 

	    #Kernsel Density Estimation 
	    xx = np.arange(-5,5,step_size)
	    iqr = scipy.stats.iqr(level)
	    a_factor = 0.9*np.min([iqr/1.34,ls]) #Silverman's Rule 
	    bandwidth = a_factor*len(level)**(-1/5)
	    if bandwidth == 0: 
	        bandwidth = 0.2
	    kde_mean = getKernelDensityEstimation(level, xx, bandwidth, 'gaussian')
	    level = [x for x in level if x >= (kde_mean)-0.5]
	    level = [x for x in level if x <= (kde_mean)+0.5]
	    try: 
	        kde_std = np.std(level)
	    except: 
	        kde_std = 0
	    kde_n = len(level)


	    #Write Summary 
	    maout.append(ma)
	    miout.append(mi)
	    kxout.append(kmx)
	    lsout.append(ls)
	    coout.append(cov)
	    lmout.append(lm) 
	    mnout.append(mn)
	    mkout.append(kde_mean)
	    nkout.append(kde_n)
	    skout.append(kde_std)


kmer_output['KXmer']=kxout
kmer_output['Coverage']=coout
kmer_output['Mean level']=lmout
kmer_output['Std level']=lsout 
kmer_output['Median level']=mnout
kmer_output['Max level']=maout
kmer_output['Min level']=miout
kmer_output['Median level']=mnout
kmer_output['KDE Mean level']=mkout
kmer_output['KDE Std level']=skout
kmer_output['KDE Coverage']=nkout
kmer_output.to_csv(output_file) 


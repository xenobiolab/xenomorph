
##############################################################
"""
xr_kmer_rescale.py 
Function: Calculates new scaling parameters using the Thiel-Sen estimator or other statistic. 
Rescaling is required for signal normalization performed using remora API to compare extracted levels to kmer models. 
Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. A. Marchand

Updated: 8/17/23

"""
##############################################################
##############################################################
import pandas as pd
import numpy as np
import itertools
import seaborn as sns
import os
from string import ascii_lowercase
import scipy.stats as stats
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import statistics
import warnings
from xr_params import *
##############################################################


#This is set in xr_params. Median suggested since it is more robust. 
if rescale_metric =='mean': 
    compare = 'Mean level'
elif rescale_metric =='median':
    compare = 'Median level'

stdlev = 'Std level'
#stdlev = 'Std level'




#Get reference model - Should match flow cell chemistry and use same scaling as XNA models. 
reference_model_path = 'models/libv2/ATGC_libv2_FLG001.csv'
kmer_reference_model = pd.read_csv(reference_model_path, sep = ',')

#Filter refernece model to only contain kmers with ATGC 
kmer_reference_model = kmer_reference_model[kmer_reference_model['KXmer'].str.contains('^[AGC]+$')]








x=[]
y=[]
xe =[]
ye =[]
kmer_n = []
xyd=[]


if 1<0: 
    for i in range(0,len(kmer_list)):
        #print(kmer_list)
        try: 
            #kmer_x = kmer_list[i][xpos-1]+'B'+kmer_list[i][xpos+1:]
            kmer_x = kmer_list[i].replace(base_x,base_y)
            km1= float(kmer_set_1[kmer_set_1['KXmer']==kmer_list[i]][compare])
            km2= float(kmer_set_2[kmer_set_2['KXmer']==kmer_x][compare])

            #Get standard deviations used
            st1= float(kmer_set_1[kmer_set_1['KXmer']==kmer_list[i]][stdlev])
            st2= float(kmer_set_2[kmer_set_2['KXmer']==kmer_x][stdlev])

            #Store values
            x.append(km1)
            y.append(km2)
            xe.append(st1)
            ye.append(st2)
            kmer_n.append(kmer_list[i])
            xydiff = km1 - km2
            #Cohen's d = xydiff / std
            xyd.append(xydiff)
        except:
            n=0





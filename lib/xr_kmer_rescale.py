
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
import sys
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



#Rescaling metric is set in xr_params. Median suggested since it is more robust. 
if rescale_metric =='mean': 
    rescale_level = 'Mean level'
elif rescale_metric =='median':
    rescale_level = 'Median level'

#Get raw data model estimate - Unscaled
raw_model_path = sys.argv[1]
kmer_raw_model = pd.read_csv(raw_model_path, sep = ',')

#Get reference model - Should match flow cell chemistry and use same scaling as XNA models. 
reference_model_path = rescale_reference_model_path
kmer_reference_model = pd.read_csv(reference_model_path, sep = ',')

#Filter raw and refernece model to only contain kmers with ATGC 
kmer_reference_model = kmer_reference_model[kmer_reference_model['KXmer'].str.contains('^[ATGC]+$')]
kmer_raw_model = kmer_raw_model[kmer_raw_model['KXmer'].str.contains('^[ATGC]+$')]

print('Xenomorph Status - [Rescale] ' +str(len(kmer_raw_model))+' kmers available for calculating global rescale parameters')

#Merge dataframes for estimate calculations
merged_kmers = kmer_raw_model [['KXmer', rescale_level]].merge(kmer_reference_model[['KXmer', rescale_level]], on='KXmer', how='left')
x=merged_kmers[rescale_level+'_x']
y=merged_kmers[rescale_level+'_y']



if rescale_save_plot == True: 
    fig, ax=plt.subplots()
    plt.errorbar(x,y, fmt ='o', solid_capstyle='projecting', capsize=5, color ='indigo')
    plt.xlabel('Raw kmer estimates')
    plt.ylabel('Model kmer estimates')
    plt.axline((1, 1), slope=1, color="black", linestyle="--")
    plt.axline((0.6, 1), slope=1, color="red", linestyle="--")
    plt.axline((1.4, 1), slope=1, color="red", linestyle="--")
    plt.xlim([-5, 5])
    plt.ylim([-5, 5])
    #plt.axis('square')
    ax.set_box_aspect(1)
    plt.savefig(sys.argv[2], format='pdf')

#Calculate test statistics 
if rescale_method =='Poly':
    theta = np.polyfit(x, y, 1)
    print(f'Xenomorph Status - [Rescale] Polyfit used to calculate new scaling parameters: m = '+str(theta[0])+'   b = '+str(theta[1])) 

elif rescale_method =='Thiel-Sen':
    thiel = stats.theilslopes(y, x)
    print(f'Xenomorph Status - [Rescale] Thiel-sen estimatate used to calculate new scaling parameters: m = '+str(thiel[0])+'   b = '+str(thiel[1])) 


new_rescale = thiel[0]
new_reshift = thiel[1]


print(f'Xenomorph Status - [Rescale] Updating global scaling and global shift in lib/xr_params.py.') 
# Read the content of params.txt
with open('lib/xr_params.py', 'r') as file:
    lines = file.readlines()

# Update the 'rescale' parameter in the content
for i, line in enumerate(lines):
    if line.startswith('global_rescale'):
        lines[i] = f'global_rescale = {new_rescale}\n'
        
    if line.startswith('global_reshift'):
        lines[i] = f'global_reshift = {new_reshift}\n'


# Write the updated content back to params.txt
with open('lib/xr_params.py', 'w') as file:
    file.writelines(lines)



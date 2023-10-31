########################################################################
########################################################################
"""
xr_gen_null.py 

Used to generate null call files from ATGC level files 


Cite us: 

H. Kawabe, C. Thomas, S. Hoshika, Myong-Jung Kim, Myong-Sang Kim, L. Miessner, N. Kaplan, J. M. Craig, 
J. Gundlach, A. Laszlo,  S. A. Benner, J. A. Marchand. "Enzymatic Synthesis and Nanopore Sequencing of 12-Letter Supernumerary DNA" 
Nature Communications. 14. (2023). DOI: 10.1038/s41467-023-42406-z 

Updated: 7/2/23

"""

#######################################################################



#import everything
import re
import os
import sys
import subprocess

from alive_progress import alive_bar; import time
import pod5
import pysam
import numpy as np
import pandas as pd
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

from remora import io, refine_signal_map, util
from time import sleep
from tqdm import tqdm
from csv import writer
from Bio import SeqIO
from Bio.Seq import Seq
from xm_tools import *
from xr_tools import *
from xm_params import *
from xr_params import *



#Handle input arguments 
##Pod5 file input is handled by coverting raw fast5 to pod5 using pod5 tools
read_level_file = sys.argv[1]
read_level = pd.read_csv(read_level_file, sep=',')

print('Xenomorph [Status] - Generating null read file by swapping using the following rule set:')
hardswap = ['GP','CZ','A-', 'T-']
print(hardswap)


#Loop through reads, set up bar
with alive_bar(len(read_level), force_tty=True) as bar: 
    for i in range(0,len(read_level)): 
        bar()

        #Check for all hard swap bases 
        for h in range(0,len(hardswap)): 

            #Check for hard swap base present in the focus position column
            if read_level.iloc[i]['read_xna']==hardswap[h][0]:

                #Swap base
                read_level.at[i,'read_xna'] = hardswap[h][1]

                #Currently hard coded to the heptamer level
                read_level.at[i,'read_xna_sequence'] = read_level.at[i,'read_xna_sequence'][0:3]+hardswap[h][1]+read_level.at[i,'read_xna_sequence'][4:]

#Drop all rows demarked by a '-' swap
read_level = read_level[read_level['read_xna'] != '-']

#Reindex after dropping 
read_level= read_level.reset_index(drop=True)

#Save read level file
print('Xenomorph [Status] - Saving null level file to'+read_level_file.replace('levels.csv','null_levels.csv'))
read_level.to_csv(read_level_file.replace('levels.csv','null_levels.csv'))
print('Xenomorph [Status] - Done generating null set from ATGC only reads.')


########################################################################
########################################################################
"""
kmer_convert.py

Title: Unpublished work

By: J. A. Marchand

Converts Xenomorph kmer table to ONT kmer table 

Updated: 8/5/23
"""
########################################################################
########################################################################

import pandas as pd
import os
import glob
import pysam
from pathlib import Path




kmer_file = '/home/marchandlab/xenomorph/xenomorph-xemora/models/libv2/P_libv2_FLG001.csv'
kmers = pd.read_csv(kmer_file)
level_type = 'Mean level'


#/home/marchandlab/xenomorph/xenomorph-xemora/models/xemora/ATGC_mean_9.4.1.csv



kmer_convert_file = open('kmers/4mer_9.4.1_PZ.csv','w')

for i in range(0,len(kmers)): 
    content = kmers.iloc[i]['KXmer']+' '+str(kmers.iloc[i][level_type])
    print(content)
    kmer_convert_file.write((content)+'\n')


kmer_convert_file.close()


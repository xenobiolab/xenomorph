########################################################################
########################################################################
"""
resquiggle.py

Description: This script is part of the Xombo file handling package. Used
to run resquiggle (Tombo) analysis on a set of reads. Automatically called
in the Xenomorph pipeline. 
For more information: python xombo.py resquiggle -h

Title: XNA tailing enables nanopore sequencing of a 12-letter genetic code

By: H. Kawabe, C. Thomas, A. Laszlo, S. Hoshika, L. Miessner, J. M. Craig, 
J. Gundlach, Myong-Jung Kim, Myong-Sang Kim, S. A. Benner, J. A. Marchand

Updated: 1/15/23
"""
########################################################################
########################################################################


import os 
import sys 
from xm_params import *

#Inputs 
input_folder = sys.argv[1].replace('//','/')
reference_fasta = sys.argv[2].replace('//','/')

if len(sys.argv)==4: 
    corr_group_name= sys.argv[3]
else: 
    corr_group_name=''

########################################################################
#Parameters for resquiggle
max_scaling_iterations = 20 #(default 20) 
signal_matching_score  = 6 #(default 6) 
max_read_length = 5000 #(default 5000) 
########################################################################

###NOTE PARAMETERS BEING USED - IN BETA 

print('Xenomorph Status - [Preprocess] Resquiggling fast5 reads with '+corr_group_name+'-based segmenation. Resquiggle stored as a meta file in input directory.') 
corr_group = 'RawGenomeCorrected_000'+corr_group_name

try: 

    exec ='tombo resquiggle '+input_folder+' '+reference_fasta+' --processes 8 --max-scaling-iterations '+str(max_scaling_iterations)+' --signal-matching-score '+str(signal_matching_score)+' --overwrite  --ignore-read-locks --num-most-common-errors 3 --q-score '+str(qscore_filter)+' --sequence-length-range 0 '+str(max_read_length)+'  --include-event-stdev --corrected-group '+corr_group+' --signal-align-parameters '+signal_align_params+' --segmentation-parameters '+segmentation_params
#print(exec)
    os.system(exec) 


except: 
    print('Error: Resquiggle failed. Check error logs and ensure tombo is installed properly.') 





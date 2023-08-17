########################################################################
########################################################################
"""
xemora_pipe.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. A. Marchand

Updated: 3/2/23
"""
########################################################################
########################################################################


import os
import glob
import sys
from pathlib import Path
from lib.xr_tools import *
from lib.xr_params import *



############################################################
#Training paths
#working_dir = '~/xenomorph/xemora-beta/xemora-test/230309_Xemora_ZC_train' 
#xna_fast5_dir = '~/xenomorph/xemora-beta/xemora-test/230309_Xemora_PZ_GC_train/230308_PZ_Xemora_train_110_fast5s'
#xna_ref_fasta = 'ref/P_Xemora_train.fa'
#dna_fast5_dir = '~/xenomorph/xemora-beta/xemora-test/230309_Xemora_PZ_GC_train/230308_GC_Xemora_train_110_fast5s'
#dna_ref_fasta = 'ref/G_Xemora_train.fa'

working_dir = '~/xenomorph/xemora-beta/xemora-test/230523_signal_extract' 
xna_fast5_dir = '~/xenomorph/xemora-beta/xemora-test/230315_PZGClibv4/230304_PZ_libv4_1k'
xna_ref_fasta = 'ref/ref_libv2_PZ_CxDx-.fa'
dna_fast5_dir = '~/xenomorph/xemora-beta/xemora-test/230315_PZGClibv4/230310_GC_libv4_1k'
dna_ref_fasta = 'ref/ref_libv2_PZ_CxDx-.fa'
ext_dir = '~/xenomorph/xemora-beta/xemora-test/230309_Xemora_PZ_GC_train/chunks/training_chunks.npz'
############################################################

############################################################
#Base calling paths 
#bc_working_dir = '~/xenomorph/xemora-beta/xemora-test/230309_Xemora_ZC_train/z-basecall' 
#bc_fast5_dir = '~/xenomorph/xemora-beta/xemora-test/230309_Xemora_PZ_GC_train/PZ_test'
#bc_xna_ref_fasta = 'ref/P_Xemora_train.fa'
#bc_model_file = '~/xenomorph/xemora-beta/xemora-test/230309_Xemora_ZC_train/model/model_best.pt'

bc_working_dir = '~/xenomorph/xemora-beta/xemora-test/230315_PZGClibv4/g-basecall-3_extval' 
bc_fast5_dir = '~/xenomorph/xemora-beta/xemora-test/230309_Xemora_PZ_GC_train/GC_Test'
bc_xna_ref_fasta = 'ref/G_Xemora_train.fa'
bc_model_file = '~/xenomorph/xemora-beta/xemora-test/230315_PZGClibv4/model/model_best.pt'

############################################################

############################################################
train_model = True
basecall_reads = False
############################################################

#conda activate nanopore-re

#Train dataset with xemora train
if train_model ==True: 
    cmd = 'python xemora.py train -w '+working_dir+' -f '+xna_fast5_dir +' '+dna_fast5_dir+' -r '+xna_ref_fasta+' '+dna_ref_fasta
    os.system(cmd)


#Basecall fast5 directory 
if basecall_reads==True: 
    cmd = 'python xemora.py basecall -w '+bc_working_dir+' -f '+bc_fast5_dir+' -r '+bc_xna_ref_fasta+' -m '+bc_model_file 
    os.system(cmd)


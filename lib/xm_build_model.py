########################################################################
########################################################################
"""
xm_build_model.py

Description: This tool is used to build kmer models using Tombo's API. 
Since Tombo API does not let you specify various quality control features, 
this tool is not used in the publication. Nonetheless, this is a useful tool 
since it can quickly generate models. The script to generate kmer models
using tombo needs to be edited below. 


Title: Synthesis and Sequencing of 12-Letter Supernumerary DNA

By: H. Kawabe, C. Thomas, S. Hoshika, Myong-Jung Kim, Myong-Sang Kim, L. Miessner, J. M. Craig, 
J. Gundlach, A. Laszlo,  S. A. Benner, J. A. Marchand

Updated: 2/14/23
"""
########################################################################
########################################################################


import os
import tombo
from Bio.Seq import Seq
import sys
import numpy as np
import pandas as pd
import itertools
import os
import h5py
from tombo import tombo_helper, tombo_stats, resquiggle
from xm_params import *
from xm_tools import *


#Input fast5 single directory 
fast5_dir = ''

#Model XNA base abbreviation 
model_base = 'B'

#Substitution base for XNA - required to manually enter it here 
sub_base='A'

#Output model filename  (model.csv)
model_name = 'XNA_FLG001_Model.csv'

#Input reference xfasta (with standard bases in the sequence, and XNA position in the header) 
reference_fasta = ''

#Tombo corrected group name 
corrected_group = 'RawGenomeCorrected_000'






def build_model(fast5_dir, model_base, model_name, reference_fasta, corrected_group, sub_base):

    #Generate a bed file 
    bed_file_name = reference_fasta.replace('.fa','.bed')
    bed_cmd = 'python lib/xm_xfasta2bed.py '+reference_fasta+' '+model_base
    os.system(bed_cmd)


    kmer_output_file = 'models/'+model_name+'.csv'
    #substitution_base = xna_base_rc(model_base,confounding_pairs)+':1'
    substitution_base = sub_base+':1'


    #Estimate model parameters
    print('Xenomorph Status - [Model Building] Using Tombo to extract Kmer specific motifs for alt model estimation')
    cmd = 'tombo build_model estimate_motif_alt_reference --alternate-model-filename '+model_name+' --alternate-model-name '+model_base+' --motif-description '+substitution_base+' --fast5-basedirs '+fast5_dir+' --corrected-group '+corrected_group+' --valid-locations-filename '+bed_file_name+' --minimum-test-reads 1  --processes 8 --upstream-bases 1 --downstream-bases 2'
    os.system(cmd)

    #Extract estimates 
    print('Xenomorph Status - [Model Building] Extracting kmer values for use with xenomorph')
    modalt=h5py.File(model_name, 'r')
    kxout=[] 
    mout=[]
    sout=[]

    for line in modalt['model']:
        kmer = str(line[0].decode("utf-8"))
        pos = line[1]
        kmer=kmer[0:pos]+model_base+kmer[pos+1:]
        kxout.append(kmer)
        mout.append(line[2])
        sout.append(line[3])


    kmer_output=pd.DataFrame()
    kmer_output['KXmer']=kxout
    kmer_output['Median level']=mout
    kmer_output['Std level']=sout 
    kmer_output.to_csv(kmer_output_file)

    #Completion
    print('Xenomorph Status - [Model Building] Kmer model saved in '+kmer_output_file)


#Execute model building
build_model(fast5_dir, model_base, model_name, reference_fasta, corrected_group, sub_base)


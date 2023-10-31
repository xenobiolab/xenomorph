########################################################################
########################################################################
"""
test.py

Tests proper installation and function of Xenomorph tools. 
Uses sample reads and sample reference file found in sample/ directory. 
Results output to test/results

See example usage: 
python test.py


Cite us: 

H. Kawabe, C. Thomas, S. Hoshika, Myong-Jung Kim, Myong-Sang Kim, L. Miessner, N. Kaplan, J. M. Craig, 
J. Gundlach, A. Laszlo,  S. A. Benner, J. A. Marchand. "Enzymatic Synthesis and Nanopore Sequencing of 12-Letter Supernumerary DNA" 
Nature Communications. 14. (2023). DOI: 10.1038/s41467-023-42406-z 

Updated: 9/07/23
"""
########################################################################
########################################################################


import os
import sys


### PZa Dataset Testing new model that Hinako made
wdir = 'sample/results'
fast5 = 'sample/fast5'
ref = 'sample/reference/PZ_sample_reference.fa'
out_file_prefix = '/sample'#Output file parameters

################################################
#Level generation
run_preprocess = True
run_null_gen = False

#Basecalling
run_null_level = False
run_morph = True
run_stats = True
run_global_morph = True

#Generate model from a level file
run_model_gen = False

#Morph model
model = 'ATGCPZ'
################################################

#0 Set up File prefixes 
out_pre = out_file_prefix+'_levels.csv'
out_bcpr = out_file_prefix+'_basecalls_per-read.csv'
out_bcglobal = out_file_prefix+'_basecalls_per-sequence.csv'


#1. Preprocess file using xemora preprocessto get a level file 
if run_preprocess == True:
    cmd = 'python xenomorph.py preprocess -w '+wdir+' -x -b -f '+fast5+' -r '+ref+' -o '+out_pre
    os.system(cmd)

#1.5 Generate a null dataset if required
if run_null_gen == True: 
    cmd = 'python lib/xr_gen_null.py '+ wdir+out_pre
    os.system(cmd) 

#1.5 Swap level file name for the null version. This will make downstream null
if run_null_level ==True: 
    out_pre = out_pre.replace('levels.csv','null_levels.csv')

#2. Testing whether we can run morph command following level extraction
if run_morph == True:
    cmd = 'python xenomorph.py morph -l '+wdir+out_pre+' -m '+model+' -o '+wdir+out_bcpr
    os.system(cmd) 

#3. Stats on a morph run 
if run_stats == True:
    #cmd = 'python lib/xm_stats.py  '+wdir+out_bcpr+' '+wdir+'/stats_log.csv' 
    cmd = 'python xenomorph.py stats -i '+wdir+out_bcpr
    os.system(cmd)


#4 Global morph (sigal averaged )
if run_global_morph == True:
    cmd = 'python xenomorph.py morph -l '+wdir+out_pre+' -m '+model+' -g -o '+wdir+out_bcglobal
    os.system(cmd) 

#5 Model generation commands - Two steps
if run_model_gen== True:
    #First, extract all kmer signals
    cmd = 'python lib/xm_extract_levels.py '+wdir+out_pre
    os.system(cmd) 

    #Take statistics on each kmer set to generate a model file. Mean, median, and a KDE mean estimate are generated. 
    cmd = 'python lib/parse_kmer.py '+wdir+out_pre.replace('_levels.csv','_kmers.csv')+' '+wdir+out_pre.replace('_levels.csv','_summary.csv')
    os.system(cmd) 







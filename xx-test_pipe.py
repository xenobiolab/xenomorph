import os


#pipeline to test the integration of xemora preprocessing into xenomorph

################################################

#directories (working, fast5,fasta)
wdir = '/home/marchandlab/xenomorph/xenomorph-xemora/xx-test/230628_PZ_libv2_INTcheck'
fast5 = '/home/marchandlab/xenomorph/xenomorph-xemora/xx-test/230628_PZ_libv2_INTcheck/230304_PZ_libv4_1k_fast5s'
ref = '/home/marchandlab/xenomorph/xenomorph-xemora/ref/ref_libv2_PZ_CxDx.fa'

#Output file parameters
out_file_prefix = 'FB_P_All'

#Morph model
model = 'ATGCPZ'

out_pre = out_file_prefix+'_FLG001'+'_levels.csv'
out_bcpr = out_file_prefix+'_FLG001_bc_global.csv'

################################################
run_preprocess = True
run_morph =False
run_stats = False

#1. Preprocess file using xemora preprocessto get a level file 
if run_preprocess == True:
    cmd = 'python xenomorph.py preprocess -w '+wdir+' -f '+fast5+' -r '+ref+' -o '+out_pre
    os.system(cmd)

#2. Testing whether we can run morph command following level extraction
if run_morph == True:
    cmd = 'python xenomorph.py morph -l '+wdir+out_pre+' -m '+model+' -o '+wdir+out_bcpr
    os.system(cmd) 

#3. Stats on a morph run 
if run_stats == True:
    cmd = 'python lib/xm_stats.py  '+wdir+out_bcpr+' '+wdir+'stats_log.csv'
    os.system(cmd) 




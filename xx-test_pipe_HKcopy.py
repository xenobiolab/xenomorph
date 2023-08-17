import os


#pipeline to test the integration of xemora preprocessing into xenomorph

################################################

#directories (working, fast5,fasta)
#wdir = '/home/marchandlab/xenomorph/xenomorph-xemora/xx-test/230628_PZ_libv2_INTcheck'
#fast5 = '/home/marchandlab/xenomorph/xenomorph-xemora/xx-test/230628_PZ_libv2_INTcheck/230304_PZ_libv4_1k_fast5s'
#ref = '/home/marchandlab/xenomorph/xenomorph-xemora/ref/ref_libv2_PZ_CxDx.fa'



#
wdir = '/home/marchandlab/xenomorph/xenomorph-xemora/xx-test/230702_PZ_libv2_full'
fast5 = '/home/marchandlab/DataAnalysis/Marchand/230125_PZ_libv2_AB/20230125_1810_MN41475_ANL798_9c8050d8/fast5'
ref = '/home/marchandlab/xombo/reference_sequences/ref_libv2_PZ_AxBx.fa'


#PZ Model Building
wdir = '/home/marchandlab/xenomorph/xenomorph-xemora/xx-test/230702_PZ_libv2_full'
fast5 = '/home/marchandlab/DataAnalysis/Marchand/230125_PZ_libv2_AB/20230125_1810_MN41475_ANL798_9c8050d8/fast5'
ref = '/home/marchandlab/xombo/reference_sequences/ref_libv2_PZ_AxBx.fa'

#JV Model Building
wdir = '/home/marchandlab/xenomorph/xenomorph-xemora/xx-test/230702_JV_libv2_full'
fast5 = '/home/marchandlab/DataAnalysis/Marchand/221216_JV_libv2/20221216_1806_MN41475_ANQ772_516c93a1/fast5'
ref = '/home/marchandlab/xombo/reference_sequences/ref_libv2_JV_CxDx-.fa'



#ATGC Model Building
wdir = '/home/marchandlab/xenomorph/xenomorph-xemora/xx-test/230702_ATGC_libv2_full'
fast5 = '/home/marchandlab/DataAnalysis/Marchand/221122_blunt_libv2/'
ref = '/home/marchandlab/xombo/reference_sequences/ref_full_ABxAB.fa'
out_file_prefix = '/FB_ATGC_All'#Output file parameters



#PZ Model Building
wdir = '/home/marchandlab/xenomorph/xenomorph-xemora/xx-test/230703_XPCR'
fast5 = '/home/marchandlab/DataAnalysis/Marchand/230626_XPCR_minION_1M/20230626_1856_MN41475_FAU81971_ef846d69/fast5'
ref = '/home/marchandlab/xombo/reference_sequences/XPCR_full.fasta'
out_file_prefix = '/XPCR_PZ_All' #Output file parameters



#Morph model
model = 'ATGCPZ'

out_pre = out_file_prefix+'_FLG001'+'_levels.csv'
out_bcpr = out_file_prefix+'_FLG001_bc.csv'

################################################
run_preprocess = True
run_morph = False
run_stats = False
run_model_gen = False

#1. Preprocess file using xemora preprocessto get a level file 
if run_preprocess == True:
    cmd = 'python xenomorph.py preprocess -w '+wdir+' -f '+fast5+' -r '+ref+' -o '+out_pre
    os.system(cmd)

#2. Testing whether we can run morph command following level extraction
if run_morph == True:
    cmd = 'python xenomorph.py morph -l '+wdir+out_pre+' -m '+model+' -o '+wdir+'/'+out_bcpr
    print(cmd)
    os.system(cmd) 

#3. Stats on a morph run 
#Does not work because headers are slightly different 
if run_stats == True:
    cmd = 'python lib/xm_stats.py  '+wdir+out_bcpr+' '+wdir+'/stats_log.csv'
    os.system(cmd) 


#4. If ge
if run_model_gen== True:
    cmd = 'python lib/xm_extract_levels.py '+wdir+out_pre
    os.system(cmd) 


    cmd = 'python lib/parse_kmer.py '+wdir+out_pre.replace('_levels.csv','_kmers.csv')+' '+wdir+out_pre.replace('_levels.csv','_summary.csv')
    os.system(cmd) 




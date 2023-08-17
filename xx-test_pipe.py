import os

#Workflow integration for Xenomorph-Xemora

#JV Dataset 1 - Model Building
wdir = '/home/marchandlab/xenomorph/xenomorph-xemora/xx-test/230702_JV_libv2_full'
fast5 = '/home/marchandlab/DataAnalysis/Marchand/221216_JV_libv2/20221216_1806_MN41475_ANQ772_516c93a1/fast5'
ref = '/home/marchandlab/xombo/reference_sequences/ref_libv2_JV_CxDx-.fa'


#PZ Dataset 1 - Model Building
wdir = '/home/marchandlab/xenomorph/xenomorph-xemora/xx-test/230702_PZ_libv2_GC'
fast5 = '/home/marchandlab/DataAnalysis/Marchand/230125_PZ_libv2_AB/20230125_1810_MN41475_ANL798_9c8050d8/fast5'
ref = '/home/marchandlab/xombo/reference_sequences/ref_libv2_PZ_AxBx.fa'
out_file_prefix = '/FB_PZ_GC'#Output file parameters


###ATGC Blunt Dataset 1 -  Model Building 
wdir = '/home/marchandlab/xenomorph/xenomorph-xemora/xx-test/230702_ATGC_blunt2_libv2_full'
fast5 = '/home/marchandlab/DataAnalysis/Marchand/221122_blunt_libv2/'
ref = '/home/marchandlab/xombo/reference_sequences/ref_full_ABxAB.fa'
out_file_prefix = '/FB_ATGC_All'#Output file parameters


###ATGC Blunt Dataset 2 --Model testing
wdir = '/home/marchandlab/xenomorph/xenomorph-xemora/xx-test/ATGC_Model_Testing'
fast5 = '/home/marchandlab/DataAnalysis/Marchand/221215_blunt_libv2/20221215_1717_MN41475_ANV913_5ca08b71/fast5'
ref = '/home/marchandlab/xombo/reference_sequences/ref_full_ABxAB.fa'
out_file_prefix = '/ATGC_Model_Testing_67'#Output file parameters



###ATGC Blunt Dataset 2 --Model testing
wdir = '/home/marchandlab/xenomorph/xenomorph-xemora/xx-test/ATGC_Model_Testing'
fast5 = '/home/marchandlab/DataAnalysis/Marchand/221215_blunt_libv2/20221215_1717_MN41475_ANV913_5ca08b71/fast5'
ref = '/home/marchandlab/xombo/reference_sequences/ref_full_ABxAB.fa'
out_file_prefix = '/ATGC_Model_Testing_Rescale'#Output file parameters



### PZ Dataset 2 -- Model testing
wdir = '/home/marchandlab/Dev/xenomorph-xemora/xx-test/PZ_Model_Testing'
fast5 = '/home/marchandlab/DataAnalysis/Marchand/221124_PZ_libv2_200k/20221124_1855_MN37138_AMY698_27041bd8/fast5'
ref = '/home/marchandlab/Dev/xombo/reference_sequences/ref_libv2_PZ_CxDx.fa'
out_file_prefix = '/PZ_Model_Testing'#Output file parameters


################################################
#Level generation
run_preprocess = True
run_null_gen = False

#Basecalling
run_null_level = False
run_morph = False
run_stats = False
run_global_morph = False 

#Generate model from a level file
run_model_gen = False


#Morph model
model = 'ATGCPZ'
################################################



#0 Set up File prefixes 
out_pre = out_file_prefix+'_FLG001'+'_levels.csv'
out_bcpr = out_file_prefix+'_FLG001_bc.csv'
out_bcglobal = out_file_prefix+'_FLG001_bc_global.csv'

#1. Preprocess file using xemora preprocessto get a level file 
if run_preprocess == True:
    cmd = 'python xenomorph.py preprocess -w '+wdir+' -f '+fast5+' -r '+ref+' -o '+out_pre
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
    cmd = 'python xenomorph.py morph -l '+wdir+out_pre+' -m '+model+' -o '+wdir+'/'+out_bcpr
   # print(cmd)
    os.system(cmd) 

#3. Stats on a morph run 
if run_stats == True:
    cmd = 'python lib/xm_stats.py  '+wdir+out_bcpr+' '+wdir+'/stats_log.csv'
    os.system(cmd) 

#4 Global morph (sigal averaged )
if run_global_morph == True:
    cmd = 'python xenomorph.py morph -l '+wdir+out_pre+' -m '+model+' -g -o '+wdir+'/'+out_bcglobal
    #print(cmd)
    os.system(cmd) 


#5 Model generation commands - Two steps
if run_model_gen== True:
    #First, extract all kmer signals
    cmd = 'python lib/xm_extract_levels.py '+wdir+out_pre
    os.system(cmd) 

    #Take statistics on each kmer set to generate a model file
    cmd = 'python lib/parse_kmer.py '+wdir+out_pre.replace('_levels.csv','_kmers.csv')+' '+wdir+out_pre.replace('_levels.csv','_summary.csv')
    os.system(cmd) 




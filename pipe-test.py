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


#PZ Dataset 1 - Model Building
wdir = '/home/marchandlab/Dev/xenomorph-xemora/xx-test/230702_PZ_libv2_GC'
fast5 = '/home/marchandlab/DataAnalysis/Marchand/230125_PZ_libv2_AB/20230125_1810_MN41475_ANL798_9c8050d8/fast5'
ref = '/home/marchandlab/Dev/xombo/reference_sequences/ref_libv2_PZ_AxBx.fa'
out_file_prefix = '/FB_PZ_GC'#Output file parameters


###ATGC Blunt Dataset 2 --Model testing
wdir = '/home/marchandlab/Dev/xenomorph-xemora/xx-test/ATGC_Model_Testing'
fast5 = '/home/marchandlab/DataAnalysis/Marchand/221215_blunt_libv2/20221215_1717_MN41475_ANV913_5ca08b71/fast5'
ref = '/home/marchandlab/Dev/xombo/reference_sequences/ref_full_ABxAB.fa'
out_file_prefix = '/ATGC_Model_Testing_Rescale'#Output file parameters


### PZ Dataset 2 -- Xenomorph testing on 10k reads- Tombo
wdir = '/home/marchandlab/Dev/xenomorph-xemora/xx-test/PZ_Model_Testing_Tombo'
fast5 = '/home/marchandlab/DataAnalysis/Marchand/221124_PZ_libv2_10k/fast5'
ref = '/home/marchandlab/Dev/xombo/reference_sequences/ref_libv2_PZ_CxDx.fa'
out_file_prefix = '/PZ_Model_Testing'#Output file parameters

### PZ Dataset 2 -- Model testing (Remora testing)
wdir = '/home/marchandlab/Dev/xenomorph-xemora/xx-test/PZ_Model_Testing'
fast5 = '/home/marchandlab/DataAnalysis/Marchand/221124_PZ_libv2_200k/20221124_1855_MN37138_AMY698_27041bd8/fast5'
ref = '/home/marchandlab/Dev/xombo/reference_sequences/ref_libv2_PZ_CxDx.fa'
out_file_prefix = '/PZ_Model_Testing'#Output file parameters

### PZa Dataset Testing new model that Hinako made
wdir = '/home/marchandlab/Dev/xenomorph-xemora/xx-test/PZa_Model_Testing'
fast5 = '/home/marchandlab/DataAnalysis/Kawabe/230822_PZa_libv4_FLG001/20230822_1101_MN41475_APU320_58f6c3e2/fast5'
ref = '/home/marchandlab/Dev/xombo/reference_sequences/ref_libv2_PZ_CxDx.fa'
out_file_prefix = '/PZa_Model_Testing'#Output file parameters

###  r10.4.1 PZn model build  (Done) 
wdir = '/home/marchandlab/Dev/xenomorph-xemora/xx-test/PZn_r10.4.1_Model_Building'
fast5 = '/home/marchandlab/DataAnalysis/Kawabe/230725_PZ_lib_v4_r10/20230725_1220_MN37138_APH167_a204cb54/fast5'
ref = '/home/marchandlab/DataAnalysis/Kawabe/ref/ref_libv2_PZ_CxDx-.fa'
out_file_prefix = '/PZn_Model_Building'#Output file parameters

###  r10.4.1 XK model build (done)
wdir = '/home/marchandlab/Dev/xenomorph-xemora/xx-test/XK_r10.4.1_Model_Building'
fast5 = '/home/marchandlab/DataAnalysis/Kawabe/230801_XK_libv4/20230801_1658_MN41475_APH278_da76761b/fast5'
ref = '/home/marchandlab/DataAnalysis/Kawabe/ref/ref_libv2_XK_CxDx-.fa'
out_file_prefix = '/XK_Model_Building'#Output file parameters

###  r10.4.1 JV model build  (Needed)
wdir = '/home/marchandlab/Dev/xenomorph-xemora/xx-test/JV_r10.4.1_Model_Building'
fast5 = '/home/marchandlab/DataAnalysis/Kawabe/230801_JV_libv4/20230801_1636_MN37138_APG963_e8c6c162/fast5'
ref = '/home/marchandlab/DataAnalysis/Kawabe/ref/ref_libv2_JV_CxDx-.fa'
out_file_prefix = '/JV_Model_Building'#Output file parameters

###  r10.4.1 PZa model build (Needed) 
wdir = '/home/marchandlab/Dev/xenomorph-xemora/xx-test/PZa_r10.4.1_Model_Building'
fast5 = '/home/marchandlab/DataAnalysis/Kawabe/230821_PZa_libv4_FLG114/20230821_2012_MN37138_APP955_7fa44073/fast5'
ref = '/home/marchandlab/DataAnalysis/Kawabe/ref/ref_libv2_PZ_CxDx-.fa'
out_file_prefix = '/PZa_Model_Building'#Output file parameters


### PZa Dataset Testing new model that Hinako made
wdir = '/home/marchandlab/Dev/xenomorph-xemora/xx-test/PZa_Model_Testing'
fast5 = '/home/marchandlab/DataAnalysis/Kawabe/230822_PZa_libv4_FLG001/20230822_1101_MN41475_APU320_58f6c3e2/fast5'
ref = '/home/marchandlab/Dev/xombo/reference_sequences/ref_libv2_PZ_CxDx.fa'
out_file_prefix = '/PZa_Model_Testing'#Output file parameters



### ATGC 10.4.1 model building 
wdir = '/home/marchandlab/Dev/xenomorph-xemora/xx-test/ATGC_10.4.1_Model_Building'
fast5 = '/home/marchandlab/DataAnalysis/Kawabe/231006_blunt_libv4_r10/20231005_1245_MN37138_AQL738_2675cb8f/fast5'
ref = '/home/marchandlab/DataAnalysis/Kawabe/ref/ref_full_CDxCD.fa'
out_file_prefix = '/ATGC_10.4.1_Model_Building'#Output file parameters

### BSn 10.4.1 model building 
wdir = '/home/marchandlab/Dev/xenomorph-xemora/xx-test/BSn_10.4.1_Model_Building'
fast5 = '/home/marchandlab/DataAnalysis/Kawabe/231003_BSn_libv4_FLG114/20231003_1555_MN37138_AQK018_63a7330b/fast5'
ref = '/home/marchandlab/DataAnalysis/Kawabe/ref/ref_libv2_BS_CxDx-.fa'
out_file_prefix = '/BSn_10.4.1_Model_Building'#Output file parameters


###  r10.4.1 PZn model build  (Done) 
wdir = '/home/marchandlab/Dev/xenomorph-xemora/xx-test/PZn_r10.4.1_Model_Building'
fast5 = '/home/marchandlab/DataAnalysis/Kawabe/230725_PZ_lib_v4_r10/20230725_1220_MN37138_APH167_a204cb54/fast5'
ref = '/home/marchandlab/DataAnalysis/Kawabe/ref/ref_libv2_PZ_CxDx-.fa'
out_file_prefix = '/PZn_Model_Building'#Output file parameters



### ATGC 10.4.1 model building 
wdir = '/home/marchandlab/Dev/xenomorph-xemora/xx-test/ATGC_10.4.1_Model_Building'
fast5 = '/home/marchandlab/DataAnalysis/Kawabe/231006_blunt_libv4_r10/20231005_1245_MN37138_AQL738_2675cb8f/fast5'
ref = '/home/marchandlab/DataAnalysis/Kawabe/ref/ref_full_CDxCD.fa'
out_file_prefix = '/ATGC_Null_Benchmarking'#Output file parameters

################################################
#Level generation
run_preprocess = True
run_null_gen = True

#Basecalling
run_null_level = True
run_morph = True
run_stats = True
run_global_morph = True

#Generate model from a level file
run_model_gen = False

#Morph model
model = 'ATGCPZ'
flowcell = 'FLG114'
################################################

#0 Set up File prefixes 
out_pre = out_file_prefix+'_'+flowcell+'_levels.csv'
out_bcpr = out_file_prefix+'_'+flowcell+'_bc.csv'
out_bcglobal = out_file_prefix+'_'+flowcell+'_bc_global.csv'


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
    cmd = 'python xenomorph.py morph -l '+wdir+out_pre+' -m '+model+' -o '+wdir+out_bcpr
    os.system(cmd) 

#3. Stats on a morph run 
if run_stats == True:
    cmd = 'python lib/xm_stats.py  '+wdir+out_bcpr+' '+wdir+'/stats_log.csv'
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







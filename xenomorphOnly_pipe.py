
import os

'''


wdir = 'xenomorph_testing'

fast5 = '~/DataAnalysis/Marchand/220805_BS_libv2_extracted'
ref = '~/xombo/reference_sequences/ref_full_xm.fa' 
out_pre = '220814_output_bc_ATGCBS_test.csv' 
out_bc = 'xenomorph_testing/BS_FLG002_4K-Global.csv' 
model = 'models/ATGCBS_FLG002.csv'
model = 'models/ATGCBS_FLG002.csv'




wdir = 'xenomorph_testing'
fast5 = '~/DataAnalysis/Marchand/220729_BS_libv2_minION_5000'
ref = '~/xombo/reference_sequences/ref_full_xm.fa' 
out_pre = 'BS_AC_temp.csv' 
out_bc = 'xenomorph_testing/AC_BS_temp.csv' 
model = 'models/ATGCBS_FST001_AG.csv'



wdir = 'xenomorph_testing'
fast5 = '~/DataAnalysis/Marchand/220729_BS_libv2_minION_5000'
ref = '~/xombo/reference_sequences/ref_full_xm.fa' 
out_pre = 'BS_AT_temp.csv' 
out_bc = 'xenomorph_testing/AG_BT_temp.csv' 
model = 'models/ATGCBS_FST001_AT.csv'



wdir = 'xenomorph_testing'
fast5 = '~/DataAnalysis/Marchand/220729_BS_libv2_minION_5000'
ref = '~/xombo/reference_sequences/ref_full_xm.fa' 
out_pre = 'BS_GC_temp.csv' 
out_bc = 'xenomorph_testing/BS_GC_temp.csv' 
model = 'models/ATGCBS_FST001_BgSc.csv'




wdir = 'xenomorph_testing'
fast5 = '~/DataAnalysis/Marchand/220805_ATGC_libv2_minION_5000'
ref = '~/xombo/reference_sequences/ref_full_xm.fa' 
out_pre = 'ATGC_BgSc_BC_Test.csv' 
out_bc = 'xenomorph_testing/ATGC_BgSc_BC_Test.csv' 
model = 'models/ATGCBS_FST001_BgSc.csv'
#model = 'models/ATGCPZ_FLG002_5k.csv'



wdir = 'xenomorph_testing'
fast5 = '~/DataAnalysis/Marchand/220729_BS_libv2_minION_5000'
ref = '~/xombo/reference_sequences/ref_full_xm.fa' 
out_pre = 'BS_GC_temp.csv' 
out_bc = 'xenomorph_testing/BS_GC_temp.csv' 
model = 'models/ATGCBS_FST001_BgSc.csv'





wdir = 'xenomorph_testing'
fast5 = '~/DataAnalysis/Marchand/22801_PZ_libv2_minION_5000'
ref = '~/xombo/reference_sequences/ref_full_pm.fa' 
out_pre = 'PZ_Out_Test.csv' 
out_bc = 'xenomorph_testing/PZ_BC_Test.csv' 
model = 'models/ATGCPZ_FST001.csv'
#model = 'models/ATGCPZ_FLG002_5k.csv'



wdir = '~/xenomorph/xenomorph-beta/220812_temp_HKtest/2203_PZ_tailing'
fast5 = '~/xenomorph/xenomorph-beta/220812_temp_HKtest/2203_PZ_tailing/fast5_multi'
ref = '~/xenomorph/xenomorph-beta/220812_temp_HKtest/2203_PZ_tailing/PZ_HPA_HPB_fasta_s.fa'
out_pre = 'PZ_Out_Test.csv' 
out_bc = '~/xenomorph/xenomorph-beta/220812_temp_HKtest/PZ_BC_Test.csv' 
model = 'models/ATGCPZ_FST001.csv'
#model = 'models/ATGCPZ_FLG002_5k.csv'
'''
'''

wdir = '~/xenomorph/xenomorph-dev/xenomorph/xenomorph_testing'
fast5 = '~/DataAnalysis/Marchand/220912_JV_libv2'
ref = '~/xombo/reference_sequences/ref_full_JV.fa' 
out_pre = '220912_JV_libv2_JV_levels.csv' 
out_bc = '~/xenomorph/xenomorph-dev/xenomorph/xenomorph_testing/JV_FLG001_morphout.csv' 
model = 'ATGCJV'


wdir_2 = '~/xenomorph/xenomorph-dev/xenomorph/xenomorph_testing'
fast5_2 = '~/DataAnalysis/Marchand/220913_XK_libv2'
ref_2 = '~/xombo/reference_sequences/ref_full_XK.fa' 
out_pre_2 = '220913_XK_libv2_XK_levels.csv' 
out_bc_2 = '~/xenomorph/xenomorph-dev/xenomorph/xenomorph_testing/XK_FLG001_morphout.csv' 
model_2 = 'ATGCXK'
'''

#wdir = '~/xenomorph/xenomorph-dev/xenomorph/xenomorph_testing/220917_BSc_FLG001_param'
#fast5 = '~/DataAnalysis/Marchand/220915_BSc_libv2'
#ref = '~/xombo/reference_sequences/ref_full_BS.fa' 
#out_pre = '220915_BSc_libv2_BSc_levels_FLG.csv' 
#out_bc = '~/xenomorph/xenomorph-dev/xenomorph/xenomorph_testing/BSc_FLG001_morphout_FLG.csv' 
#model = 'ATGCBS'


wdir = '~/xenomorph/xenomorph-dev/xenomorph/xenomorph-xemora/xx-test/230610test'
fast5 = '/home/marchandlab/xenomorph/xemora-beta/xemora-test/230309_Xemora_PZ_GC_train/230308_GC_Xemora_train_110_fast5s'
ref = '/home/marchandlab/xenomorph/xemora-beta/ref/G_Xemora_train.fa'
out_pre = '230610_levels.csv' 
out_bc = '~/xenomorph/xenomorph-dev/xenomorph/xenomorph_testing/JV_FLG001_morphout_AT.csv' 
model = 'ATGCJV'
'''
wdir_2 = '~/xenomorph/xenomorph-dev/xenomorph/xenomorph_testing/220916_JVXK_ATsubst'
fast5_2 = '~/DataAnalysis/Marchand/220913_XK_libv2'
ref_2 = '~/xombo/reference_sequences/ref_full_XK.fa' 
out_pre_2 = '220913_XK_libv2_XK_levels_AT.csv' 
out_bc_2 = '~/xenomorph/xenomorph-dev/xenomorph/xenomorph_testing/XK_FLG001_morphout_AT.csv' 
model_2 = 'ATGCXK'
'''


#Preprocess
if True: 
    cmd = 'python xenomorph.py preprocess -w '+wdir+' -f '+fast5+' -r '+ref+' -b -x -o  '+out_pre
#    cmd = 'python xenomorph.py preprocess -w '+wdir+' -f '+fast5+' -r '+ref+' -x -o  '+out_pre
#cmd = 'python xenomorph.py preprocess -w '+wdir+' -f '+fast5+' -r '+ref+' -o '+out_pre
    os.system(cmd) 

    print("1st preprocessing finished. \n")

#print("Starting preprocess for XK. \n")
if False: 
    cmd = 'python xenomorph.py preprocess -w '+wdir_2+' -f '+fast5_2+' -r '+ref_2+' -b -x -o  '+out_pre_2
#    cmd = 'python xenomorph.py preprocess -w '+wdir+' -f '+fast5+' -r '+ref+' -x -o  '+out_pre
#cmd = 'python xenomorph.py preprocess -w '+wdir+' -f '+fast5+' -r '+ref+' -o '+out_pre
    os.system(cmd) 

    print("2nd preprocessing finished. \n")
#Morph 

if False: 
#cmd = 'python xenomorph.py morph -l '+wdir+'/'+out_pre+ ' -m '+model+' -o '+out_bc
    cmd = 'python xenomorph.py morph -l '+wdir+'/'+out_pre+ ' -m '+model+' -o '+out_bc
    os.system(cmd)

#Qucik stats
if False:

    cmd = 'python xenomorph.py stats -i '+out_bc
    os.system(cmd)

#Global stats
if False:

    cmd = 'python xenomorph.py stats -i '+out_bc 
    os.system(cmd)


#python xombo.py subset ~/DataAnalysis/Marchand/220801_PZ_libv2_minION_single ~/xenomorph/xenomorph-dev/xenomorph/xenomorph_testing/22801_PZ_libv2_minION_5000_single 5000
#python xombo.py subset ~/DataAnalysis/Marchand/220805_ATGC_libv2_minION_single ~/xenomorph/xenomorph-dev/xenomorph/xenomorph_testing/220805_ATGC_libv2_minION_5000_single 5000
--------------------------------------------




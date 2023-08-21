import os

################################################
#Workflow integration for Xenomorph-Xemora
wdir = '/Xenomorph_Testing'
fast5 = '/path_to_fast5_directory/fast5'
ref = '/path_to_FASTA_reference_with_XNA_bases_in_sequence/xna_reference.fa'
out_file_prefix = '/output_prefix'#Output file parameters

#Level generation
run_preprocess = True

#Basecalling
run_morph = True
run_stats = True
run_global_morph = True

#Generate model from a level file
run_model_gen = True

#Morph model
model = 'ATGCPZ'

#File prefixes
out_pre = out_file_prefix+'_levels.csv'
out_bcpr = out_file_prefix+'_per-read_basecalls.csv'
out_bcglobal = out_file_prefix+'_global_basecalls.csv'
################################################


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







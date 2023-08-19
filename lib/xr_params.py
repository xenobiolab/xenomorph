########################################################################
########################################################################
"""
xr_params.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. A. Marchand

Updated: 8/16/23 
"""

#######################################################################

import numpy as np
from remora import io, refine_signal_map, util


#############################REMORA PARAMETERS####################

#Kmer level table file path 
level_table = '/home/marchandlab/Dev/xenomorph-xemora/kmers/4mer_9.4.1.csv'

#Signal extraction type (norm, dac, pa; default = 'norm')
signal_type = 'norm'

#Set min match score. Sequences should have a match quality greather than match_score to pass filter. (default = 0)
min_match_score=0

#Skip indexing bam file for getting quality scores and mapping scores. Will default to QS 15, Map score of 10
skip_qscore_extract  = True

#Re-basecall pod5 file. Required if new reference files are being used. 
basecall_pod = False

#Re-generate BAM files for reference-based basecalling.
regenerate_bam = False

#Re-generate training or basecalling chunks.
regenerate_chunks = False

#generate a .bai file for your bam file -- you need to do this the first time you run analysis
gen_bai = False

#Override XNA position detection and specify position. Required for ATGC-only signal extraction or if reference file does not specify XNA position. 
force_extract_position = False

#If force_extract_position is true, specify position in reference. 
extract_pos = 67

#Preprocess a maximum number of reads (default = 0 == all reads) 
max_num_reads = 1000


######################LEVEL RESCALING######################

#If true, use manual rescale. If false, automatically calculate rescaling parameters. 
manual_rescale_override = True

#Manually set kmer rescale slope
manual_rescale = 1.4826

#Manually set kmer shift
manual_reshift = 0

#Global rescale is automatically updated after auto-scaling is performed
global_rescale = 1.4541330366934333

#Global reshift is automatically updated after auto-scaling is performed
global_reshift = 0.028562533484942365

#Path to model used as ground truth for rescaling 
rescale_reference_model_path = 'models/libv2/ATGC_libv2_FLG001.csv'

#Rescale metric (mean or median. Default = median)
rescale_metric = 'median'

#Rescale method calculation (Thiel-Sen, polyfit. Default = 'Thiel-Sen')
rescale_method ='Thiel-Sen'

#Rescale max number of reads to use
rescale_max_num_reads = 200

#Number of levels before and after to extract surrounding an XNA (default = 3) 
rescale_xmer_boundary = 30

#Number of bases before and after XNA that are required in a matching read (default = 30 alt) 
rescale_xmer_padding = 35

#Show rescale plot 
rescale_save_plot = True

######################REMORA SIGNAL REFINER######################
#Set up SigMapRefiner
sig_map_refiner = refine_signal_map.SigMapRefiner(
    kmer_model_filename=level_table,
    scale_iters= 20,
    do_rough_rescale = True, 
    do_fix_guage= True
)
################################################################


######################XFASTA GENERATION######################

#Fasta2x - write sequences to xfasta even if they do no contain XNAs. Default = False 
write_no_xna_seq = True

#Fasta2x - Write gaps in place of XNAs in fasta reference file for null testing
write_gaps = False


############################################################
##Model Training and Basecalling Parameters

#Modified base in Fasta sequence you wish to train model or use model to basecall
mod_base = 'J'
mod_rev_base = 'F'

#Most similar substituted canonical base you will be comparing against 
#can_base = 'G'
#can_rev_base = 'C'


#Extent of Kmer content (-,+) to store for model training
kmer_context ='4 4' 

#Extent of chunk context (centered around modified base) 
chunk_context = '50 50' 

#Proportion of reads to use for validation 
val_proportion = 0.01

############################################################
#Guppy Base caller configuration

#Path to guppy basecaller
basecaller_path ='~/ont-guppy/bin/guppy_basecaller' 

#GPU enabled 
device_type = 'cuda:all' 

#Config file 
guppy_config_file = 'dna_r9.4.1_450bps_hac.cfg' 



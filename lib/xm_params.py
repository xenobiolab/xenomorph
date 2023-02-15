########################################################################
########################################################################
"""
xm_params.py

Description: Here you can modify parameters used for processing, extraction,
alignment, and basecalling. Default values 

Title: Synthesis and Sequencing of 12-Letter Supernumerary DNA

By: H. Kawabe, C. Thomas, A. Laszlo, S. Hoshika, L. Miessner, J. M. Craig, 
J. Gundlach, Myong-Jung Kim, Myong-Sang Kim, S. A. Benner, J. A. Marchand

Updated: 1/15/23
"""
########################################################################
########################################################################




import numpy as np
############################################################
######################PARAMETER############################
############################################################

#Standard basepairs written in 'purine pyrimidine' order
standard_base_pairs = ['AT','GC']

#Convert this to set
standard_bases = np.concatenate(list(list(i) for i in standard_base_pairs))

#Alternative basepairs written in 'purine pyrimidine' order
xna_base_pairs = ['BS','PZ','JV','XK','QW', 'ER']

#Seperate to segment data as non-complementary pairs 
xna_segmentation_model_sets = ['B','S','PZ','JV','X','K']

#Most similar canonical pair for each XNA
confounding_pairs =  ['BA','SA','PG','ZC','JC','VG','XA','KG'] # Verified optimal pairs

#Alternative basepairs written in 'purine pyrimidine' order
xna_bases = np.concatenate(list(list(i) for i in xna_base_pairs))

######################XFASTA GENERATION######################

#Fasta2x - write sequences to xfasta even if they do no contain XNAs. Default = False 
write_no_xna_seq = True

#Fasta2x - Write gaps in place of XNAs in fasta reference file for null testing
write_gaps = False


######################LEVEL EXTRACTION#######################

#Bases after XNA to use for kmer model (5'->3') (default = 2) 
kmer_plus = 2

#Bases before XNA to use for kmer model (5'->3') (default = 1) 
kmer_minus = 1

#Number of levels before and after to extract surrounding an XNA (default = 3) 
xmer_boundary = 3

#Number of bases before and after XNA that are required in a matching read (default = 30 alt) 
#xmer_padding = 15
xmer_padding = 15


######################RESQUIGGLE PARAMS#######################

#Default: 4.2 4.2 300 1500 20.0 40 750 2500 250
signal_align_params = '4.2 4.2 300 1500 20.0 40 750 2500 250'


#Segmentation parameters DNA default 5 3 1 5
segmentation_params = '5 3 1 5' #5 3 1 5

#Currently manually set
skip_sequence_rescaling = False 


######################Quick stat/stats params#######################
#Quick stat params 
qscore_filter = 15#(default 15) 

#Signal filter
signal_filter = 1.5#(default 1.1)

#Minimum number of reads required for calculating a concensus value
concensus_stat_filter = 10


############################################################
######################CONFIGURATION#########################
############################################################

###PREPROCESSING 
#Guppy Base caller - Note: Needs to be manually configured. 
basecaller_path ='~/ont-guppy/bin/guppy_basecaller' 

#The following parameters are not changed. HAC basecalling is fixed. 
flowcell_type = 'FLO-FLG001' #FLO-FLG001 #FLO-MIN106
kit_type = 'SQK-LSK110' 
device_type = 'cuda:all' 
config_file = 'dna_r9.4.1_450bps_hac.cfg' 

#Use reference file in basecalling alignment
use_reference_to_basecall=False

#Assign Fastq from 'pass', 'fail', or 'both' basecall files 
read_assign = 'both'

####MORPH PARAMETERS 
#Nmer model to use (nnnxnnn, nxnn, nnxnn etc) 
nmer_model = 'nnnxnnn'

#Used to specify where the zero position of the kmer is relative to the kmer boundaries (x = 0) 
kmer_mask = ['nxnn']

#Optional kmer weights to use (default unweighted: [ 1, 1, 1 ,1] 
kmer_weights = [1, 1, 1 ,1] 

#Kmer models to use (e.g. [4, 5, 6] will use 4kmer, 5kmer, and 6kmer models in ML calculations 
kmer_sizes = [len(kmer_mask)] 

#Set mu as either: 'KDE Mean level', 'Median level',or 'Mean level'. Basecalling performance can vary depending on accuracy of measurement. 
mu = 'Mean level'

#Set sigma as 'Global-Median', 'Global-Mean' for global mean (recommended),  'Kmer' for kmer-specific, or set as float for fixed
sigma = 0.4

#Likelihood ratio (Normal; Outlier-Robust; Outlier-Robust)
likelihood_ratio = 'Outlier-Robust'

#Threshold for rejecting or accepting LLR call in global_morph (default = 0)
likelihood_threshold = 0

#Max reads for morph
max_reads = 0

#Max reads for null morph
max_reads_null = 0

#####Outlier robust parameters#####
#Changes position and magnitude of max/min of outlier-robust function
Sf =  4#4

#Changes magnitude of max/min of outlier-robust function but not position 
Sf2 = 3#3 

#Changes magnitude of max/min of outlier-robust function but not position 
Sp = 0.3#0.3 

#Bases for performing alternative hypothesis: Options: all [ATGCXY], standard [ATGCX], pyrpur [AGX/TCY], confounding[GX/CY]
alt_base_type = 'all'




########################################################################
########################################################################
"""
xm_params.py

Description: Here you can modify parameters used for processing, extraction,
alignment, and basecalling. Additional parameter file is available for 
remora API usage (lib/xr_params.py). 

Title: Synthesis and Sequencing of 12-Letter Supernumerary DNA

By: H. Kawabe, C. Thomas, A. Laszlo, S. Hoshika, L. Miessner, J. M. Craig, 
J. Gundlach, Myong-Jung Kim, Myong-Sang Kim, S. A. Benner, J. A. Marchand

Updated: 8/20/23
"""
########################################################################
########################################################################
import numpy as np


############## SEGMENTATION MODE ##############
#Segmentation mode - either use Tombo (from Xenomorph v1.0) or Remora (Implemented of v1.5). Refer to documentation for more information about which is more appropriate for your use.
segmentation_mode = 'tombo'


############## ESTABLISH DNA AND XNA RULES ##############
#Standard basepairs written in 'purine pyrimidine' order
standard_base_pairs = ['AT','GC']

#Convert this to set
standard_bases = np.concatenate(list(list(i) for i in standard_base_pairs))

#Alternative basepairs written in 'purine pyrimidine' order
xna_base_pairs = ['BS','PZ','JV','XK','QW', 'ER', 'NN']

#Seperate to segment data as non-complementary pairs 
xna_segmentation_model_sets = ['B','S','PZ','JV','X','K','N']

#Most similar canonical pair for each XNA
confounding_pairs =  ['BA','SA','PG','ZC','JC','VG','XA','KG','NN'] # Verified optimal pairs

#Alternative basepairs written in 'purine pyrimidine' order
xna_bases = np.concatenate(list(list(i) for i in xna_base_pairs))


############## PREPROCESSING PARAMETERS ##############
####### PREPROCESSING - XFASTA FILE GENERATION
#Regenerate xfasta file(s).
regenerate_xfasta = True

#Fasta2x - write sequences to xfasta even if they do no contain XNAs. Default = False 
write_no_xna_seq = True

#Fasta2x - Write gaps in place of XNAs in fasta reference file for null testing
write_gaps = False


####### PREPROCESSING - GUPPY BASECALLING
#Guppy Base caller - Note: Needs to be manually configured. 
basecaller_path ='~/ont-guppy/bin/guppy_basecaller' 

#Mininum qscore set on guppy basecalling for pass/fail (default = 9) 
guppy_min_qscore = 9 

#GPU device (default = cuda:all)
device_type = 'cuda:all' 

#Guppy config file (default = dna_r9.4.1_450bps_hac.cfg)
guppy_config_file = 'dna_r9.4.1_450bps_hac.cfg' 

#Guppy align type flag: full, coarse, or auto (default = full)
guppy_align_type = 'full'

#Use reference file in basecalling alignment (default = true)
use_reference_in_basecall=True

#Use basecalled reads from 'pass', 'fail', or 'both' basecall files for segmentation (default = both)
read_assign = 'both'



####### PREPROCESSING - TOMBO SEGMENTATION
#Default: 4.2 4.2 300 1500 20.0 40 750 2500 250
signal_align_params = '4.2 4.2 300 1500 20.0 40 750 2500 250'

#Segmentation parameters DNA default 5 3 1 5
segmentation_params = '5 3 1 5' #5 3 1 5

#Currently manually set
skip_sequence_rescaling = False 


####### PREPROCESSING - LEVEL EXTRACTION
#Bases after XNA to use for kmer model (5'->3') (default = 2) 
kmer_plus = 2

#Bases before XNA to use for kmer model (5'->3') (default = 1) 
kmer_minus = 1

#Number of levels before and after to extract surrounding an XNA (default = 3) 
xmer_boundary = 3

#Number of bases before and after XNA that are required in a matching read (default = 30 alt) 
xmer_padding = 15



############## MORPH PARAMETERS ##############
####### READ FILTERING AND STAT CALCULATIONS
#Quick stat params 
qscore_filter = 9#(default = 9) 

#Signal filter
signal_filter = 1.5 #1.5#(default =1.1 for Tombo; default = 100 for remora)

#Minimum number of reads required for calculating a concensus value (both per-read and global)
concensus_stat_filter = 20


####MORPH PARAMETERS 
#Nmer model to use (nnnxnnn, nxnn, nnxnn etc) 
nmer_model = 'nnnxnnn'

#Used to specify where the zero position of the kmer is relative to the kmer boundaries (x = 0) 
kmer_mask = ['nxnn']

#Optional kmer weights to use (default unweighted: [1, 1, 1 ,1])
kmer_weights = [1, 1, 1 ,1] 

#Kmer models to use (e.g. [4, 5, 6] will use 4kmer, 5kmer, and 6kmer models in ML calculations 
kmer_sizes = [len(kmer_mask)] 

#Set mu (per-read) as either: 'KDE Mean level', 'Median level',or 'Mean level'. Basecalling performance can vary depending on accuracy of measurement. 
mu = 'Median level'

#Set global mu (per-sequence) as either: 'mean or median". Used for concensus recall calculations at the sequence level. 
mu_global = 'Median'

#Set sigma as 'Global-Median', 'Global-Mean' (recommended),'Kmer' (for kmer-specific), or as float for fixed
sigma = 'Global-Mean'

#Max reads to basecall using morph command (default = 0 == all reads)
max_reads = 20000

#Max reads for null morph (beta analysis feature only)
max_reads_null = 0


####### LIKELIHOOD RATIO CALCULATIONS
#Likelihood ratio calculation (Normal or Outlier-Robust; default = Outlier-Robust)
likelihood_ratio = 'Outlier-Robust'

#Threshold for rejecting or accepting LLR call in global_morph (default = 0)
likelihood_threshold = 0

#Changes position and magnitude of max/min of outlier-robust function (default = 4). For ORLLR only. 
Sf =  4

#Changes magnitude of max/min of outlier-robust function but not position (default = 3). For ORLLR only. 
Sf2 = 3

#Changes magnitude of max/min of outlier-robust function but not position (default = 0.3). For ORLLR only. 
Sp = 0.3 

#Bases for performing alternative hypothesis: Options: all [ATGCXY], standard [ATGCX], pyrpur [AGX/TCY], confounding[GX/CY]
alt_base_type = 'confounding'




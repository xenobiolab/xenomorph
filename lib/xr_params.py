########################################################################
########################################################################
"""
xr_params.py 
Param file for remora segmentation support to Xenomorph. For all other
parameters required for Xenomorph, see lib/xm_params.py


Updated: 9/4/23 
"""

#######################################################################

import sys
import numpy as np

############## PREPROCESSING PARAMETERS ##############
####### PREPROCESSING - REMORA SEGMENTATION (Advanced)

#Signal extraction type (norm, dac, pa; default = 'norm')
signal_type = 'norm'

#Set min match score. Sequences should have a match quality greather than match_score to pass filter. (default = 0)
min_match_score=0

#Skip indexing bam file for getting quality scores and mapping scores. Will default to QS 15, Map score of 10
skip_qscore_extract  = True

#Override XNA position detection and specify position. Required for ATGC-only signal extraction or if reference file does not specify XNA position. 
force_extract_position = False

#If force_extract_position is true, specify position in reference. Advanced use only. 
extract_pos = 67

#Preprocess a maximum number of reads (default = 0 == all reads) 
max_num_reads = 100000




####### PREPROCESSING - REMORA RESCALING (Advanced)
#If true, use manual rescale. If false, automatically calculate rescaling and reshift parameters. 
manual_rescale_override = True

#Manually set kmer rescale slope (this is pa to norm scaling factor)
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

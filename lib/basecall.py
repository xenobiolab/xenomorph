########################################################################
########################################################################
"""
basecall.py

Description: This script is part of the Xombo file handling package. 
Calls guppy for initial fast5 file processing. Default model is set 
to r.9.4.1 DNA HAC. Path to guppy needs to be specified in xm_params.py. 
Basecalling is required for initial fast5-sequence assignment prior to 
signal-to-sequence alignment according to the Tombo/Remora pipeline. 

Title: Synthesis and Sequencing of 12-Letter Supernumerary DNA

By: H. Kawabe, C. Thomas, S. Hoshika, Myong-Jung Kim, Myong-Sang Kim, L. Miessner, J. M. Craig, 
J. Gundlach, A. Laszlo,  S. A. Benner, J. A. Marchand

Updated: 8/19/23
"""
########################################################################
########################################################################

import os 
import sys 
from xm_params import *

#Inputs 
input_folder = os.path.normpath(sys.argv[1])
output_folder = os.path.normpath(sys.argv[2])
reference_file = os.path.normpath(sys.argv[3])

# If output basecall folder doesn't exist, then create it.
CHECK_FOLDER = os.path.isdir(output_folder)
if not CHECK_FOLDER:
    os.makedirs(output_folder)
    print("Xenomorph Status - [Preprocess] Creating fastq output directory : "+output_folder)


# Run guppy basecall with parameters set in xm_params.py
try: 
    if segmentation_mode.lower() == 'tombo': 
        if use_reference_in_basecall==True: 
            print("Xenomorph Status - [Preprocess] Performing guppy basecalling using reference xfasta file: "+reference_file+".")
            cmd = basecaller_path+' -i '+input_folder+' -s '+ output_folder + '  --device '+device_type+' -a '+reference_file+' -c '+guppy_config_file+' --min_qscore '+str(guppy_min_qscore)+' --align_type '+guppy_align_type
            
        else: 
            print("Xenomorph Status - [Preprocess] Performing guppy basecalling without reference alignment.")
            cmd = basecaller_path+' -i '+input_folder+' -s '+ output_folder + ' --device '+device_type +' -c '+guppy_config_file+' --min_qscore '+str(guppy_min_qscore)
            
    if segmentation_mode.lower() == 'remora': 
        print("Xenomorph Status - [Preprocess] Performing guppy basecalling using reference xfasta file for segmentation: "+reference_file+".")
        cmd=(os.path.expanduser(basecaller_path)+' -i '+input_folder+' -s '+output_folder+' -c '+guppy_config_file+' --min_qscore '+str(guppy_min_qscore)+
        ' --device '+device_type +' --align_type '+guppy_align_type+' --bam_out --index --moves_out -a '+reference_file)

    #Execute guppy to basecall
    os.system(cmd) 
    
except: 
    print('Xenomorph Status - [Error] Basecalling with guppy failed. Ensure guppy is properly installed and configured.') 
    sys.exit(1)




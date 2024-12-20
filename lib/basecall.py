########################################################################
########################################################################
"""
basecall.py

Description: This script is part of the Xombo file handling package. 
Calls guppy for initial fast5 file processing. Default model is set 
to r.9.4.1 DNA HAC. Path to guppy needs to be specified in xm_params.py. 
Basecalling is required for initial fast5-sequence assignment prior to 
signal-to-sequence alignment according to the Tombo/Remora pipeline. 

Cite us: 
H. Kawabe, C. Thomas, S. Hoshika, Myong-Jung Kim, Myong-Sang Kim, L. Miessner, N. Kaplan, J. M. Craig, 
J. Gundlach, A. Laszlo,  S. A. Benner, J. A. Marchand. "Enzymatic Synthesis and Nanopore Sequencing of a 12-Letter Supernumerary DNA" 
Nature Communications. 14. (2023). DOI: 10.1038/s41467-023-42406-z 


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


#Set guppy_config_file if empty or auto
if guppy_config_file == '' or guppy_config_file =='auto': 
	if flowcell_version == '9.4.1': 
		guppy_config_file = 'dna_r9.4.1_450bps_hac.cfg'
	if flowcell_version == '10.4.1': 
		guppy_config_file = 'dna_r10.4.1_e8.2_400bps_sup.cfg' 


# Run guppy basecall with parameters set in xm_params.py
try: 

    if use_gpu == True: 
        device_configuration = ' --device '+device_type 
    else: 
        device_configuration = ' --cpu_threads_per_caller '+str(cpu_threads)
    
    if segmentation_mode.lower() == 'tombo': 
        if use_reference_in_basecall==True: 
            print("Xenomorph Status - [Preprocess] Performing guppy basecalling using reference xfasta file: "+reference_file+".")
            cmd = basecaller_path+' -i '+input_folder+' -s '+ output_folder + device_configuration +' -a '+reference_file+' -c '+guppy_config_file+' --min_qscore '+str(guppy_min_qscore)+' --align_type '+guppy_align_type
            
        else: 
            print("Xenomorph Status - [Preprocess] Performing guppy basecalling without reference alignment.")
            cmd = basecaller_path+' -i '+input_folder+' -s '+ output_folder + device_configuration +' -c '+guppy_config_file+' --min_qscore '+str(guppy_min_qscore)
            
    if segmentation_mode.lower() == 'remora': 
        print("Xenomorph Status - [Preprocess] Performing guppy basecalling using reference xfasta file for segmentation: "+reference_file+".")
        cmd=(os.path.expanduser(basecaller_path)+' -i '+input_folder+' -s '+output_folder+' -c '+guppy_config_file+' --min_qscore '+str(guppy_min_qscore)+
        device_configuration+' --align_type '+guppy_align_type+' --bam_out --index --moves_out -a '+reference_file)

    #Execute guppy to basecall
    os.system(cmd) 
    
except: 
    print('Xenomorph Status - [Error] Basecalling with guppy failed. Ensure guppy is properly installed and configured.') 
    sys.exit(1)




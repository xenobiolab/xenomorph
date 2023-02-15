########################################################################
########################################################################
"""
basecall.py

Description: This script is part of the Xombo file handling package. 
Calls guppy for initial fast5 file processing. Default model is set 
to r.9.4.1 DNA HAC. Path to guppy needs to be specified in xm_params.py. 
Basecalling is required for initial fast5-sequence assignment prior to 
signal-to-sequence alignment according to the Tombo pipeline. 

Title: Synthesis and Sequencing of 12-Letter Supernumerary DNA

By: H. Kawabe, C. Thomas, S. Hoshika, Myong-Jung Kim, Myong-Sang Kim, L. Miessner, J. M. Craig, 
J. Gundlach, A. Laszlo,  S. A. Benner, J. A. Marchand

Updated: 2/14/23
"""
########################################################################
########################################################################
import os 
import sys 
from xm_params import *

#Inputs 
input_folder = sys.argv[1].replace('//','/')
output_folder = sys.argv[2].replace('//','/')



#Parameters
# If folder doesn't exist, then create it.
CHECK_FOLDER = os.path.isdir(output_folder)
if not CHECK_FOLDER:
    os.makedirs(output_folder)
    print("Xenomorph Status - [Preprocess] Creating fastq output directory : "+output_folder)




try: 
	if use_reference_to_basecall==True: 
		reference_file = sys.argv[3].replace('//','/')
		print("Xenomorph Status - [Preprocess] Basecalling with reference xfasta file "+sys.argv[3]+".")
		exec = basecaller_path+' -i '+input_folder+' -s '+ output_folder + '  --device '+device_type+' -a '+sys.argv[3]+' -c dna_r9.4.1_450bps_hac.cfg'
		os.system(exec) 
	else: 

		print("Xenomorph Status - [Preprocess] De novo basecalling samples.")
		exec = basecaller_path+' -i '+input_folder+' -s '+ output_folder + '  --device '+device_type +' -c dna_r9.4.1_450bps_hac.cfg'
		os.system(exec) 


except: 
	print('Error: Basecalling with guppy failed. Ensure guppy is properly installed and configured.') 
	exit()




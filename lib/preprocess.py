########################################################################
########################################################################
"""
preprocess.py

Description: This script is part of the Xombo file handling package. Used
to preprocess reads for Tombo analysis. Assigns fastq basecalls to fast5 
single read files. Automatically called in the Xenomorph pipeline. 
For more information: python xombo.py preprocess -h

Cite us: 

H. Kawabe, C. Thomas, S. Hoshika, Myong-Jung Kim, Myong-Sang Kim, L. Miessner, N. Kaplan, J. M. Craig, 
J. Gundlach, A. Laszlo,  S. A. Benner, J. A. Marchand. "Enzymatic Synthesis and Nanopore Sequencing of a 12-Letter Supernumerary DNA" 
Nature Communications. 14. (2023). DOI: 10.1038/s41467-023-42406-z 

Updated: 2/14/23
"""
########################################################################
########################################################################

import os 
import sys 
from xm_params import * 


#Inputs 
input_folder = sys.argv[1]+'/'.replace('//','/')
fastq_folder = sys.argv[2]+'/'.replace('//','/')


# If folder doesn't exist, then create it.
if read_assign == 'pass' or read_assign=='both': 
    passfail = 'pass' 
    print("Xenomorph Status - [Preprocess] Using tombo preprocess to assign /pass fastq basecalls to fast5 files.")
    try: 

	    exec = 'tombo preprocess annotate_raw_with_fastqs --fast5-basedir '+input_folder+' --fastq-filenames '+fastq_folder+passfail+'/*.fastq --sequencing-summary-filenames '+fastq_folder+'sequencing_summary.txt --overwrite'
	    os.system(exec) 


    except: 
	    print('Error: Failed to run [tombo preprocess]. Ensure Tombo is installed correctly.')
	    exit()



if read_assign == 'fail' or read_assign=='both': 
    passfail = 'fail' 
    print("Xenomorph Status - [Preprocess] Using tombo preprocess to assign /fail fastq basecalls to fast5 files.")
    fail_folder = fastq_folder+passfail
    if os.path.exists(fail_folder)==True: 
        try: 
	        exec = 'tombo preprocess annotate_raw_with_fastqs --fast5-basedir '+input_folder+' --fastq-filenames '+fastq_folder+passfail+'/*.fastq --sequencing-summary-filenames '+fastq_folder+'sequencing_summary.txt --overwrite'
	        os.system(exec) 

        except: 
	        print('Error: Failed to run [tombo preprocess]. Ensure Tombo is installed correctly.') 




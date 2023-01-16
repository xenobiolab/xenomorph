########################################################################
########################################################################
"""
split.py

Description: This script is part of the Xombo file handling package. Used 
to split multi-fast5 into single-fast5 directory. Required for tombo handling. 
For more information: python xombo.py split -h


Title: XNA tailing enables nanopore sequencing of a 12-letter genetic code

By: H. Kawabe, C. Thomas, A. Laszlo, S. Hoshika, L. Miessner, J. M. Craig, 
J. Gundlach, Myong-Jung Kim, Myong-Sang Kim, S. A. Benner, J. A. Marchand

Updated: 1/15/23
"""
########################################################################
########################################################################


#Split
import pandas as pd 
import os 
import sys 

input_folder = sys.argv[1].replace('//','/')
output_folder = sys.argv[2].replace('//','/')



print('Splitting mutli-fast5 files into single-fast5 files') 


# If folder doesn't exist, then create it.
CHECK_FOLDER = os.path.isdir(output_folder)
if not CHECK_FOLDER:
    os.makedirs(output_folder)
    print("Creating output folder : ", output_folder)





#Run multi to single fast5 
try: 
	exec = 'multi_to_single_fast5 --input_path '+input_folder+' --save_path '+output_folder+'/temp --recursive'
	os.system(exec) 


except: 
	print('Error: Split requires multi_to_single_fast5 installation. Ensure it is properly installed.') 






print('Merging split fast5 files into a single directory ') 


#Move files
exec = 'find '+output_folder+'/temp -maxdepth 2 -type f -print0 | xargs -0 mv -t '+output_folder
os.system(exec) 

#Delete temp folder
exec = 'rm -r '+output_folder+'/temp' 
os.system(exec)





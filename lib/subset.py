########################################################################
########################################################################
"""
subset.py

Description: This script is part of the Xombo file handling package. Used 
to subset a larger single fast5 directory into a user-specified smaller 
directory for testing, prototyping, and troubleshooting. 
For more information: python xombo.py subset -h


Title: XNA tailing enables nanopore sequencing of a 12-letter genetic code

By: H. Kawabe, C. Thomas, A. Laszlo, S. Hoshika, L. Miessner, J. M. Craig, 
J. Gundlach, Myong-Jung Kim, Myong-Sang Kim, S. A. Benner, J. A. Marchand

Updated: 1/15/23
"""
########################################################################
########################################################################



#Split
import os 
import sys 

#Inputs 
input_folder = sys.argv[1].replace('//','/')
output_folder = sys.argv[2].replace('//','/')
n_files = str(sys.argv[3])

mode = 'copy' #Defualt copy, but can be split

# If folder doesn't exist, then create it.
CHECK_FOLDER = os.path.isdir(output_folder)
if not CHECK_FOLDER:
	os.makedirs(output_folder)
	print("Creating output folder : ", output_folder)


if mode == 'copy': 
	print('Subsetting fast5 files - copying '+str(n_files)+ ' fast5 files to new directory: '+output_folder)
	try: 
	    exec = 'find '+input_folder+' -maxdepth 1 | head -'+n_files+'  |  xargs cp -p -t '+output_folder
	    os.system(exec) 
	    print('Successfully moved files. Safe to ignore pipe error.')


	except: 
	    print('Error: Moving folders failed. Ensure paths are correct and that proper permissions have been set.') 


if mode == 'split': 
	try: 
	    print('Warning: Subset is in move mode, not copy mode')
	    exec = 'find '+input_folder+' -name '+"'*.fast5'"+' -maxdepth 1 | head -'+n_files+'  |  xargs mv -t '+output_folder
	    os.system(exec) 
	    print(exec)
	    print('Successfully moved files. Safe to ignore pipe error.')


	except: 
	    print('Error: Moving folders failed. Ensure paths are correct and that proper permissions have been set.') 





########################################################################
########################################################################
"""
xr_tools.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. A. Marchand

Updated: 3/2/23
"""
########################################################################
########################################################################

import pandas as pd
import os
import glob
import pysam
from pathlib import Path




def filter_primary_alignments_file(input_bam_file, output_bam_file):
    # Open the BAM file
    samfile = pysam.AlignmentFile(input_bam_file, "rb")

    # Create a new BAM file for primary alignments
    primary_alignments_file = pysam.AlignmentFile(output_bam_file, "wb", header=samfile.header)

    # Iterate through each read in the BAM file
    for read in samfile:
        # Check if the alignment is primary (not secondary)
        if not read.is_secondary and not read.is_supplementary:
            # If it is primary, write it to our new file
            primary_alignments_file.write(read)

    # Close the BAM files
    samfile.close()
    primary_alignments_file.close()

    # Sort the output BAM file
    pysam.sort("-o", output_bam_file, output_bam_file)

    # Create an index for the output BAM file
    pysam.index(output_bam_file)

    # Re-open the new BAM file in read mode and return it
    return pysam.AlignmentFile(output_bam_file, "rb")




def phred_to_qscore(phred_string):
    total_quality = 0
    try: 
        for char in phred_string:
            total_quality += ord(char) - 33
            avg_qscore = total_quality / len(phred_string)
    except: 
        avg_qscore = 0 
    return avg_qscore


def get_read_scores(bam, read_id):
    for read in bam.fetch():
        if read.query_name == read_id:
            mapping_quality_score = read.mapping_quality
            qscore_quality_score = read.qual
            return mapping_quality_score, phred_to_qscore(qscore_quality_score)
    raise ValueError(f"No read with ID {read_id} found in {bam_file}")

def index_read_scores(bam, score_df):
    for read in bam.fetch():
        read_id = read.query_name 
        mapping_quality_score = read.mapping_quality
        qscore_quality_score = phred_to_qscore(read.qual)
        scores = {'Read ID':read_id, "match_score": mapping_quality_score , "q-score": qscore_quality_score}
        score_df.loc[len(score_df)] = scores
    return score_df


def fetch_xna_pos(xm_header):
    pos=xm_header[xm_header.find('XPOS[')+5:-1].split('-')
    xpos = [x.split(':') for x in pos]
    return xpos


def xna_base_rc(xna_base, xna_bp): 
	for x in xna_bp: 
		if xna_base in x:
			if len(x)>1:
			    xx=list(x)
			    xx.remove(xna_base)
			    xrc=xx[0]
			else:
			   xrc=False
	return xrc


def check_xfasta_format(xfasta_file,standard_bases): 
    xfasta_header=False
    xna_in_sequence = False 
    with open(xfasta_file, "r") as fh:
        for line in fh:
            #Get header
            if line[0]=='>':
                header = line
                if 'XPOS[' in header: 
                    xfasta_header=True
                
            #Get sequence
            if line[0]!='>':
            
                #Make all upper case
                uline = line.upper()
                
                #Look for non-standard base
                diff = list(set(uline.replace('\n',''))-(set(standard_bases)))
                
                #This sequence contains an XNA
                if len(diff)>0:
                    xna_in_sequence = True 
    if xfasta_header==True and xna_in_sequence == False: 
        return True 
    else: 

        return False 



#Scans directory and subdirectory to get proper fast5 file path. Not explicitly handled with pod5 commands
def get_fast5_subdir(fast5_dir): 
    path = os.path.normpath(os.path.expanduser(fast5_dir))
    if os.path.exists(path):
        fast5_files = list(Path(path).rglob("*.fast5" ))
        if len(fast5_files)>0:
            fast5_subdir = os.path.dirname(fast5_files[0])
            print('Xenomorph - [Status] Found '+str(len(fast5_files))+' fast5 files in '+fast5_dir)
            return fast5_subdir
        else: 
            print('Xenomorph - [Status] Could not find Fast5 files in specified directory. Check .fast5 exist.')
            return False
    else: 
        print('Xenomorph - [Status] Could not find Fast5 directory. Check path')
        return False


#Check if working directory exists, if not create it. 
def check_make_dir(directory):
    directory = os.path.expanduser(directory)
    if not os.path.isdir(directory):
        os.makedirs(directory)
        print('Xenomorph - [Status] Required directory not found. Creating directory: '+directory)
    return directory


#Fast5 to pod5 conversion
def cod5_to_fast5(fast5_input, pod5_output):
    #for update remora #cmd = 'pod5 convert fast5 '+fast5_input+'/*.fast5 -o '+pod5_output
    cmd = 'pod5 convert fast5 '+fast5_input+'/*.fast5 '+ '-o '+pod5_output
    os.system(cmd)

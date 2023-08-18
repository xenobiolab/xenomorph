########################################################################
########################################################################
"""
xr_get_levels_pos.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. A. Marchand

Updated: 7/2/23 

*****Change log 230702*****
#Created xr_get_levels_pos.py to get levels in a position-specified format. 
#Removed redundant BAM processing code 
#Fixed bug in read count assignment and changed var to num_reads instead of num_file 
#Added optional qscore extraction to properly fall back if needed
#Annotated sig map refiner 
#Added code comments
#Removed j variable as a counter, calculate pass/fail from output dataframes instead
#Added read_map_start and read_map_end info to the level output file 
#Modularized to handle more than one XNA per sequence
#Added setting to calculate mean or trimmean when extracting signal levels (now set in xr_params)
#Remove hardcoding for xmer boundary. Now set with xm_params.py (xmer_boundary) 
#Changed "heptamer lvl" variable to "xna_region" so it is generalized to any bounding region set by xna_boundary variable
#Removed hard coding of XNA region extraction, XNA signal extraction. Now generalized and set with xm_params 
#Replaced "can_df" with "bed_df" since it is a more descriptive variable name. 
#Fixed hardcoding of XNA relative position, now uses fasta header from Xenomorph. 
#Bed dataframe is no longer used since its only purpose was XNA position, which is now given to us from sequence header. 
#Removed sequence length hard coding by looking up length of the alignment for each read (rather than assuming a fixed length)
#Removed pre-indexing to calculate qscores. Qscores are now calculated by extracted from io_read objects
#Removed pre indexing for map score calculations. Will use match quality scores that are already calculated. 
#Updated to now have correct output headers

"""

#######################################################################



#import everything
import re
import os
import sys
import subprocess

from alive_progress import alive_bar; import time
import pod5
import pysam
import numpy as np
import pandas as pd
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

from remora import io, refine_signal_map, util
from time import sleep
from tqdm import tqdm
from csv import writer
from Bio import SeqIO
from Bio.Seq import Seq
from xm_tools import *
from xr_tools import *
from xm_params import *
from xr_params import *



#Handle input arguments 
##Pod5 file input is handled by coverting raw fast5 to pod5 using pod5 tools
pod5_path = sys.argv[1]+'/'.replace('//','/')

##Bam file is generated from basecalling raw data with a reference alignment and --move_out 
bam_path = sys.argv[2]

##Bed file contains location of XNAs in the raw reference sequence
bed_path = sys.argv[3]

##Path to the reference sequence in fasta format (XNAs converted back to standard DNA for fasta handling) 
fasta_path = sys.argv[4]

##Output folder location 
output_folder = sys.argv[5]

#Import fasta file and bed file 
##Fasta file contains the reference sequences 
fasta_file = pysam.FastaFile(fasta_path)

#Generate output file dataframe
column_names = ["Read_ID","reference_sequence", "read_xna_strand", "read_signal_match_score", "read_q-score","read_reference_locus" ,"read_start", "read_end", "read_xna_position","read_xna", "read_xna_sequence", "read_levels", "read_sd", "read_dwell"]
output_summary = pd.DataFrame(columns = column_names)

#Index fasta reference for future lookup
fasta_ref_dict = SeqIO.index(fasta_path, "fasta")

####POD file operations 
#Import pod5 data into python using pod5.Reader 
pod5_fh = pod5.Reader(pod5_path)


####BAM file operations 
#Index bam file using setting in gen_bai from lib/xr_params.py. Required if input bam file is unindexed. 
if gen_bai==True:
    cmd = 'samtools index ' + bam_path
    os.system(cmd)

#Import bam data using pysam and get header information 
bam_fh = pysam.AlignmentFile(bam_path, 'rb')
headers = bam_fh.header
reference_info = headers['SQ']
contigs = headers.references 
num_reads = pysam.AlignmentFile(bam_path, 'rb').count()



####Sig Map Refiner########
#Set up SigMapRefiner.  
##Loaded 6-mer table with 3 central position. (#ATGC Model Building
##Rough re-scaling will be executed. 
##Signal mapping refinement will be executed using the dwell_penalty refinement method (band half width: 5). 
##Short dwell penalty array set to [8.  4.5 2. ]
sig_map_refiner = refine_signal_map.SigMapRefiner(
    kmer_model_filename=level_table,
    scale_iters=0,
    do_fix_guage=True,
)
############################

enable_rescale = True
if len(sys.argv)==7: 
    print('--------------------------------')
    print('Warning - Rescale parameters is enabled. This estimate is data set specific')



    #Rescale for ATGC
    rescale  = 1.4746563490046836
    reshift = 0.01980030555804524


    print('Rescale = '+str(rescale)+'      reshift = '+str(reshift))
    print('--------------------------------')
else:
    rescale = 1
    reshift = 0 

############################

#Num read overwrite
if max_num_reads >0: 
    num_reads = max_num_reads

count =0 
if 1>0: 

    #Set up progress bar
    with alive_bar(int(num_reads), force_tty=True) as bar: 
        for i in range(0,int(num_reads)):

            #Increment progress bar meter
            bar()

            #Begin read procress
            #if 1>0: 
            try:

                #Iterate through bam read alignments
                bam_read = next(bam_fh)

                #Fetch pod5 data for corresponding bam read 
                pod5_read = next(pod5_fh.reads(selection=[bam_read.query_name]))

                #Load raw data and alignment into io object
                io_read = io.Read.from_pod5_and_alignment(pod5_read, bam_read)

                #Get header of reference sequence that read mapped to
                reference_sequence_header = io_read.ref_reg.ctg

                #Get reference sequence taht read mapped to
                reference_sequence = fasta_ref_dict[reference_sequence_header].seq

                #Get length of reference sequence/alignment
                len_sequence = len(reference_sequence)

                #Set reference region to focus on alignment. Currently hardcoded to be entirety of alignment length 
                ref_reg = io.RefRegion(ctg=io_read.ref_reg.ctg, strand=io_read.ref_reg.strand, start=0, end=len_sequence)

                #Refine raw signal scaling around the focus region using sig_map_aligner 
                io_read.set_refine_signal_mapping(sig_map_refiner, ref_mapping=True)

                #Get name of raw read 
                read_ID = io_read.read_id


                #Get identity of strand: (+) or (-) 
                strand = io_read.ref_reg.strand

                #Get read map start position
                read_map_start = int(io_read.ref_reg.start)

                #Get read map end position
                read_map_end = int(io_read.ref_reg.end)

                #Get match quality of the full read alignment 
                match_quality=io_read.full_align["map_quality"]

                #Get qscore from read phred score
                read_qs = phred_to_qscore(io_read.full_align["qual"])

                #Get XNA positions from reference read header
                xna_base_pos = fetch_xna_pos(reference_sequence_header)

                #Get identity of XNA base
                xna_base=reference_sequence[extract_pos]

                #Get position of XNA 
                xna_pos=extract_pos

                #Extract read levels if read passes quality filters set in xr_params.py 
                if int(match_quality) >= min_match_score:
                    if strand == '+':

                        #Extract signal level and perform arithmetic computation (mean or trimmean, set in xr_params)
                        read_dms = io_read.compute_per_base_metric("dwell_trimmean_trimsd", ref_anchored=True, signal_type = "norm",region=ref_reg, sig_map_refiner = sig_map_refiner)

                        #Get signal means for whole read
                        read_means = read_dms["trimmean"]*rescale+reshift

                        #Extract sub-region defined by xmer_padding for filtering (set in xm_params ; xmer_padding). This is used for segmentation check.
                        read_region = read_dms["trimmean"][xna_pos-xmer_padding:xna_pos+1+xmer_padding]

                        #Extract signal bounding region of interest for level output file (defualt = +/- 3. set in xm_params; xmer_boundary)
                        xna_region_mean = read_means[xna_pos-xmer_boundary:xna_pos+1+xmer_boundary]

                        #Get signal means for whole read
                        xna_region_sd = read_dms["trimsd"][xna_pos-xmer_boundary:xna_pos+1+xmer_boundary]

                        #Get signal means for whole read
                        xna_region_dwell = read_dms["dwell"][xna_pos-xmer_boundary:xna_pos+1+xmer_boundary]

                        #Extract sequence bounding region of interest for level output file (defualt = +/- 3. set in xm_params; xmer_boundary)
                        read_xna_seq = reference_sequence[xna_pos-xmer_boundary:xna_pos] + xna_base + reference_sequence[xna_pos+1:xna_pos+xmer_boundary+1]

                        #Reference locus 
                        read_ref_locus = read_xna_seq

                    elif strand == '-':

                        #Reference locus (before taking reverse complement, and before taking rc of the xna base)
                        read_ref_locus = reference_sequence[xna_pos-xmer_boundary:xna_pos] + xna_base + reference_sequence[xna_pos+1:xna_pos+xmer_boundary+1]

                        #If alignment matched (-) strand, use r.c. of XNA base 
                        try: 
                            xna_base = xna_base_rc(xna_base, xna_base_pairs) 
                        except: 
                            xna_base = xna_base_rc(xna_base, standard_base_pairs)

                        #If alignment matched (-) strand, use r.c. of reference sequence
                        reference_sequence_rc = reference_sequence.reverse_complement()

                        #Extract signal level and perform arithmetic computation (mean or trimmean, set in xr_params)
                        read_dms = io_read.compute_per_base_metric("dwell_trimmean_trimsd", ref_anchored=True, signal_type = "norm", region=ref_reg, sig_map_refiner = sig_map_refiner)

                        #Get signal means for whole read
                        read_means = read_dms["trimmean"]*rescale+reshift

                        #Extract sub-region defined by xmer_padding for filtering (set in xm_params ; xmer_padding). This is used for segmentation check.
                        read_region = read_means[len_sequence-xna_pos-xmer_padding-1:len_sequence-xna_pos+xmer_padding]

                        #Extract signal bounding region of interest for level output file (defualt = +/- 3. set in xm_params; xmer_boundary)
                        xna_region_mean = read_means[len_sequence-xna_pos-xmer_boundary-1:len_sequence-xna_pos+xmer_boundary]

                        #Get signal means for whole read
                        xna_region_sd = read_dms["trimsd"][len_sequence-xna_pos-xmer_boundary-1:len_sequence-xna_pos+xmer_boundary]

                        #Get signal means for whole read
                        xna_region_dwell = read_dms["dwell"][len_sequence-xna_pos-xmer_boundary-1:len_sequence-xna_pos+xmer_boundary]

                        #Extract sequence bounding region of interest for level output file (defualt = +/- 3. set in xm_params; xmer_boundary)
                        read_xna_seq = reference_sequence_rc[len_sequence-xna_pos-xmer_boundary-1:len_sequence-xna_pos-1]+xna_base+reference_sequence_rc[len_sequence-xna_pos:len_sequence-xna_pos+xmer_boundary]



                    #Use signals extracted from read region to check if segmentation passed around region of interest
                    if np.isnan(sum(read_region))==False: 

                        #Store output dataframe for level file generation
                        read_out = {'Read_ID':read_ID, "reference_sequence": reference_sequence_header, "read_xna_strand": strand, "read_reference_locus": read_ref_locus, "read_start": read_map_start, "read_end":read_map_end, "read_levels": xna_region_mean, "read_sd": xna_region_sd, "read_dwell": xna_region_dwell, "read_q-score": read_qs, "read_signal_match_score": match_quality,"read_xna_position": xna_pos, "read_xna": xna_base, "read_xna_sequence": read_xna_seq}


                        output_summary.loc[len(output_summary)] = read_out

                    #If there are missing values around segmentation region, skip processing read
                    else: 
                        segmentation_failed+=1 


            except:
                pass


    #Generate reporting summary 
    passed_reads = len(output_summary)
    failed_reads = int(num_reads)-passed_reads
    print("[Xemora Status] Analyzed " + str(num_reads) + " reads.\n"+ str(passed_reads) + " passed alignment and segmentation.\n" + str(failed_reads) + " reads did not align. \n" +str(segmentation_failed)+" reads failed to segment properly")


    #Save output to csv file
   # output_summary.to_csv(output_folder, sep='\t', encoding='utf-8',index=False)
    output_summary.to_csv(output_folder, encoding='utf-8',index=False)


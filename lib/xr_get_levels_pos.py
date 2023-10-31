########################################################################
########################################################################
"""
xr_get_levels_pos.py 
Remora based segmentation signal extraction for Xenomorph. 
Uses a fixed position within a read (specified in xr_params.py) for extracting
XNA signals. 

Cite us: 

H. Kawabe, C. Thomas, S. Hoshika, Myong-Jung Kim, Myong-Sang Kim, L. Miessner, N. Kaplan, J. M. Craig, 
J. Gundlach, A. Laszlo,  S. A. Benner, J. A. Marchand. "Enzymatic Synthesis and Nanopore Sequencing of 12-Letter Supernumerary DNA" 
Nature Communications. 14. (2023). DOI: 10.1038/s41467-023-42406-z 

Updated: 8/19/23 
"""

#######################################################################

#Import packages and parameter files
import re
import os
import sys
import subprocess
import pod5
import pysam
import numpy as np
import pandas as pd
from alive_progress import alive_bar; import time
from remora import io, refine_signal_map, util
from Bio import SeqIO
from Bio.Seq import Seq
from xm_tools import *
from xr_tools import *
from xm_params import *
from xr_params import *

if __name__ == '__main__': 
	#Print status
	print('Xenomorph Status - [Preprocess] Initializing level extraction for fixed position (not typical)')

	#Handle input arguments 
	##Pod5 file input is handled by coverting raw fast5 to pod5 using pod5 tools
	pod5_path = os.path.normpath(sys.argv[1])+'/'

	##Bam file is generated from basecalling raw data with a reference alignment and --move_out 
	bam_path = sys.argv[2]

	##Path to the reference sequence in fasta format (XNAs converted back to standard DNA for fasta handling) 
	fasta_path = sys.argv[3]

	##Output folder location 
	output_folder = sys.argv[4]

	#Check if reverse complement flag is present, if so use rc 
	if len(sys.argv)==6: 
		if sys.argv[5] == 'rc':
		    output_summary = pd.read_csv(output_folder)

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

	#Segmentation failed count
	segmentation_failed = 0

	#Index bam file using setting in gen_bai from lib/xr_params.py. Required if input bam file is unindexed. 
	if gen_bai==True:
		cmd = 'samtools index ' + bam_path
		os.system(cmd)

	#Generate file path for a primary alignment bam output file 
	primary_bam_path = bam_path.replace('.bam','_primary.bam')

	#Extract only reads with primary alignments and put them in a new bam file 
	bam_fh = filter_primary_alignments_file(bam_path, primary_bam_path)

	#Get headers of bam file 
	headers = bam_fh.header

	#Get info on reference sequences
	reference_info = headers['SQ']

	#Get contigs 
	contigs = headers.references 

	#Get total number of reads in primary only bam output 
	num_reads = pysam.AlignmentFile(primary_bam_path, 'rb').count()


	#Num read overwrite
	if max_num_reads >0: 
		num_reads = max_num_reads

	##Set global scaling shift and scale
	if manual_rescale_override == True:
		rescale = manual_rescale
		reshift = manual_reshift
	else: 
		rescale = global_rescale
		reshift = global_reshift

	#Perform global rescaling estimate on ATGC portions of read
	if len(sys.argv)==6: 
		if sys.argv[5]=='rescale':
		    print('Xenomorph Status - [Preprocess] Extracting kmers for calculating global scaling paramters')
		    #Reinitialize slope at 1 (no scaling) for rescale calculation
		    rescale = 1

		    #Reinitialize shift at 0 (no scaling) for rescale calculation
		    reshift = 0 
		    
		    #Number of levels before and after to extract surrounding an XNA (default = 3) 
		    xmer_boundary = rescale_xmer_boundary

		    #Number of bases before and after XNA that are required in a matching read (default = 30 alt) 
		    xmer_padding = rescale_xmer_padding
		    
		    #Maximum number of reads to process
		    num_reads = rescale_max_num_reads

	print('Xenomorph Status - [Preprocess] Rescaling using the following parameters: m = '+str(rescale)+'   b = '+str(reshift))



	#Set up progress bar
	with alive_bar(int(num_reads), force_tty=True) as bar: 
		for i in range(0,int(num_reads)):

		    #Increment progress bar meter
		    bar()

		    #Begin read procress
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

		        #Check if RC is allowed for XNA 
		        rc_allowed = xna_base_rc(xna_base, xna_segmentation_model_sets) 

            
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

		            elif strand == '-' and rc_allowed != False:
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

		                #Store to output summary dataframe
		                output_summary.loc[len(output_summary)] = read_out

		            #If there are missing values around segmentation region, skip processing read
		            else: 
		                segmentation_failed+=1 
		    except:
		        pass


	#Generate reporting summary 
	passed_reads = len(output_summary)
	failed_reads = int(num_reads)-passed_reads
	print("Xenomorph Status - [Preprocess] Analyzed " + str(num_reads) + " reads")
	print("Xenomorph Status - [Preprocess] Number of reads passed alignmnet: "+str(passed_reads))
	print("Xenomorph Status - [Preprocess] Number of reads failed to segment: "+str(segmentation_failed))
	print("Xenomorph Status - [Preprocess] Number of reads failed for unknown reason: "+str(failed_reads-segmentation_failed))

	#Save output to csv file
	output_summary.to_csv(output_folder, encoding='utf-8',index=False)


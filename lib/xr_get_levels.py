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
from xm_tools import *
from xr_tools import *
from xm_params import *
from xr_params import *

####Changes 230630
#Removed redundant BAM processing code 
#Fixed bug in read count assignment and changed var to num_reads instead of num_file 
#Added optional qscore extraction to properly fall back if needed
#Annotated sig map refiner 
#Added code comments




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

##Convert bed file to a dataframe for easier handling 
can_df = pd.read_csv(bed_path, sep='  ', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'], engine = 'python')
can_df['chrom'] = can_df['chrom'].str.replace('\\n', '')


#Generate output file dataframe
column_names = ["Read ID","Reference_seq", "strand", "read_signal_match_score", "read_q-score", "read_map_score", "read_xna", "read_xna_sequence", "read_levels"]
output_summary = pd.DataFrame(columns = column_names)

##XXXXXX## Hard coded 
##From bed file, read positions of XNAs within each reference. Currently this is hard coded since star/end can vary from read to read. 
start_pos = can_df['start'][0]
end_pos = can_df['end'][0]

##XXXXXX## Hard coded
#Get length of each sequence in the reference file
fasta_sequences = SeqIO.parse(open(fasta_path),'fasta')
for fasta in fasta_sequences:
    len_sequence =len(fasta.seq)


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

#Get quality scores and mapping scores. Note: This can be slow since it requires iterate. 
if skip_qscore_extract == False: 
    print('Xenomorph Status - [Preprocess] Re-indexing bam file for extracting q-score and map score.')
    column_names = ["Read ID","match_score", "q-score"]
    score_df = pd.DataFrame(columns = column_names)
    score_df = index_read_scores(pysam.AlignmentFile(bam_path, 'rb'), score_df)
    print('Xenomorph Status - [Preprocess] Re-indexing complete.')
else: 
    print('Xenomorph Status - [Preprocess] Skipping extraction of qscore and mapping score from bam file.')
    print('Xenomorph Status - [Preprocess] Defaulting to (q-score = 15, map_score = 10)')
    qscore = 15
    map_score = 10




####Sig Map Refiner
#Set up SigMapRefiner.  
##Loaded 6-mer table with 3 central position. (NNXNNN)
##Rough re-scaling will be executed. 
##Signal mapping refinement will be executed using the dwell_penalty refinement method (band half width: 5). 
##Short dwell penalty array set to [8.  4.5 2. ]
sig_map_refiner = refine_signal_map.SigMapRefiner(
    kmer_model_filename=level_table,
    scale_iters=0,
    do_fix_guage=True,
)



#Get number of total reads
#print(len(score_df))
#print(dir(bam_fh))
#print(num_reads)




if 1>0: 
    #basecall anchored metrics
    #bam_read = next(bam_fh)
    #pod5_read = next(pod5_fh.reads(selection=[bam_read.query_name]))





    j=0
    with alive_bar(int(num_reads), force_tty=True) as bar: #sets bar visually for progress
        for i in range(0,int(num_reads)):
        #for i in range(0,1000):
            bar()


            #if 1>0: 
            try:
                #go to the next bam and can read
                bam_read = next(bam_fh)
                pod5_read = next(pod5_fh.reads(selection=[bam_read.query_name]))



                #load read, reference region, reference sequence
                io_read = io.Read.from_pod5_and_alignment(pod5_read, bam_read)
                ref_reg = io.RefRegion(ctg=io_read.ref_reg.ctg, strand=io_read.ref_reg.strand, start=0, end=len_sequence)
                reference_sequence = io_read.ref_reg.ctg

                io_read.set_refine_signal_mapping(sig_map_refiner, ref_mapping=True)
                read_ID = io_read.read_id
                strand = io_read.ref_reg.strand
                match_quality=io_read.full_align["map_quality"]
                #set qscore and match_score
          



                if skip_qscore_extract == False: 
                    qscore = score_df.loc[score_df['Read ID'] == read_ID, 'q-score'].values[0]
                    map_score = score_df.loc[score_df['Read ID'] == read_ID, 'match_score'].values[0]

                if int(match_quality) >=match_score:
                    if strand == '+':
                        XNA = mod_base
                        read_means = io_read.compute_per_base_metric("dwell_trimmean", ref_anchored=True, region=ref_reg, sig_map_refiner = sig_map_refiner)["trimmean"]
                        heptamer_lvl = read_means[start_pos-3:end_pos+3]
                        read_xna_seq = "AGTPCAG"

                    else:
                        XNA = mod_rev_base
                        read_means = io_read.compute_per_base_metric("dwell_trimmean", ref_anchored=True, region=ref_reg, sig_map_refiner = sig_map_refiner)["trimmean"][::-1]
                        read_hep = io_read.compute_per_base_metric("dwell_trimmean", ref_anchored=True, region=ref_reg, sig_map_refiner = sig_map_refiner)["trimmean"]
                        heptamer_lvl = read_hep[len_sequence-start_pos-4:len_sequence-start_pos+3]
                        read_xna_seq = "CTGZACT"


                    read_out = {'Read ID':read_ID, "Reference_seq": reference_sequence, "strand": strand, "read_levels": heptamer_lvl, "read_q-score": qscore, "read_map_score": map_score, "read_signal_match_score": match_quality, "read_xna": XNA, "read_xna_sequence": read_xna_seq}
                    output_summary.loc[len(output_summary)] = read_out

                    read_ref_reg = io_read.extract_ref_reg(ref_reg)

            except:
                j+=1
                pass

    passed_reads = int(num_reads)-j
    failed_reads = j

    output_summary.to_csv(output_folder, sep='\t', encoding='utf-8',index=False)
    print("[Xemora Status] Analyzed " + str(num_reads) + " reads.\n"+ str(passed_reads) + " passed.\n" + str(failed_reads) + " reads did not align.")


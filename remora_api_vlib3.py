#this is a pipeline that uses Remora to extract levels for the heptamers.
#This code starts with the pod5 and bam files already generated as well as the xfasta and is meant for troubleshooting
#If you want to run the version that is integrated into Xenomorph (start with fast5 files) please use the xx-test_pipe.py file

#--------------------------------------------------------------------------------

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

#--------------------------------------------------------------------------------

# data locations
'''
mod_pod5_path = '/home/marchandlab/xenomorph/xemora-beta/xemora-test/230523_signal_extract/modified/pod5/230304_PZ_libv4_1k.pod5'
mod_bam_path = '/home/marchandlab/xenomorph/xemora-beta/xemora-test/230523_signal_extract/modified/bam/bam.bam'
mod_bed_path = '/home/marchandlab/xenomorph/xemora-beta/xemora-test/230523_signal_extract/references/P.bed'

fasta_path = '/home/marchandlab/xenomorph/xemora-beta/xemora-test/230523_signal_extract/references/xref_libv2_PZ_CxDx.fa'
output_folder= '/home/marchandlab/xenomorph/xemora-beta/xemora-test/230523_signal_extract/levels.csv'
'''

can_pod5_path = '/home/marchandlab/xenomorph/xenomorph-xemora/xx-test/230625_PZ_libv2_remoraAPItesting/modified/pod5/fast5.pod5'
can_bam_path = '/home/marchandlab/xenomorph/xenomorph-xemora/xx-test/230625_PZ_libv2_remoraAPItesting/canonical/bam/bam.bam'
bed_path = '/home/marchandlab/xenomorph/xemora-beta/xemora-test/230309_Xemora_PZ_GC_train/references/P.bed'

fasta_path = '/home/marchandlab/xenomorph/xenomorph-xemora/xx-test/230625_PZ_libv2_remoraAPItesting/references/xref_libv2_PZ_CxDx.fa'
output_folder= '/home/marchandlab/xenomorph/xenomorph-xemora/xx-test/230625_PZ_libv2_remoraAPItesting/230625PZlib_levels.csv'

gen_bai = True

#find length of sequences to input for reference region
#this could be coded better so it only takes the first sequence 
fasta_sequences = SeqIO.parse(open(fasta_path),'fasta')
for fasta in fasta_sequences:
    len_sequence =len(fasta.seq)

#if you don't have your bam file indexed please run the following code to index (remora will give an error saying you're missing index)
#this is something we're going to have to integrate in the overall pipeline as well
if gen_bai==True:
    cmd = 'samtools index ' + can_bam_path
    os.system(cmd)

#read in your pod5 and bam files
can_pod5_fh = pod5.Reader(can_pod5_path)
can_bam_fh = pysam.AlignmentFile(can_bam_path)

#import in kmer table from remora, I've left both 10.4.1 and 9.4.1 kmers below
#level_table = '/home/marchandlab/xenomorph/xemora-beta/xemora-test/230523_signal_extract/references/levels.txt'
level_table = '/home/marchandlab/xenomorph/xemora-beta/xemora-test/kmers/9.4_6mers_450bps_trimmed.txt'

#bam file and contig definition, unclear as to what exactly it does
bam_file = pysam.AlignmentFile(can_bam_path, 'rb')
contigs = bam_file.header.references


#signal map refiner, this is what helps us extract the levels
sig_map_refiner = refine_signal_map.SigMapRefiner(
    kmer_model_filename=level_table,
    scale_iters=0,
    do_fix_guage=True,
)


#find number of files
cmd = 'samtools view -c ' + can_bam_path
os.system(cmd)
file_num = subprocess.check_output(cmd, shell=True, text=True)


#basecall anchored metrics
bam_read = next(can_bam_fh)
pod5_read = next(can_pod5_fh.reads(selection=[bam_read.query_name]))


#fnd position of XNA
fasta_file = pysam.FastaFile(fasta_path)
can_df = pd.read_csv(bed_path, sep='  ', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'], engine = 'python')
can_df['chrom'] = can_df['chrom'].str.replace('\\n', '')
start_pos = can_df['start'][0]
end_pos = can_df['end'][0]
#print(start_pos)

fasta_sequences = {}
fasta_reverse_complements = {}

#index all fasta sequences so that they can be referred back when extracting the heptamer region
with open(fasta_path, "r") as file:
    sequence = ""
    sequence_name = ""

    for line in file:
        if line.startswith(">"):
            if sequence_name != "":
                fasta_sequences[sequence_name] = sequence
                fasta_reverse_complements[sequence_name] = sequence[::-1].translate(str.maketrans("ACGT", "TGCA"))

            sequence = ""
            sequence_name = line.strip()[1:]
        else:
            sequence += line.strip()

    # Store the last sequence
    fasta_sequences[sequence_name] = sequence
    fasta_reverse_complements[sequence_name] = sequence[::-1].translate(str.maketrans("ACGT", "TGCA"))

#create output file
column_names = ["Read ID","Reference_seq", "strand","read_levels", "read_signal_match_score", "read_q-score", "read_map_score", "read_xna", "read_xna_sequence"]
output_summary = pd.DataFrame(columns = column_names)

#number of reads that didn't align and also a match score used to filter reads (currently set to 0)
j=0
match_score = 0

#if you only want to process a subset of data
use_subset = True
subset = 1000

with alive_bar(int(file_num), force_tty=True) as bar: #sets bar visually for progress
    #for i in range(0,int(file_num)):
    for i in range(0,subset):
        bar()

        try:
            #go to the next bam and can read
            bam_read = next(can_bam_fh)
            pod5_read = next(can_pod5_fh.reads(selection=[bam_read.query_name]))

            #load read, reference region, reference sequence, and other parameter
            io_read = io.Read.from_pod5_and_alignment(pod5_read, bam_read)
            ref_reg = io.RefRegion(ctg=io_read.ref_reg.ctg, strand=io_read.ref_reg.strand, start=0, end=len_sequence)
            reference_sequence = io_read.ref_reg.ctg
            io_read.set_refine_signal_mapping(sig_map_refiner, ref_mapping=True)
            read_ID = io_read.read_id
            strand = io_read.ref_reg.strand
            match_quality=io_read.full_align["map_quality"]
            read_name = io_read.ref_reg.ctg
            fasta_sequence = fasta_sequences.get(read_name)
            reverse_complement = fasta_reverse_complements.get(read_name)

            #this is a potential filter by match score but it's currently set to 0
            if int(match_quality) >=match_score:

                #extract the levels of the heptamer region and the corresponding heptamer sequence for + and - strand
                if strand == '+':
                    XNA = "P"
                    read_means = io_read.compute_per_base_metric("dwell_trimmean", ref_anchored=True, region=ref_reg, sig_map_refiner = sig_map_refiner)["trimmean"]
                    heptamer_lvl = read_means[start_pos-3:end_pos+3]
                    read_xna_seq = fasta_sequence[start_pos-3:start_pos] + XNA + fasta_sequence[start_pos+1:start_pos+4]
                else:
                    XNA = "Z"
                    read_means = io_read.compute_per_base_metric("dwell_trimmean", ref_anchored=True, region=ref_reg, sig_map_refiner = sig_map_refiner)["trimmean"][::-1]
                    read_hep = io_read.compute_per_base_metric("dwell_trimmean", ref_anchored=True, region=ref_reg, sig_map_refiner = sig_map_refiner)["trimmean"]
                    heptamer_lvl = read_hep[len_sequence-start_pos-4:len_sequence-start_pos+3]
                    read_xna_seq = reverse_complement[len_sequence-start_pos-4:len_sequence-start_pos-1] + XNA + reverse_complement[len_sequence-start_pos:len_sequence-start_pos+3]

                #store all parameters 
                read_out = {'Read ID':read_ID, "Reference_seq": reference_sequence, "strand": strand, "read_levels": heptamer_lvl, "read_q-score": 10, "read_map_score": match_quality, "read_signal_match_score": 1, "read_xna": XNA, "read_xna_sequence": read_xna_seq}

                output_summary.loc[len(output_summary)] = read_out
                read_ref_reg = io_read.extract_ref_reg(ref_reg)
        except:
            j+=1
            pass

if use_subset == True:
    passed_reads = subset-j
else:
    passed_reads = int(file_num)-j
failed_reads = j

#save output and print status message
output_summary.to_csv(output_folder, index=False)
if use_subset == True:
    print("[Xemora Status] Analyzed " + str(subset) + " reads.\n"+ str(passed_reads) + " passed.\n" + str(failed_reads) + " reads did not align.")
else:
    print("[Xemora Status] Analyzed " + str(file_num) + " reads.\n"+ str(passed_reads) + " passed.\n" + str(failed_reads) + " reads did not align.")

########################################################################
########################################################################
"""
xm_get_levels.py

Description: Extracts levels from kmers of a given read, and performs
 XNA replacement based on standard library. This is part of the Xenomorph
preprocessing pipeline and is automatically called by xenomorph.py preprocess. 
Requires resquiggled dataset and reference file as inputs. 


Title: Synthesis and Sequencing of 12-Letter Supernumerary DNA

By: H. Kawabe, C. Thomas, S. Hoshika, Myong-Jung Kim, Myong-Sang Kim, L. Miessner, J. M. Craig, 
J. Gundlach, A. Laszlo,  S. A. Benner, J. A. Marchand

Updated: 2/14/23
"""
########################################################################
########################################################################




from Bio.Seq import Seq
import sys
import numpy as np
import pandas as pd
import itertools
import os
from tombo import tombo_helper, tombo_stats, resquiggle
from string import ascii_lowercase
from xm_params import *
from xm_tools import *
from alive_progress import alive_bar




#User inputs 
fast5_basedir = [os.path.normpath(sys.argv[1])]
reference_fasta = os.path.normpath(sys.argv[2])
output_folder = os.path.normpath(sys.argv[3])


#Create new output file 
output_column_names=['read_ID','reference_sequence','read_q-score','read_signal_match_score','read_signal_match_start','read_signal_match_end','read_pos_relative_to_raw','read_reference_locus','read_substituted_base','read_xna_strand','read_xna_sequence','read_xna','read_xna_position','read_levels']
output_summary = pd.DataFrame(columns = output_column_names)

#Set up corrected group read if required 
if len(sys.argv)==5: 
    corr_group= 'RawGenomeCorrected_000'+sys.argv[4]
    if sys.argv[4][len(sys.argv[4])-1]=='r':
        #Merge with previous output file  Load previous. Always true if rc taken as input. 
        output_summary = pd.read_csv(output_folder)
else: 
    corr_group= 'RawGenomeCorrected_000' 


#Parse Tombo index from previously re-squiggled set of reads
reads_index = tombo_helper.TomboReads(fast5_basedir, corrected_group=corr_group)


#Print status
print('Xenomorph Status - [Preprocess] Performing level extraction on '+corr_group+' read set.') 

#Get list of all chromosomes
clist = reads_index.get_all_cs()

#Sort list for improved reproduceabiity in troubleshooting
clist.sort()

#Loop through all alignments
with alive_bar(len(clist), force_tty=True) as bar: 
    for i in range(0,len(clist)):
        bar()

        #Get read ID
        echrm = clist[i][0]

        #Get strand
        s = clist[i][1]

        #Enforce strand selection 
        if '_Gap' not in echrm:
            try:
            #if 1>0:

                xna_base_pos = fetch_xna_pos(echrm)
                #For every xna position, extract levels


                for x in xna_base_pos:


                    #Set xna base and xna position in read alignment
                    xna_base=x[0]
                    xna_pos=int(x[1])

                    #Check if RC is allowed for XNA 
                    rc_allowed = xna_base_rc(xna_base, xna_segmentation_model_sets) 



                    if rc_allowed != False: 
                        check_strand=['+','-']
                    else: 
                        check_strand=['+']

                    #Only certain strands might be allowed depending on model segmentation rules
                    if s in check_strand: 
                        #Pick data from strand, and specify motif between 30 and 65 bp of reference
                        reg_data = tombo_helper.intervalData(chrm=echrm, start=xna_pos-xmer_padding, end=xna_pos+xmer_padding+1, strand=s)
                        

                        #Extract sequence/fa information for the mapped read
                        reg_bases = reg_data.add_seq(tombo_helper.Fasta(reference_fasta))


                        #Extract levels and force full span (set false to not force full span)
                        reg_data.add_reads(reads_index,True)


                        if len(reg_data.reads)>0:

                            reg_base_levels = reg_data.get_base_levels()              

                            reg_base_levels = reg_base_levels[xmer_padding-xmer_boundary: xmer_padding+xmer_boundary+1]
                                #Substitue XNA
                            ref_seq_full = reg_bases.seq[xmer_padding-xmer_boundary:xmer_padding]+xna_base+reg_bases.seq[xmer_padding+1:xmer_padding+1+xmer_boundary]



                            if s=='+':
                                #Set XNA
                                X = xna_base

                                #Substitue XNA
                                seq_full = reg_bases.seq[xmer_padding-xmer_boundary:xmer_padding]+X+reg_bases.seq[xmer_padding+1:xmer_padding+1+xmer_boundary]

                                #Base that was previously substituted - Store this as a column
                                sub_base = reg_bases.seq[xmer_padding]


                            if s=='-':
                                #Set XNA
                                X = xna_base_rc(xna_base, xna_base_pairs) ##This needs to be set somehow

                                #Base that was previously substituted - Store this as a column
                                sub_base = xna_base_rc(reg_bases.seq[xmer_padding], standard_base_pairs)

                                #Reverse base levels
                                reg_base_levels = reg_base_levels[::-1]

                                #Take complement of sequence
                                reg_bases_c = str(Seq(reg_bases.seq).complement())

                                #Perform substitution
                                seq_full_c = reg_bases_c[xmer_padding-xmer_boundary:xmer_padding]+X+reg_bases_c[xmer_padding+1:xmer_padding+1+xmer_boundary]

                                #Next reverse sequence with XNA substitution
                                seq_full = str(seq_full_c)[::-1]



                            #For each read, extract levels and sequence
                            output_summary_temp = pd.DataFrame(columns = output_column_names)
                            for rid in range(0,len(reg_data.reads)):
                                #Read metadata
                                read_id = reg_data.reads[rid].read_id
                                read_fn = reg_data.reads[rid].fn
                                read_qscore = reg_data.reads[rid].mean_q_score
                                read_sscore = reg_data.reads[rid].sig_match_score
                                read_start = reg_data.reads[rid].start
                                read_end = reg_data.reads[rid].end
                                read_rel = reg_data.reads[rid].read_start_rel_to_raw

                                #Sequence - level - XNA
                                read_ref_locus = ref_seq_full
                                read_sub_base = sub_base
                                read_xna_strand = s 
                                read_seq = seq_full
                                read_xna = X
                                read_xna_pos = xna_pos
                                read_level = np.transpose(reg_base_levels)[rid]
                                read_out=[read_id, echrm, read_qscore, read_sscore, read_start, read_end, read_rel, read_ref_locus,read_sub_base, read_xna_strand, read_seq, read_xna, read_xna_pos, read_level]
                                #print(read_out)
                                output_summary_temp.loc[len(output_summary_temp)]=read_out
                            output_summary  = output_summary.append(output_summary_temp, ignore_index = True)
     




            except:
                print('Error processing read. Skipping '+echrm+'.')


print('Xenomorph Status - Extracted '+str(len(output_summary))+' sequence levels with XNAs')
output_summary.to_csv(output_folder, index = False)





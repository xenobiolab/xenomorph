########################################################################
########################################################################
"""
xm_get_raw.py

Description: This is a sample script used to generate output files that 
contain unsegmented raw and non-normalized signals, unsegmented and normalized
signals, segmented and non-normalized signals. Uses tombo pipeline for 
extraction and alignment of signal events. Requires specifying region to 
extract data from. Requires user inputs. 

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
import h5py
import matplotlib.pyplot as plt
from alive_progress import alive_bar
import mappy



##############Handle user inputs ##############
fast5_basedir = [sys.argv[1].replace('//','/')] #
reference_fasta = sys.argv[2].replace('//','/') #
level_input_file = sys.argv[3].replace('//','/') #
output_folder = sys.argv[4].replace('//','/') #

##############XNA specific paramters##############
#Sequence bases upstream of heptamer
heptamer_us = 'ACTGCTGA' 

#Sequence bases downstream of heptamer 
heptamer_ds = 'TCAGCAGT' 

#Number of additional bases before and after heptamer to extract
heptamer_padding = 0 

#Tombo index corrected group (default)
corr_group= 'RawGenomeCorrected_000' 

#XNA basepairs to extract signal from 
base_pairs = ['P','Z']

#Strands to check for XNA basepairs 
strands=['+','-']
#################################################





#Print status
print('Xenomorph Status - [Raw] Performing raw level extraction on '+corr_group+' read set.') 


output_column_names=['read_ID','sequence','raw_signal', 'norm_raw_signal']
output_summary = pd.DataFrame(columns = output_column_names)


##############functino  ##############

def get_mean_signal(raw_signal, segs): 
    mean_sig = [] 
    for ii in range(0,len(segs)-1): 
        mean_sig.append(np.mean(raw_signal[segs[ii]:segs[ii+1]]))
    return(mean_sig)


##############Choose a source of fast5 files to analyze##############
#Option 1 - Get list of fast5 single files (from a folder)
fast5_fn_list = os.listdir(fast5_basedir[0])

#Option 2 - Get list of fast5 single files from a level output file (from xenomorph.py preprocess) 
level_file = pd.read_csv(level_input_file)
level_file_list = level_file['read_ID'].to_list()
fast5_fn_list = level_file_list



with alive_bar(len(fast5_fn_list), force_tty=True) as bar: 
    for i in range(0,len(fast5_fn_list)): 
        bar()
        try: 
            #Get fast5 file 
            fast5_fn = os.path.normpath(fast5_basedir[0])+'/'+fast5_fn_list[i]+'.fast5'

            #Set up fast5 reference fasta
            reference_fn = reference_fasta

            #Get fast5 data
            fast5_data = h5py.File(fast5_fn, 'r')

            #Get sample type, should be DNA for params
            seq_samp_type = tombo_helper.get_seq_sample_type(fast5_data)

            ####### Aligner 
            # prep aligner, signal model and parameters
            aligner = mappy.Aligner(reference_fn, preset=str('map-ont'), best_n=1)

            #Set up rerference for aligner 
            std_ref = tombo_stats.TomboModel(seq_samp_type=seq_samp_type)

            #Set up resquiggle params based on DNA default
            rsqgl_params = tombo_stats.load_resquiggle_parameters(seq_samp_type)


            #######extract data from FAST5
            map_results = resquiggle.map_read(fast5_data, aligner, std_ref)

            all_raw_signal = tombo_helper.get_raw_read_slot(fast5_data)['Signal'][:]
            if seq_samp_type.rev_sig:
                all_raw_signal = all_raw_signal[::-1]
            map_results = map_results._replace(raw_signal=all_raw_signal)



            #Run resquiggle (optional) 
            rsqgl_results = resquiggle.resquiggle_read(map_results, std_ref, rsqgl_params, all_raw_signal=all_raw_signal)



            ####### Individual steps
            num_events = tombo_stats.compute_num_events(
                all_raw_signal.shape[0], len(map_results.genome_seq),
                rsqgl_params.mean_obs_per_event)
            valid_cpts, norm_signal, scale_values = resquiggle.segment_signal(
                map_results, num_events, rsqgl_params)
            event_means = tombo_stats.compute_base_means(norm_signal, valid_cpts)
            dp_results = resquiggle.find_adaptive_base_assignment(
                valid_cpts, event_means, rsqgl_params, std_ref, map_results.genome_seq)
            norm_signal = norm_signal[
                dp_results.read_start_rel_to_raw:
                dp_results.read_start_rel_to_raw + dp_results.segs[-1]]


            #Raw signal 
            raw_sig = all_raw_signal[dp_results.read_start_rel_to_raw:dp_results.read_start_rel_to_raw + dp_results.segs[-1]]

            #Normalized signal
            norm_sig = norm_signal 

            #Segments for each base
            segs = resquiggle.resolve_skipped_bases_with_raw(dp_results, norm_signal, rsqgl_params)

            #sequence of each base 
            mapped_sequence = dp_results.genome_seq


            #Sets maps 1:1 with mapped_sequence. Search position in dp_results.genome_seq
            st = mapped_sequence.find(heptamer_us)+len(heptamer_us)-heptamer_padding
            en = mapped_sequence.find(heptamer_ds)+heptamer_padding


            #########Outputs#########
            #Get heptamer sequence
            #Sequence of range
            ex_seq = mapped_sequence[st:en]


            #Segmentation of range (subtract minimum to re-center) 
            ex_seg = segs[st:en+1]-segs[st]

            #Get strand it mapped to 
            strand = map_results.genome_loc.Strand

            #Get corresponding XNA base of the strand
            xna_base = base_pairs[strands.index(strand)]

            if len(ex_seq)==7+(2*heptamer_padding):
                #Raw signal for given range
                ex_raw_sig = raw_sig[segs[st]:segs[en+1]]

                #Normalized raw signal for given range (???)
                ex_nor_raw_sig = rsqgl_results.raw_signal[segs[st]:segs[en+1]]

                #Normalized signal for given range (Segmented and noramlized)
                ex_nor_sig = norm_sig[segs[st]:segs[en+1]]

                #Get position of XNA (always centered in this analysis) 
                xna_pos = int((len(ex_seq)-1)/2)

                #Generate sequence with XNA
                xna_ex_seq = ex_seq[0:xna_pos]+xna_base+ex_seq[xna_pos+1:]

                output_summary.loc[len(output_summary)] = [fast5_fn.split('/')[2], xna_ex_seq, ex_raw_sig,ex_nor_raw_sig]




        except: 
            n = 0


output_summary.to_csv(output_folder)


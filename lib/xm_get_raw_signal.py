########################################################################
########################################################################
"""
xm_get_raw_signal.py

Description: This is a sample script used to generate output files that 
contain unsegmented raw and non-normalized signals, unsegmented and normalized
signals, segmented and non-normalized signals. Uses tombo pipeline for 
extraction and alignment of signal events. Requires a level file output from 
xenomorph.py preprocess. Outputs a new file that contains the raw signal around 
XNAs. 

Sample usage: 
python lib/xm_get_raw_xna.py input_fast5_dir ref_xfasta level_file output_file
python lib/xm_get_raw_xna.py input_fast5_dir ref_xfasta level_file output_file tombo_index_suffix

Cite us: 

H. Kawabe, C. Thomas, S. Hoshika, Myong-Jung Kim, Myong-Sang Kim, L. Miessner, N. Kaplan, J. M. Craig, 
J. Gundlach, A. Laszlo,  S. A. Benner, J. A. Marchand. "Enzymatic Synthesis and Nanopore Sequencing of 12-Letter Supernumerary DNA" 
Nature Communications. 14. (2023). DOI: 10.1038/s41467-023-42406-z 

Updated: 2/23/23
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
fast5_basedir = [sys.argv[1].replace('//','/')] #This is where the fast5 files are stored, required for going back to signal extract
reference_fasta = sys.argv[2].replace('//','/') #Reference file needs to match what was used for resquiggle 
level_input_file = sys.argv[3].replace('//','/') #'xenomorph_testing/PZ_1_MIN001PZ_CDxCD_levels.csv' #This has reads extracted 
output_folder = sys.argv[4].replace('//','/') #Currently not used, but output folder

#Load  level input file into pandas dataframe
level_file = pd.read_csv(level_input_file)

#Set up corrected group read if required 
if len(sys.argv)==6: 
    corr_group= 'RawGenomeCorrected_000'+sys.argv[5]
else: 
    corr_group= 'RawGenomeCorrected_000' 

#Determine XNAs in input level file
base_pairs = list(set(level_file['read_xna'])-set(['-']))

#Print status
print('Xenomorph Status - [Extract Raw] Performing raw level extraction on '+level_input_file)

#Create new output file 
output_column_names=['read_ID','xna_base','strand','sequence','raw_signal', 'norm_raw_signal', 'read_levels']
output_summary = pd.DataFrame(columns = output_column_names)



##############Additional function##############
def get_mean_signal(raw_signal, segs): 
    mean_sig = [] 
    for ii in range(0,len(segs)-1): 
        mean_sig.append(np.mean(raw_signal[segs[ii]:segs[ii+1]]))
    return(mean_sig)


with alive_bar(len(level_file), force_tty=True) as bar: 
    for i in range(0,len(level_file)): 
        bar()
        try: 


            #Get fast5 file name
            fast_f= level_file.iloc[i]['read_ID']
            
            #Get fast5 file path
            fast5_fn = os.path.normpath(fast5_basedir[0])+'/'+fast_f+'.fast5'

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


            #Get strand it mapped to 
            strand = level_file.iloc[i]['read_xna_strand'] 

            #Sets maps 1:1 with mapped_sequence. Search position in dp_results.genome_seq
            if strand == '+':
                st = level_file.iloc[i]['read_xna_position']-level_file.iloc[i]['read_signal_match_start']-xmer_boundary
                en = level_file.iloc[i]['read_xna_position']-level_file.iloc[i]['read_signal_match_start']+xmer_boundary+1
            if strand== '-':
                st = level_file.iloc[i]['read_xna_position']-level_file.iloc[i]['read_signal_match_start']+level_file.iloc[i]['read_signal_match_start']-xmer_boundary-2
                en = level_file.iloc[i]['read_xna_position']-level_file.iloc[i]['read_signal_match_start']+level_file.iloc[i]['read_signal_match_start']+xmer_boundary-1



            #########Outputs#########
            #Sequence of range
            ex_seq = mapped_sequence[st:en]

            #Segmentation of range (subtract minimum to re-center) 
            ex_seg = segs[st:en+1]-segs[st]

            #Get corresponding XNA base of the strand
            xna_base = level_file.iloc[i]['read_xna']

            #Ensure proper segmentation signal
            heptamer_padding = (2*xmer_boundary)+1

            if len(ex_seq)==7:
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


                #Normalized levels from initial resquiggle 
                read_levels = level_file.iloc[i]['read_levels']

                output_summary.loc[len(output_summary)] = [fast5_fn.split('/')[2], xna_base, strand, xna_ex_seq, ex_raw_sig,ex_nor_raw_sig, read_levels]

        except: 
            print('Xenomorph [Get Raw] - Error processing read.')

#Generate output file
output_summary.to_csv(output_folder)


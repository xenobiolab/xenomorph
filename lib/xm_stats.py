########################################################################
########################################################################
"""
xm_stats.py

Description: Calculates per-read, concensus, and global stats of a basecalled
level file (output from xenomorph.py morph command). Read quality filter and
concensus filters can be changed in xm_params.py. 

Update: 
-Split basecalling for each read


Title: Synthesis and Sequencing of a 12-letter supernumerary DNA

By: H. Kawabe, C. Thomas, S. Hoshika, Myong-Jung Kim, Myong-Sang Kim, L. Miessner, J. M. Craig, 
J. Gundlach, A. Laszlo,  S. A. Benner, J. A. Marchand

Updated: 10/17/23
"""
########################################################################
########################################################################


import pandas as pd
import numpy as np
import itertools
import seaborn as sns
import os
import sys
import math
from collections import Counter
from itertools import compress
from string import ascii_lowercase
import scipy.stats as stats
import matplotlib.pyplot as plt
import statistics

import warnings
warnings.filterwarnings("ignore")
import multiprocessing
import random
from alive_progress import alive_bar
from xm_params import * 


############################
#Xenomorph - xm_stats 
############################

if __name__ == '__main__': 
	#Quick function to find mode, even if there is more than 1
	def mmode(blist): 
		res = []
		c_list = Counter(blist) 
		temp = c_list.most_common(1)[0][1] 
		for ele in blist:
		    if blist.count(ele) == temp:
		        res.append(ele)
		res = list(set(res))
		concensus_freq = np.round(c_list[res[0]]/sum(c_list.values())*100)/100
		return res,concensus_freq


	##Load basecall output file 
	input_per_read_bc = sys.argv[1]
	output_basecalls = pd.read_csv(input_per_read_bc, sep=',')


	#Get list of all bases found in basecalling and list of all reads
	xnas = list(set(output_basecalls['read_xna'])-set(['-']))

	all_bases = xnas+['A','T','G','C']

	#Filter reads based on threshold setting. Modify in lib/xm_params.py. 
	print('Xenomorph Status - [Stats] Filtering reads with q score > '+str(qscore_filter))

	#Additional filters not yet added
	#print('Xenomorph Status - [Stats] Filtering reads with signal score < '+str(signal_filter))
	#print('Xenomorph Status - [Stats] Filtering reads with signal score < '+str(signal_filter))
	#output_basecalls=output_basecalls[output_basecalls['read_start']>read_start_filter]
	#output_basecalls=output_basecalls[output_basecalls['read_end']<read_end_filter]
	#print('Xenomorph Status - [Stats] Unique reads:'+str(len(list(set(output_basecalls['read_ID'])))))


	output_basecalls=output_basecalls[output_basecalls['read_q-score']>qscore_filter]
	output_basecalls=output_basecalls[output_basecalls['read_signal_match_score']<signal_filter]
	output_basecalls=output_basecalls[output_basecalls['xeno_basecall'] != '-'] 
	#print('Xenomorph Status - [Stats] Unique reads (filtered):'+str(len(list(set(output_basecalls['read_ID'])))))




	#Calculate correct base calls 
	for x in range(0,len(xnas)): 
		outx = output_basecalls[output_basecalls['read_xna']==xnas[x]]
		bc_count = 0 
		for xx in range(0,len(outx)): 
		    if outx.iloc[xx]['read_xna'] ==  outx.iloc[xx]['xeno_basecall']:
		        bc_count=bc_count+1
		print('Xenomorph Status - [Stats - Summary] Per-read recall of ['+xnas[x]+'] = '+str(round(bc_count/len(outx)*1000)/10)+'%   ('+str(bc_count)+'/'+str(len(outx))+')')

	#Get list of reads that pass filter
	ref_seqs = list(set(output_basecalls['reference_sequence']))


	#Set up output dataframe
	output_column_names=['reference_sequence','n_reads','mean_q-score','mean_signal_match_score','reference_locus',
	'strand','xna_sequence','xna', 'xna_freq','concensus_basecall', 'concensus_basecall_freq', 'xna_is_conensus', 
	'model_sigma', 'model_mean', 'read_xna_pos', 'model_file']

	output_summary = pd.DataFrame(columns = output_column_names)

	#numbers to calculate consensus count
	consensus_count=0
	total_count=0

	#numbers to calculate consensus count
	consensus_count=0
	total_count= {}
	global_count = {}
	
	for i in range(0,len(ref_seqs)): 
		ref_ii = output_basecalls[output_basecalls['reference_sequence']==ref_seqs[i]]

		#Reference sequence 
		ref_i_id = ref_seqs[i] 

		#Position
		ref_i_p = list(set(ref_ii['read_xna_position']))

		for p in range(0,len(ref_i_p)): 

		    ref_i = ref_ii[ref_ii['read_xna_position']==ref_i_p[p]]
		    #Strands sequenced (+,-)
		    ref_s = list(set(ref_i['read_xna_strand']))


		    for j in range(0,len(ref_s)): 

		        #Strand
		        ref_strand=ref_s[j]

		        #Reads for strand 
		        ref_i_s = ref_i[ref_i['read_xna_strand']==ref_strand]

		        #Number of reads sequenced for this strand and reference
		        ref_i_n = len(ref_i_s) 

		        #Read reference locus
		        ref_i_loc = ref_i_s.iloc[0]['read_reference_locus']

		        #Read reference sequence (= to locus if (+) strand, reverse complement of locus for (-) strand)
		        ref_i_seq = ref_i_s.iloc[0]['read_xna_sequence']

		        #xna pos
		        ref_i_xnapos = ref_i_p[p]

		        #Xenobase prior 
		        ref_i_x = ref_i_s.iloc[0]['read_xna']

		        #XNA prior frequency 
		        ref_i_x_freq = ref_i_s['xeno_basecall'].to_list().count(ref_i_x)
		        ref_i_x_freq = np.round(ref_i_x_freq/len(ref_i_s['xeno_basecall'].to_list())*100)/100

		        #Basecalls made
		        ref_i_bc = ref_i_s['xeno_basecall']

		        #Concensus and fraction of reads called to censensus 
		        ref_i_concensus, ref_i_concensus_freq = mmode(ref_i_s['xeno_basecall'].to_list())

		        #Check if XNA is one of the conensus base calls
		        ref_i_x_is_concensus = ref_i_concensus_freq == ref_i_x_freq

		        #Average q-score
		        ref_i_q = np.mean(ref_i_s['read_q-score'].to_list())

		        #Average signal score
		        ref_i_sc = np.mean(ref_i_s['read_signal_match_score'].to_list())

		        #model sigma used 
		        ref_i_sig = ref_i_s.iloc[0]['model_sigma']

		        #model mean used
		        ref_i_mean = ref_i_s.iloc[0]['model_mean']

		        #model used
		        ref_i_model = ref_i_s.iloc[0]['model_file']

		        ref_out = [ref_i_id, ref_i_n, ref_i_q, ref_i_sc,ref_i_loc, ref_strand, ref_i_seq, ref_i_x, ref_i_x_freq, ref_i_concensus, ref_i_concensus_freq, ref_i_x_is_concensus, ref_i_sig, ref_i_mean, ref_i_xnapos, ref_i_model]

		        output_summary.loc[len(output_summary)]=ref_out

		                
		        if ref_i_n >=concensus_stat_filter: 

		            if ref_i_x not in global_count:
		                global_count[ref_i_x]=0
		                
		            try: 
		                total_count[ref_i_x] = total_count[ref_i_x]+1
		            except: 
		                total_count[ref_i_x] = 1
		            consensus_count +=1

		            if ref_i_x_is_concensus ==True:
		                try: 
		                    global_count[ref_i_x] = global_count[ref_i_x]+1
		                except: 
		                    global_count[ref_i_x] = 1
		                consensus_count +=1

	outfn = input_per_read_bc.replace('.csv','_per-read_concensus')+'.csv'
	output_summary.to_csv(outfn, index=False)

	#Print output summary of global morph result
	for xna_key in global_count: 
		print("Xenomorph Status - [Stats - Summary] Per-read consensus (n > "+str(concensus_stat_filter)+") of ["+xna_key+"] = "+str(round((global_count[xna_key]/total_count[xna_key])*100,1))+ "%     "+str(global_count[xna_key]) + "/" + str(total_count[xna_key]))


	print('Xenomorph Status - [Stats] Saving global summary of stats to '+outfn)




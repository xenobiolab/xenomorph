########################################################################
########################################################################
"""
xenomorph.py 

Description: Xenomorph - XNA sequencing pipeline. Manages scripts for 
integrating nanopore preprocessing with nanopore xenonucleotide basecalling. 

Usage: python xenomorph.py -h for help 

Title: Synthesis and Sequencing of a 12-letter supernumerary DNA

By: H. Kawabe, C. Thomas, S. Hoshika, Myong-Jung Kim, Myong-Sang Kim, L. Miessner, J. M. Craig, 
J. Gundlach, A. Laszlo,  S. A. Benner, J. A. Marchand

Updated: 2/14/23
"""
########################################################################
########################################################################

import argparse, textwrap
import os 
import sys 
from lib.xm_tools import *
from lib.xm_params import *


parser = argparse.ArgumentParser(
        usage='"python xenomorph.py  [-h] {preprocess, discover, morph, de-novo, fasta2x}',
        formatter_class=argparse.RawTextHelpFormatter,
	description= textwrap.dedent('''



                                                          _     
                                                         | |    
__  __  ___  _ __    ___   _ __ ___    ___   _ __  _ __  | |__  
\ \/ / / _ \| '_ \  / _ \ | '_ ` _ \  / _ \ | '__|| '_ \ | '_ \ 
 >  < |  __/| | | || (_) || | | | | || (_) || |   | |_) || | | |
/_/\_\ \___||_| |_| \___/ |_| |_| |_| \___/ |_|   | .__/ |_| |_|
                                                  | |           
                                                  |_|           


********** Xenomorph : An XNA sequencing toolkit *********

Xenomorph is a suite of tools used for nanopore sequencing of alternative basepairs (XNAs). This toolkit incorporates ONT-workflows to preprocess fast5 files and extract signal levels for kmers in a sequence. Models parameterized on XNA basepairs can then be used to test if signal levels match XNA pairs. Xenomorph relies on kmer models that were parameterized using libraries of XNA-containing DNA. The general pipeline consists of two steps: 1) preprocessing fast5 reads and extracting level information and 2) basecalling using a selected kmer model. 


Xenomorph command groups (additional help available within each command group):
	preprocess	[Preprocess] fast5 reads with reference fasta containing XNAs for 'morph' or 'de-novo'
	morph		[Basecall] Use kmer levels identified with preprocess to basecall XNA position based on alternative hypothesis testing (per read)
	stats		[Utility] Calculate global concensus basecalls from per-read output file and generate summary output
	fasta2x		[Utility] Converts fasta file with XNAs in sequence (e.g. BSPZKXJV) to fasta file with XNA positional information in header
	models		[Utility] View summary of active and inactive kmer models that can be used for basecalling or activate models
         '''))


subparsers = parser.add_subparsers(dest='subparsers')



#Preprocess
parser_preprocess = subparsers.add_parser('preprocess', help='[-w working_dir] [-f fast5_dir] [-r reference_fasta]')
parser_preprocess.add_argument('-f',metavar ='[fast5_dir]', type=str,required = True, help='Input directory containing multi-fast5 folders.')
parser_preprocess.add_argument('-w',metavar = '[working_dir]', type=str,required = True, help='Working directory for storing analysis pipeline temp files and outputs.')
parser_preprocess.add_argument('-r',metavar = '[reference_fasta]', type=str, required = True, help='Fasta (.fa, .fasta) file of sequence or sequences with XNAs (e.g. BSPZKXJV) in sequence.')
parser_preprocess.add_argument('-t',metavar = '[tombo_index_id]', type=str, required = True, help='Store resquiggle results with the following ID (default 000)')
parser_preprocess.add_argument('-o',metavar = '[output_file]', type=str, required = False, help='Specify path of output file (.csv) for level extraction (optional; default = working_dir/level_output_summary.csv)')
parser_preprocess.add_argument('-b',action = 'store_true', help='Use guppy to add basecall reads to file (optional).')
parser_preprocess.add_argument('-x',action = 'store_true', help='Overwrite all temp analysis files and force re-analysis (optional).')



#Morph 
parser_morph = subparsers.add_parser('morph', help='[-l level_summary_file] [-m model_input] [-o output_summary_file] [-g global_basecalling]')
parser_morph.add_argument('-l',metavar = '[level_summary_file]', type=str,required = True, help='Working directory for storing analysis pipeline temp files and outputs.')
parser_morph.add_argument('-m',metavar = '[model_file]', type=str, required = True, help='Single letter abbreviations (e.g. ATGCBS) or specify kmer model file (.csv).')
parser_morph.add_argument('-o',metavar = '[output_file]', type=str, help='Output file name (optional; default output_summary_basecalls.csv)')
parser_morph.add_argument('-g',action = 'store_true', help='Perform global basecalling (across all reads) instead of per-read base calling.')
parser_morph.add_argument('-n',action = 'store_true', help='Perform global null basecalling using reads without XNAs.')


#Stats
parser_stats = subparsers.add_parser('stats', help='[-i input_bc_file]')
parser_stats.add_argument('-i',metavar ='[input_bc_file]', type=str,required = True, help='Calculate global concensus basecall from per-read output file (generated by morph).')


#Fasta2x
parser_fasta2x = subparsers.add_parser('fasta2x', help='[-i fasta_input] [-o xfasta_output]')
parser_fasta2x.add_argument('fasta_input',metavar ='[fasta_input]', type=str, help='Input fasta file containing XNAs in sequence.')
parser_fasta2x.add_argument('xfasta_output', metavar ='[xfasta_output]', type=str, help='Output fasta file with XNA positional information in header and gapped sequences.')


#Models
parser_models = subparsers.add_parser('models', help='[-s active/inactive/all] [-a base_abbreviation]')
parser_models.add_argument('-a',metavar ='[abbreviation]', type=str,required = False, help='Activate model params for XNA with the given abbreviation (e.g. -a Bn to activate Bn, deactivate Bc).')
parser_models.add_argument('-s',metavar = '[active/inactive/all]', type=str, required = False, help='Show model descriptions of active, inactive, or all kmer models.')



args = parser.parse_args()


args.subparsers


#Print help if no arguments provided
if args.subparsers==None: 
	parser.print_help()
	sys.exit(0)

#Parse arguments 
if args.subparsers == 'fasta2x' :
	os.system('python lib/xm_fasta2x_rc.py '+args.fasta_input+' '+args.xfasta_output)


if args.subparsers == 'morph': 

	if os.path.exists(args.l)==False:
		print("Xenomorph Status - [Error] Level-extracted input file not found. Check path or run 'xenomorph preprocess' to generate")
		sys.exit()



	if sum(1 for line in open(args.l)) <=1:
		print("Xenomorph Status - [Error] No reads found in level extraction input file. Ensure 'xenomorph preprocess' ran successfully")
		sys.exit()


	level_summary_file = args.l

	#load model fn
	if '.csv' in args.m: 
		model_fn = args.m
	elif all(elem in 'ATGCBSPZJVKX' for elem in args.m):
		print('Xenomorph Status - [Morph] Compiling model from input model base abbreviations.')

		#Parse model, get paths of each model file component
		mf = parse_model_files(args.m,True)

		#Compile model from file paths into a master kmer model
		km = compile_model(mf, args.m) 

		#Save model 
		model_fn=os.path.normpath(args.l).replace(os.path.basename(args.l),'')+'active_kmer_model.csv'
		km.to_csv(model_fn, index=False)

	if args.o:
		output_fn=args.o
	else: 
		output_fn=os.path.normpath(args.l).replace(os.path.basename(args.l),'')+'/output_summary_basecalled.csv'


	if args.g: 
		if args.n:
			print("Xenomorph Status - [Morph Global Null] Testing null hypothesis of XNA basecalls by grouping reads that map to the same region (global).")
			cmd =  'python lib/global_null_morph.py '+model_fn+' '+level_summary_file+' '+output_fn
			os.system(cmd)
		else:
			print("Xenomorph Status - [Global Morph] Testing alternative hypothesis of XNA basecalls by grouping reads that map to the same region (global).")
			cmd =  'python lib/global_morph.py '+model_fn+' '+level_summary_file+' '+output_fn
			os.system(cmd)
	elif args.n: 
		print("Xenomorph Status - [Morph Null] Testing null hypothesis of XNA basecalls by grouping reads that map to the same region (global).")
		cmd =  'python lib/null_morph.py '+model_fn+' '+level_summary_file+' '+output_fn
		os.system(cmd)
	else: 
		print("Xenomorph Status - [Morph] Testing alternative hypothesis of XNA basecalls on individual reads (per-read).")
		cmd =  'python lib/morph.py '+model_fn+' '+level_summary_file+' '+output_fn
		os.system(cmd)


elif args.subparsers == 'preprocess': 

	working_dir = os.path.normpath(args.w)

	#Create a working directory if it does not existW
	CHECK_FOLDER = os.path.isdir(args.w)
	if not CHECK_FOLDER:
		os.makedirs(working_dir)
		print("Xenomorph Status - [Preprocess] Creating output working directory : "+working_dir)




	#Convert reference fasta to xfasta
	if os.path.exists(args.r)==False:
		print("Xenomorph Status - [Error] Fasta file input not found.")
		sys.exit()
	if not check_xfasta_format(args.r, standard_bases): 
	    xfasta=os.path.normpath(args.w)+'/xfasta_'+os.path.basename(args.r)
	    if not os.path.exists(xfasta) or args.x:
		    print("Xenomorph Status - [Preprocess] Performing fasta to xfasta conversion on input reference fasta. Saving output "+xfasta+'.')
		    cmd =  'python lib/xm_fasta2x_rc.py '+os.path.normpath(args.r)+' '+xfasta
		    os.system(cmd) 
	    else: 
		    print("Xenomorph Status - [Preprocess] xfasta file found in working directory. Skipping xfasta conversion.")
	else:
		xfasta=args.r
		print("Xenomorph Status - [Preprocess] Input reference fasta already in xfasta format. Skipping xfasta conversion.")

	#check for creation of a reverse fasta file - required for some XNAs for segmentation 
	if os.path.exists(xfasta[0:xfasta.find('.fa')]+'_rc.fa')==True:
		print("Xenomorph Status - [Preprocess] xfasta reverse complement created for handling special bases.")
		xfasta_rc = xfasta[0:xfasta.find('.fa')]+'_rc.fa'

		###Multi to single and merge single folders
	fast5_sing_dir=os.path.normpath(working_dir)+'/'+os.path.basename(args.f)+'_single'
	CHECK_FAST5_SINGLE = os.path.isdir(fast5_sing_dir)
	if not CHECK_FAST5_SINGLE:
		print("Xenomorph Status - [Preprocess] Converting multi fast5 to single fast5 using ONT tools.")
		cmd =  'python xombo.py split '+os.path.normpath(args.f)+' '+fast5_sing_dir
		os.system(cmd) 
	else: 
		print("Xenomorph Status - [Preprocess] Single fast5 directory found in working directory. Skipping multi-to-single conversion.")


		#Check that fast5 multi-to-single worked properly. If not, exit with error message to check paths. 
	CHECK_FAST5_SINGLE_N = sum('.fast5' in s for s in os.listdir(fast5_sing_dir))
	if CHECK_FAST5_SINGLE_N > 0: 
		print("Xenomorph Status - [Preprocess] Found "+str(CHECK_FAST5_SINGLE_N)+' single fast5 files in directory')
	else: 
		print("Xenomorph Status - [Error] Single fast5 files not found. Check paths to ensure multi-to-single fast5 conversion is working properly.")
		sys.exit()

	###Rebasecall - required for resquiggling and assignment
	if args.b == True: 
		fastq_dir=os.path.normpath(args.w)+'/'+os.path.basename(args.f)+'_fastq'
		print("Xenomorph Status - [Preprocess] Using guppy to basecall single fast5 files.")
		if use_reference_to_basecall==True: 
			cmd =  'python xombo.py basecall '+fast5_sing_dir+' '+fastq_dir+' -r '+xfasta
		else: 
			cmd =  'python xombo.py basecall '+fast5_sing_dir+' '+fastq_dir
		os.system(cmd) 

		###Assign reads
		print("Xenomorph Status - [Preprocess] Using tombo preprocess to assign fastq basecalls to fast5 files.")
		cmd = 'python xombo.py preprocess '+fast5_sing_dir+' '+fastq_dir
		os.system(cmd) 
	else: 
		print("Xenomorph Status - [Preprocess] No basecall flag (-b) entered. Proceeding with resquiggle.")



	#Set up index group from -t flag
	if args.t: 
		index_group = args.t
		rc_index_group = args.t+'r'
	else: 
		index_group = ''
		rc_index_group = 'r'

	#Resquiggle if needed
	tombo_index= os.path.normpath(args.w)+'/.'+os.path.basename(fast5_sing_dir)+'.RawGenomeCorrected_000'+index_group+'.tombo.index'
	if not os.path.exists(tombo_index) or args.x:
		print("Xenomorph Status - [Preprocess] Using tombo to resquiggle fast5 reads using xfasta reference.")
		cmd = 'python xombo.py resquiggle '+fast5_sing_dir+' '+os.path.normpath(xfasta)+' -g '+index_group
		os.system(cmd) 
		if os.path.exists(xfasta[0:xfasta.find('.fa')]+'_rc.fa')==True:
			print("Xenomorph Status - [Preprocess] Using tombo to resquiggle fast5 reads using rc-xfasta reference.")
			cmd = 'python xombo.py resquiggle '+fast5_sing_dir+' '+os.path.normpath(xfasta_rc)+' -g '+rc_index_group
			os.system(cmd) 
	else: 
		print("Xenomorph Status - [Preprocess] Resquiggle file found. Skipping and proceeding to level extraction. Use (-x) flag to force resquiggle with new xfasta file.")



	if args.o:
		level_output_fn=os.path.normpath(args.w)+'/'+args.o

	else: 
		level_output_fn= os.path.normpath(args.w)+'/output_levels_summary.csv'


	print("Xenomorph Status - [Preprocess] Performing level extraction surrounding XNA locations.")
	cmd = 'python lib/xm_get_levels.py '+fast5_sing_dir+' '+os.path.normpath(xfasta)+' '+level_output_fn+' '+index_group
	os.system(cmd) 

	if os.path.exists(xfasta[0:xfasta.find('.fa')]+'_rc.fa')==True:
		print("Xenomorph Status - [Preprocess] Performing level extraction on surrounding XNA locations with reverse set.")
		cmd = 'python lib/xm_get_levels.py '+fast5_sing_dir+' '+os.path.normpath(xfasta_rc)+' '+level_output_fn+' '+rc_index_group
		os.system(cmd) 
	print("Xenomorph Status - [Preprocess] Saving output level file to "+level_output_fn)
	print("Xenomorph Status - [Preprocess Complete] Use 'xenomorph.py morph' for alternative hypothesis testing on XNA positions.")



#Model handling 
elif args.subparsers == 'models': 

	if args.a: 
		try: 

			activate_model(args.a)
			model_summary('all')
			print("Xenomorph Status - [Models] "+args.a+" successfully set as active model within group. Other models sharing same single letter code inactivated.")
		except: 
			print("Xenomorph Status - [Error] Error activating model. Check proper abbreivation was used, or manually activate by editing models/config_models.csv.")

	if args.s: 
		model_summary(args.s)


elif args.subparsers == 'stats': 

	try: 
		cmd = 'python lib/xm_stats.py '+args.i
		os.system(cmd)
	except: 
		print("Xenomorph Status - [Error] Stats failed. Ensure file path to basecalled output.csv is correct, or run xenomorph morph to generate new basecall file.")









##############################
##Xombo - An XNA toolkit  ###
##############################

import argparse, textwrap
import os 
import sys 
from lib.xm_tools import *
from lib.xr_tools import *
from lib.xm_params import *
from lib.xr_params import *

parser = argparse.ArgumentParser(
        usage='"python xenomorph.py  [-h] {preprocess, morph, extract, models, fasta2x}',
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
	extract		[Utility] Extracts raw signal in region associated with XNA bases. 
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
#parser_preprocess.add_argument('-t',metavar = '[tombo_index_id]', type=str, required = True, help='Store resquiggle results with the following ID (default 000)')
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

#Extract
parser_extract = subparsers.add_parser('extract', help='[-w working_dir] [-f fast5_dir] [-r reference_fasta] [-l level_summary_file] [-o output_raw_file]')
parser_extract.add_argument('-w',metavar = '[working_dir]', type=str,required = True, help='Working directory for storing analysis pipeline temp files and outputs.')
parser_extract.add_argument('-f',metavar ='[fast5_dir]', type=str,required = True, help='Input directory containing multi-fast5 folders.')
parser_extract.add_argument('-r',metavar = '[reference_fasta]', type=str, required = True, help='Fasta (.fa, .fasta) file of sequence or sequences with XNAs (e.g. BSPZKXJV) in sequence.')
parser_extract.add_argument('-l',metavar = '[level_summary_file]', type=str,required = True, help='Output file (.csv) from preprocess that contains normalized levels of each read.')
parser_extract.add_argument('-o',metavar = '[output_file]', type=str, required = True, help='Output file name for raw reads (optional; default output_levels_raw.csv)')



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



elif args.subparsers == 'extract': 
	print("Xenomorph Status - [Extract] Extracting raw signal in region associated with XNA level measurements.")
	working_dir = os.path.normpath(args.w)
	CHECK_FOLDER = os.path.isdir(args.w)
	if not CHECK_FOLDER:
		print("Xenomorph Status - [Error] Working directory "+args.w+" not found.")
		print("Xenomorph Status - [Error] Run xenomorph.py preprocess first to generate directories.")


	if os.path.exists(args.l)==False:
		print("Xenomorph Status - [Error] Could not locate file: "+args.l)
		print("Xenomorph Status - [Error] Level-extracted input file not found. Check path or run 'xenomorph preprocess' to generate.")
		sys.exit()

	level_summary_file = args.l


	if not check_xfasta_format(args.r, standard_bases): 
		xfasta=os.path.normpath(args.w)+'/xfasta_'+os.path.basename(args.r)

	if args.o:
		output_fn=args.o


	fast5_sing_dir=os.path.normpath(working_dir)+'/'+os.path.basename(args.f)+'_single'

	xfasta=os.path.normpath(args.w)+'/xfasta_'+os.path.basename(args.r)


	cmd = 'python lib/xm_get_raw_signal.py '+fast5_sing_dir+' '+xfasta+' '+level_summary_file+' '+output_fn
	os.system(cmd)



elif args.subparsers == 'preprocess': 

    #Set working directory path 
	working_dir = os.path.normpath(args.w)
	
	#Create a working directory if it does not exist
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
		    print("Xenomorph Status - [Preprocess] Using fasta2x_rc.")
		    cmd =  'python lib/xm_fasta2x_rc.py '+os.path.normpath(args.r)+' '+xfasta
		    os.system(cmd) 
	    else: 
		    print("Xenomorph Status - [Preprocess] xfasta file found in working directory. Skipping xfasta conversion.")
	else:
		xfasta=args.r
		print("Xenomorph Status - [Preprocess] Input reference fasta already in xfasta format. Skipping xfasta conversion.")

	#Check for creation of a reverse fasta file - required for some XNAs for segmentation 
	if os.path.exists(xfasta[0:xfasta.find('.fa')]+'_rc.fa')==True:
		print("Xenomorph Status - [Preprocess] xfasta reverse complement created for handling special bases.")
		xfasta_rc = xfasta[0:xfasta.find('.fa')]+'_rc.fa'
		xfasta_rc_dir = os.path.join(args.w+'/'+os.path.basename(xfasta_rc))





    ###### Remora preprocessing integration 
    #Convert fast5 to pod5 format 
	pod5_dir = os.path.normpath(args.w)+'/'+os.path.basename(args.f)+'_pod5'
	pod5_dir = check_make_dir(pod5_dir)
	check_pod5_dir = os.path.join(pod5_dir+'/'+os.path.basename(args.f)+'.pod5')
	if os.path.isfile(check_pod5_dir)==False: 
		cod5_to_fast5(get_fast5_subdir(args.f), os.path.join(pod5_dir + '/'+os.path.basename(args.f))+'.pod5')
	else: 
		print('Xenomorph Status - [Preprocess] POD5 file for modified base found. Skipping POD5 coversion')

    #Basecall pod5 files 
	fastq_dir = os.path.normpath(args.w)+'/'+os.path.basename(args.f)+'_fastq'
	if os.path.isdir(fastq_dir)==False or basecall_pod == True: 
		fastq_dir = check_make_dir(fastq_dir)
		cmd=os.path.expanduser(basecaller_path)+' -i '+pod5_dir+' -s '+fastq_dir+' -c '+guppy_config_file+' -x auto --bam_out --index --moves_out -a '+os.path.join(args.w+'/'+os.path.basename(xfasta))
		cmd=os.path.expanduser(basecaller_path)+' -i '+pod5_dir+' -s '+fastq_dir+' -c '+guppy_config_file+' --min_qscore 9 --align_type full -x auto --bam_out --index --moves_out -a '+os.path.join(args.w+'/'+os.path.basename(xfasta))
		os.system(cmd)
	else:
		print('Xenomorph Status - [Preprocess] Skipping POD5 basecalling for modified bases.')

		###Merge Bam files 
	bam_dir = os.path.normpath(args.w)+'/'+os.path.basename(args.f)+'_bam'
	check_bam_dir=bam_dir+'.bam'
	if os.path.isfile(bam_dir+'.bam') == False or regenerate_bam == True: 

		# commented out lines below added to merge both pass and fail into one master bam file

		#merge pass bam files
		#cmd = 'samtools merge '+os.path.join(bam_dir+'pass.bam'+' '+os.path.join(fastq_dir,'pass/*.bam -f'))
		#os.system(cmd)

		#merge fail bam files
		#cmd = 'samtools merge '+os.path.join(bam_dir+'fail.bam'+' '+os.path.join(fastq_dir,'fail/*.bam -f'))
		#os.system(cmd)

		#merge all bam files
		#cmd = 'samtools merge '+os.path.join(bam_dir+'.bam'+' '+os.path.join(args.w + '/*.bam -f'))


		cmd = 'samtools merge '+os.path.join(bam_dir+'.bam'+' '+os.path.join(fastq_dir,'pass/*.bam -f'))
		print('Xenomorph Status - [Preprocess] Merging modified BAM files.')
		os.system(cmd)
	else:
		print('Xenomorph Status - [Preprocess] Skipping merging modified BAM files.')


		###Bed file generation 
	xfasta_dir = os.path.join(args.w+'/'+os.path.basename(xfasta))
	if os.stat(xfasta_dir).st_size == 0: 
		print('Xenomorph Status - [Preprocess] Empty xfasta file generated. Check that XNA bases were present in sequence of input fasta file.')
		sys.exit()
	print('Xenomorph Status - [Preprocess] Generating bed file for modified base.')
	cmd = 'python lib/xr_xfasta2bed.py '+os.path.join(xfasta_dir)+' '+os.path.join(args.w,mod_base+'.bed ' +mod_base+' '+mod_base)
	bed_dir = os.path.join(args.w,mod_base+'.bed')
	os.system(cmd)


	if args.o:
		level_output_fn=os.path.normpath(args.w)+'/'+args.o

	else: 
		level_output_fn= os.path.normpath(args.w)+'/output_levels_summary.csv'

	#Remora segmentation will require rescaling (e.g., MAD to SD) for comparison with Tombo models
	if manual_rescale_override == False:
		#Make rescale directory if it does not exist.
		CHECK_FOLDER = os.path.isdir(working_dir+'/rescale')
		if not CHECK_FOLDER:
			os.makedirs(working_dir+'/rescale')
		#Run level extract on raw data 
		if force_extract_position == True: 
			##For fixed position extract 
			cmd = 'python lib/xr_get_levels_pos.py '+check_pod5_dir+' '+ check_bam_dir + ' '+ bed_dir + ' '  +xfasta_dir+' '+working_dir+'/rescale/rescale_raw_levels.csv rescale'
			os.system(cmd) 
		else:
			#For standard pipeline
			cmd = 'python lib/xr_get_levels.py '+check_pod5_dir+' '+ check_bam_dir + ' '+ bed_dir + ' '  +xfasta_dir+' '+working_dir+'/rescale/rescale_raw_levels.csv rescale'
			os.system(cmd) 
		
		#Covert raw level file to kmer extract csv
		cmd = 'python lib/xm_extract_levels.py '+working_dir+'/rescale/rescale_raw_levels.csv'
		os.system(cmd) 

		cmd = 'python lib/parse_kmer.py '+working_dir+'/rescale/rescale_raw_kmers.csv'+' '+working_dir+'/rescale/rescale_raw_kmer_model.csv ATGC'
		os.system(cmd) 

		cmd = 'python lib/xr_kmer_rescale.py '+working_dir+'/rescale/rescale_raw_kmer_model.csv '+working_dir+'/rescale/rescale_kmer_comparison.pdf'
		os.system(cmd) 

	#After rescaling paramters are calculated, continue with pipeline
	if force_extract_position == True: 
		#For fixed position extract
		print("Xenomorph Status - [Warning] Overriding XNA kmer extraction position. Forced positional extract is set to true (not typical).")
		print("Xenomorph Status - [Warning] Change this setting in lib/xr_params.py")
		print("Xenomorph Status - [Preprocess] Performing level extraction surrounding the specified position location.")
		cmd = 'python lib/xr_get_levels_pos.py '+check_pod5_dir+' '+ check_bam_dir + ' '+ bed_dir + ' '  +xfasta_dir+' '+level_output_fn
		os.system(cmd) 

	else: 
		#For standard pipeline
		print("Xenomorph Status - [Preprocess] Performing level extraction surrounding XNA locations.")
		cmd = 'python lib/xr_get_levels.py '+check_pod5_dir+' '+ check_bam_dir + ' '+ bed_dir + ' '  +xfasta_dir+' '+level_output_fn
		os.system(cmd) 

	if os.path.exists(xfasta[0:xfasta.find('.fa')]+'_rc.fa')==True:
		print("Xenomorph Status - [Preprocess] Performing level extraction on surrounding XNA locations with reverse set.")
		cmd = 'python lib/xr_get_levels.py '+check_pod5_dir+' '+ check_bam_dir +' '+ bed_dir + ' ' +xfasta_rc_dir+' '+level_output_fn
		os.system(cmd) 
	print("Xenomorph Status - [Preprocess] Saving output level file to "+os.path.normpath(level_output_fn))
	print("Xenomorph Status - [Preprocess Complete] Use 'xenomorph.py morph' for alternative hypothesis testing on XNA positions.")
    #####End Remora Section





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









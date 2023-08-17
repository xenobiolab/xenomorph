########################################################################
########################################################################
"""
xemora.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. A. Marchand

Updated: 3/2/23
"""
########################################################################
########################################################################

import argparse, textwrap
import os 
import sys 
from lib.xr_tools import *
from lib.xr_params import *


parser = argparse.ArgumentParser(
        usage='"python xemora.py  [-h] {train, basecall}',
        formatter_class=argparse.RawTextHelpFormatter,
	description= textwrap.dedent('''



______________________________________________________________________________
______________________________________________________________________________
_/\/\____/\/\_________________________________________________________________
___/\/\/\/\______/\/\/\____/\/\/\__/\/\______/\/\/\____/\/\__/\/\__/\/\/\_____
_____/\/\______/\/\/\/\/\__/\/\/\/\/\/\/\__/\/\__/\/\__/\/\/\/\________/\/\___
___/\/\/\/\____/\/\________/\/\__/\__/\/\__/\/\__/\/\__/\/\________/\/\/\/\___
_/\/\____/\/\____/\/\/\/\__/\/\______/\/\____/\/\/\____/\/\________/\/\/\/\/\_
______________________________________________________________________________
______________________________________________________________________________


********** Xemora : An XNA sequencing neural network trainer *********

Xemora is a tool pipeline used for nanopore sequencing of alternative basepairs (XNAs) that latches onto Remora 2.0 (ONT). This toolkit incorporates ONT-workflows to preprocess fast5, pod5, fastq, bam, and bed files for training a remora model. XNA sequence handling is done using the xFASTA format. Therefore, FASTA reference inputs should contain XNA bases in sequence lines. Models for basecalling can be trained on specific set of sequences (rather than entire sequence space). For optimal implementation, Xemora should be trained on at least two datasets, with and without the XNA substitutions. Xemora models can be exported and used as remora models for guppy, dorado, or bonito basecalling. Alternatively, Xemora can handle reference-based XNA identification directly. The general pipeline consists of two steps: 1) Training xemora on a reference set of reads with and without XNAs 2) basecalling using the trained model. 


Xemora command groups (additional help available within each command group):
	train		[Train] a xemora model on a set of input fast5 reads, localized to reference fasta containing XNAs. 
	basecall	[Basecall] a fast5 reads around XNA using a previously trained xemora model. 
         '''))


subparsers = parser.add_subparsers(dest='subparsers')



#Train
parser_train = subparsers.add_parser('train', help='[-w output_dir] [-f xna_fast5_dir dna_fast5_dir] [-r xna_ref_fasta dna_ref_fasta]')
parser_train.add_argument('-w',metavar = '[output_dir]', type=str,required = True, help='Path to output directory for storing output intermediates, temp files, and models.')
parser_train.add_argument('-f',metavar ='[xna_fast5_dir dna_fast5_dir]', nargs=2, type=str,required = True, help='Path to input directories containing multi-fast5 folders of XNA-containing sequence and DNA-containing sequence.')
parser_train.add_argument('-r',metavar = '[xna_ref_fasta dna_ref_fasta]', type=str, nargs=2, required = True, help='Path to FASTA (.fa, .fasta) file of sequence or sequences with XNAs (e.g. BSPZKXJV) in the position focused for training. XNA should be in both xna_reference_fasta and dna_reference_fasta. ')


#Basecall 
parser_basecall = subparsers.add_parser('basecall', help='[-w output_dir] [-f xna_fast5_dir] [-r xna_ref_fasta] [-m model_file]')
parser_basecall.add_argument('-w',metavar = '[output_dir]', type=str,required = True, help='Path to output directory for storing output intermediates and basecall results.')
parser_basecall.add_argument('-f',metavar ='[xna_fast5_dir]', type=str,required = True, help='Path to input directories containing multi-fast5 folders of XNA-containing sequence.')
parser_basecall.add_argument('-r',metavar = '[xna_ref_fasta]', type=str, required = True, help='Path to FASTA (.fa, .fasta) file of sequence or sequences with XNAs (e.g. BSPZKXJV). Should be same sequence context as xemora model training.')
parser_basecall.add_argument('-m',metavar = '[model_file]', type=str, required = True, help='Path to trained xemora model (e.g. model/model_best.pt)')



args = parser.parse_args()
args.subparsers
exit_flag = False

#Print help if no arguments provided
if args.subparsers==None: 
    parser.print_help()
    sys.exit(0)



if args.subparsers == 'train': 
    if os.path.exists(os.path.expanduser(args.f[0]))==False or os.path.exists(os.path.expanduser(args.f[1]))==False:
        print('Xemora [ERROR] - Training requires two input fast5 datasets. At least one of the two file paths was not valid. Check to ensure path to fast5 directory is correct.')
        exit_flag = True 

    if os.path.exists(os.path.expanduser(args.r[0]))==False or os.path.exists(os.path.expanduser(args.r[1]))==False:
        print('Xemora [ERROR] - Training requires two input reference fasta files: One to designate XNA-containing reads and the other with reads that contain a standard base substitution. At least one of the two file paths was not valid. Check to ensure path to fast5 directory is correct.')
        exit_flag = True 


    if exit_flag == False: 
        cmd = 'python lib/xr_train.py '+args.w+' '+args.f[0]+' '+args.r[0]+' '+args.f[1]+' '+args.r[1]
        os.system(cmd)
    else: 
        print('Xemora [ERROR] - At least one file path not properly set. Xemora exiting.')
        sys.exit()



#Guppy paths and stuff - add that to xr-parms
if args.subparsers == 'basecall': 
    if os.path.exists(os.path.expanduser(args.f))==False:
        print('Xemora [ERROR] - Input Fast5 directory path not found or is not valid. Check to ensure path to fast5 directory is correct.')
        exit_flag = True 

    if os.path.exists(os.path.expanduser(args.r))==False:
        print('Xemora [ERROR] - Xemora basecalling requires identifying locations of possible XNAs by giving it a reference fasta file with XNA bases (e.g., BSPZXKJV) in the sequence fields.')
        exit_flag = True 

    if os.path.exists(os.path.expanduser(args.m))==False:
        print('Xemora [ERROR] - Valid model file (e.g. model/best_model.pt) not provided. Run "xemora train" to generate a model for basecalling.')
        exit_flag = True 


    if exit_flag == False: 
        cmd = 'python lib/xr_basecall.py '+args.w+' '+args.f+' '+args.r+' '+args.m
        os.system(cmd)
    else: 
        print('Xemora [ERROR] - At least one file path not properly set. Xemora basecaller exiting.')
        sys.exit()





########################################################################
########################################################################
"""
xombo.py 

Description: Xombo - An XNA toolkit. Manages scripts for integrating 
various ONT processing pipelines required for signal extraction and
signal alignment of raw nanopore reads. 

Cite us: 

H. Kawabe, C. Thomas, S. Hoshika, Myong-Jung Kim, Myong-Sang Kim, L. Miessner, N. Kaplan, J. M. Craig, 
J. Gundlach, A. Laszlo,  S. A. Benner, J. A. Marchand. "Enzymatic Synthesis and Nanopore Sequencing of 12-Letter Supernumerary DNA" 
Nature Communications. 14. (2023). DOI: 10.1038/s41467-023-42406-z 

Updated: 2/14/23
"""
########################################################################
########################################################################

import argparse, textwrap
import os 
import sys 



parser = argparse.ArgumentParser(
        usage='"python xombo.py  [-h] {split, basecall, subset, preprocess, resquiggle, get-level}',
        formatter_class=argparse.RawTextHelpFormatter,
	description= textwrap.dedent('''

********** Xombo *********

Xombo is a suite of tools built from Tombo to streamline processing and level extraction of fast5 files. 

Xombo command groups (additional help available within each command group):
	split		Split multi-fast5 into single fast5 format and merge directory
	basecall 	Basecall single fast5 directory using guppy 
	subset		Subset a portion of single fast5 files to a new directory
	preprocess	Add fastq basecalls to single fast5 files 
	resquiggle 	Resquiggle single fast5 files using Tombo 
	plot-motif	Plot resquiggle aligned to reference mapping using Tombo
	get-level 	Get mean levels of resquiggled data centered around desired 4-kmer
         '''))





subparsers = parser.add_subparsers(dest='subparsers')

parser_split = subparsers.add_parser('split', help='[fast5_dir] [output_dir]')
parser_split.add_argument('fast5_dir',metavar ='[fast5_dir]', type=str, help='Input directory containing multi-fast5 folders.')
parser_split.add_argument('output_dir', metavar ='[output_dir]', type=str, help='Output directory for storing single-fast5 files')

parser_basecall = subparsers.add_parser('basecall', help='[single_fast5_dir] [fastq_output_dir]')
parser_basecall.add_argument('fast5_dir',metavar ='[fast5_dir]', type=str, help='Input directory containing multi-fast5 folders.')
parser_basecall.add_argument('output_dir', metavar ='[output_dir]', type=str, help='Output directory for storing fastq files')
parser_basecall.add_argument('-r',metavar = '[ref_fasta]', type=str, required = False, help='Basecall using a reference fasta file')


parser_subset = subparsers.add_parser('subset', help='[fast5_dir] [output_dir] [n]')
parser_subset.add_argument('fast5_dir', metavar = '[fast5_dir]', type=str, help='Input directory containing single-fast5 files.')
parser_subset.add_argument('output_dir', metavar = '[output_dir]', type=str, help='Directory to move copy of fast5 files. ')
parser_subset.add_argument('n_files', metavar = '[n]', type=int, help='Number of files to move to new directory.')

parser_preprocess = subparsers.add_parser('preprocess', help='[fast5_dir] [fastq_dir]')
parser_preprocess.add_argument('fast5_dir',metavar = '[fast5_dir]', type=str, help='Input directory containing single-fast5 files.')
parser_preprocess.add_argument('fastq_dir',metavar = '[fastq_dir]', type=str, help='Input directory containing fastq basecall files.')

parser_resquiggle = subparsers.add_parser('resquiggle', help='[fast5_dir] [fasta_ref]')
parser_resquiggle.add_argument('fast5_dir',metavar = '[fast5_dir]', type=str, help='Input directory containing basecalled single-fast5 files.')
parser_resquiggle.add_argument('fasta_ref',metavar = '[fasta_ref]', type=str, help='Fasta reference file for mapping.')
parser_resquiggle.add_argument('-g',metavar = '[raw_corr_group]', type=str, required = False, help='Raw corrected group name for storing resquiggle index. Default RawGenomeCorrected_000.')

parser_preprocess = subparsers.add_parser('plot-motif', help='[fast5_dir] [fasta_ref] -m [motif]')
parser_preprocess.add_argument('fast5_dir',metavar = '[fast5_dir]', type=str, help='Input directory containing resquiggled single-fast5 files.')
parser_preprocess.add_argument('fasta_ref',metavar = '[fasta_ref]', type=str, help='Fasta reference file for mapping.')

parser_preprocess.add_argument('-w',metavar = '[window_size]', type=str, help='Window size of base pairs to plot (default = 50)')
parser_preprocess.add_argument('-r',metavar = '[num_regions]', type=str, help='Number of regions to plot (default = 25)')
parser_preprocess.add_argument('-o',metavar = '[output_file]', type=str, help='Output file name e.g. dir/file.pdf. Default plots/tombo_results.motif_centered.pdf')

parser_preprocess = subparsers.add_parser('get-level', help='[fast5_dir] [fasta_ref] [kmer_params]')
parser_preprocess.add_argument('fast5_dir',metavar = '[fast5_dir]', type=str, help='Input directory containing resquiggled single-fast5 files.')
parser_preprocess.add_argument('fasta_ref',metavar = '[fasta_ref]', type=str, help='Fasta reference file for mapping.')
parser_preprocess.add_argument('kmer_params',metavar = '[kmer_params]', type=str, help='File containing parameters for kmer extraction and processing.')
parser_preprocess.add_argument('-o',metavar = '[output_file]', type=str, help='Output file name e.g. dir/levels.csv.')





args = parser.parse_args()


args.subparsers


#Print help if no arguments provided
if args.subparsers==None: 
	parser.print_help()
	sys.exit(0)

#Parse arguments 
if args.subparsers == 'split' :
	os.system('python lib/split.py '+args.fast5_dir+' '+args.output_dir)

elif args.subparsers == 'basecall' and args.r: 
	os.system('python lib/basecall.py '+args.fast5_dir+' '+args.output_dir+' '+args.r)

elif args.subparsers == 'basecall' and not args.r: 
	os.system('python lib/basecall.py '+args.fast5_dir+' '+args.output_dir)

elif args.subparsers == 'subset': 
	os.system('python lib/subset.py '+args.fast5_dir+' '+args.output_dir+' '+str(args.n_files))

elif args.subparsers == 'preprocess': 
	os.system('python lib/preprocess.py '+args.fast5_dir+' '+args.fastq_dir)

elif args.subparsers == 'resquiggle': 
	raw_group = ''
	if args.g:
		raw_group = args.g
	os.system('python lib/resquiggle.py '+args.fast5_dir+' '+args.fasta_ref+' '+raw_group)







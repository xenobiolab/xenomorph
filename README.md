
                                                                  _     
                                                                 | |    
        __  __  ___  _ __    ___   _ __ ___    ___   _ __  _ __  | |__  
        \ \/ / / _ \| '_ \  / _ \ | '_ ` _ \  / _ \ | '__|| '_ \ | '_ \ 
         >  < |  __/| | | || (_) || | | | | || (_) || |   | |_) || | | |
        /_/\_\ \___||_| |_| \___/ |_| |_| |_| \___/ |_|   | .__/ |_| |_|
                                                          | |           
        v1.0                                              |_|           
                                                                        


# Xenomorph : An XNA sequencing toolkit 

## About 
Xenomorph is a suite of tools used for nanopore sequencing of alternative 'xeno'-basepairs (XNAs). This toolkit incorporates ONT-workflows to preprocess FAST5 files and extract signal levels for kmers in a sequence. Models parameterized on XNA basepairs can then be used to test if signal levels match XNA pairs. Xenomorph relies on kmer models that were parameterized using libraries of XNA-containing DNA. The general pipeline consists of two steps: 1) preprocessing FAST5 reads and extracting level information and 2) basecalling using a selected kmer model. This version of Xenomorph was built and tested on Oxford Nanopore Technologies r.9.4.1 flow cells (Flongle or MinION). Xenomorph will continue development and updating models to track the latest releases of Nanopore chemistry. 

This public repository is maintained by the XenoBiology Research Group at the University of Washington, Department of Chemical Engineering. 

## Sample nanopore sequences with XNA basepairs 
This toolkit was created to work with a series of non-standard nucleotides that can form the basis of an expanded genetic alphabet (up to 12 letters). In addition to the standard base pairs (A:T, G:C), the XNA models described in this work allow for single xenonucleotide detection of specific forms of B:S, P:Z, X:K, and J:V. A sample multi-FAST5 file of raw nanopore read data containing PZ insertion can be found in /example_reads/PZ_sample_reads.fast5. Reference file containing all possible sequences in this file, including no-PZ insertion ("Gap") is found in example_reads/PZ_esample_reference.fa. Nanopore XNA sequences used in generation of this work containing BSPZXKJV insertions can be downloaded from the Sequence Reads Archive (SRA Bioproject: PRJNA932328). More information can be found with the associate publication.

## Dependencies
Xenomorph requires ONT tools (ont-fast5-api, tombo, guppy), minimap2 (mappy), and various python packages. A full list of dependencies can be found in the xenomorph-env.yml document. To use conda for installing dependencies, simply load xenomorph-env.yml into a new environment. Xenomorph was built and tested on Ubuntu 18.04 and 20.04 with an Nvidia GTX 3060 GPU. 

        conda env create -f xenomorph-env.yml

To enter xenomorph conda environment, then use: 

        conda activate xenomorph-env


## Xenomorph command groups 
	preprocess	[Preprocess] FAST5 reads with reference FASTA containing XNAs for 'morph' or 'de-novo'
	morph		[Basecall] Use kmer levels identified with preprocess to basecall XNA position based on alternative hypothesis testing (per read)
	extract		[Utility] Extracts raw signal in region associated with XNA bases. 
	stats		[Utility] Calculate global concensus basecalls from per-read output file and generate summary output
	fasta2x		[Utility] Converts FASTA file with XNAs in sequence (e.g. BSPZKXJV) to FASTA file with XNA positional information in header
	models		[Utility] View summary of active and inactive kmer models that can be used for basecalling or activate models


### Preprocess 
Preprocess takes raw nanopore multi-FAST5 files alongside a reference FASTA with XNAs and outputs extracted normalized signal levels for the region surrounding the XNA. Reference FASTA file input should have XNA bases at locations where you want to perform hypothesis testing, which is automatically converted to xfasta format in this pipeline. In this pipeline, raw multi-FAST5 files are first converted to single-FAST5 format then basecalled de novo using guppy. Guppy basecalls are then assigned to the single-fast5 reads. Tombo is then used to match raw signal to bases using the resquiggle command. Using positional header stored in the xFASTA format, normalized signals surrounding XNA positions are then extracted and stored for each read. Any non-ATGC base can be used at this stage for level extraction. However, alternative hypothesis testing basecalling will require using base abbreviations found in kmer model files (accessible with xenomorph models). Intermediate files are output in working directory. Following preprocess command, the output summary file can be used as an input for alternative hypothesis testing with the xenomorph.py morph command. 

        python xenomorph.py preprocess -w [working directory] -f [fast5 diretory] -r [reference_fasta] 
                Optional flags: 
                -b force_rebasecall_with_guppy
                -x force_reprocessing
                -o [output_summary_file_name.csv]

### Morph
Morph uses alternative hypothesis testing to determine if the expected XNA is truly the best basecall for that position, given the levels extracted from 'xenomorph preprocess'. Morph input file is the output summary file (.csv) generated from the 'xenomorph.py preprocess' command. A kmer model is required as an input. Pre-calculated kmer models can be accessed via 'xenomorph.py models' and specified using the -m flag (e.g. -m ATGCBS). Alternatively, custom models can be input by providing a .csv model as an input for the -m flag. Additional model parameters can be changed in the lib/xm_params.py file. The output of the morph command is a .csv file contains most likely base given the possible bases at that given position on a per read basis. Most likely base is determined by maximum likelihood estimation. Log-likelihood ratios in the output file compare the XNA basecall to the next most likely (XNA prefered if >0). Outlier-robust log-likelihood is the default statistic used for determining basecalls, and can be modified in lib/xm_params. 

        python xenomorph.py morph -l [level_summary_file.csv] -m [kmer_model.csv OR ATGCXY]
                Optional flags: 
                -o [output_basecalled_summary_file_name.csv]

### Extract
Extract performs a resquiggle on each read to extract the raw signal (unsegmented) around the XNA. This feature is intended for method development, as it allows users to re-analyze signal associated with XNA regions. The boundaries for the raw signal are determined by the xmer_boundary parameter in the lib/xm_params.py file, and should match the settings used during preprocessing. The output of this command is a .csv file which contains individual reads on each row with their associated raw signal (unsegmented, not normalized), normalized signal (unsegmented), and segmented normalized signal. Segmented normalized signals are the same as found in the level output file generated by xenomorph preprocess. Since region of raw signal extraction is specified by the preprocess level file, users will need to run xenomorph.py preprocess first to generate a level file (e.g., level_summary_file.csv) and provide this as an input to xenomorph.py extract. 

        python xenomorph.py morph -w [working directory] -f [fast5 diretory] -r [reference_fasta] -l [level_summary_file.csv] -o [output_file_name]

### Stats
Stats command takes the per read output basecall file (.csv) generated from morph and calculates a global summary. Of note, stats identifies consensus base called for each reference read. If more than one base is called at the same frequency, both are output. 

        python xenomorph.py stats -i [output_basecalled_summary_file.csv]

### Models
Linked kmer models can be accessed using the 'xenomorph.py models' command. As a general paradigm for modular design, xenomorph models are considered to be fully orthogonal and completely modular. Xenomorph comes built with empirically determined models for standard bases (ATGC) and a expanded junimoji alphabet (BSPZKXJV). Models used for the morph command can be specified by using the '-m' flag and inputing the desired alphabet (e.g. ATGCBSJV). Note, ATGC should always be specified as part of this alphabet since xenomorph is designed for ground-up alphabet expansion (ATGC+) as opposed to top-down (BSPZ-only). Unlike the standard bases, there are various variants of xenobases that are not mutually orthogonal (e.g. N-glycoside S and C-glycoside S). To select different chemical variants of expanded bases, the '-a' [abbreviation] flag can be used. This command will set the abbreviated base as the active XNA variant and disable alterntative variants. For example, 'xenomorph.py models -s Sn' will set the N-glycoside S model as the active basecalling model. To view which models are available, single letter codes, abbreviations, and active/inactive state use the '-s all' flag. Models can be linked, edited, added, or deleted by modifying the models/config_model.csv parameter file. 


        python xenomorph.py models 
                Optional flags: 
                -s [active/inactive/all]
                -a [base_abbreviation_to_activate]

### xFASTA format 
Many tools used to read and manipulate nucleic acid sequences are not built to handle non ATGC bases. In an effort to streamline how we handle XNAs in a 'FASTA-like format', xenomorph was built to handle non-standard bases in FASTA sequences (e.g. BSPZJVKX). The xFASTA file format (.fa) stores XNA positional information in the header of each sequence. xFASTA files can be generated from standard Fasta files that contain non-ATGC bases (e.g. BSPZJVKX) in the sequence. xFASTA files are automatically generated in the standard xenomorph preprocessing workflow. The fasta2x command is provided for utility, but generally not required. Note that XNA bases in the input sequence will be replaced for a standard ATGC. Signal matching to sequence is highly sensitive what base is chosen as the replacement base. As default, XNA bases are replaced as followed: B>G, S>C, P>G, Z>G. Base substitution settings can be modified in lib/xm_params.py by changing the paired base in the confounding_base variable. 


        Standard FASTA format with XNA in sequence
                >header_for_a_read
                ATGGCAACAGGATGABAAGGACGTA

        xFASTA format with XNA information stored in header. (B replaced with G in sequence)
                >header_for_a_read+X_POS[B:18]
                ATGGCAACAGGATGAGAAGGACGTA

### Modifying parameters
Preprocess, morph, and stats use various parameters that can be tuned or modified for desired application. Though default parameters are recommended for most applications, parameters can be modified by editing lib/xm_params.py. 



## Examples

#### Preprocessing with basecalling and resquiggle (full) 


        1) Preprocess a raw nanopore run (convert, basecall, resquiggle, and level extract) 
            python xenomorph.py preprocess -w path/to/output/directory -f path/to/fast5/folder -r path/to/fasta.fa -o preprocess_output_file.csv -b -x 

#### Preprocessing without basecalling and resquiggle (partial) 


        1) Preprocess a raw nanopore run but skip basecalling and resquiggle if previous performed (level extract only). Requires having run full preprocessing pipeline previously
            python xenomorph.py preprocess -w path/to/output/directory -f path/to/fast5/folder -r path/to/fasta.fa -o preprocess_output_file.csv

#### Basecall a preprocessed file 


        1) Using prebuilt XNA model for PZ 
            python xenomorph.py morph -l path/to/preprocess_output_file.csv -m ATGCPZ -o path/to/basecall_output_file.csv
        2) Using prebuilt XNA model for PZBS 
            python xenomorph.py morph -l path/to/preprocess_output_file.csv -m ATGCBSPZ -o path/to/basecall_output_file.csv
        3) Using prebuilt XNA model for PZBSJVKX 
            python xenomorph.py morph -l path/to/preprocess_output_file.csv -m ATGCBSPZJVKX -o path/to/basecall_output_file.csv
        4) Using custom kmer-model file (should be .csv with format similar to models/*.csv files)
            python xenomorph.py morph -l path/to/preprocess_output_file.csv -m path/to/kmer/model.csv -o path/to/basecall_output_file.csv


#### Preprocess, basecall, generate global summary (end-to-end pipeline)


        1) Preprocess a raw nanopore run (convert, basecall, resquiggle, and level extract) for ATGCBS containing reads
            python xenomorph.py preprocess -w path/to/output/directory -f path/to/fast5/folder -r path/to/fasta.fa -o preprocess_output_file.csv -b -x 
        2) Perform alternative hypothesis testing on a per-read basis to obtain basecalls 
            python xenomorph.py morph -l path/to/preprocess_output_file.csv -m ATGCBS -o path/to/basecall_output_file.csv
        3) Generate global summary file from the per-read output basecall file 
            python xenomorph.py stats -i path/to/basecall_output_file.csv


#### Preprocess and extract raw signal around XNA regions


        1) Preprocess a raw nanopore run (convert, basecall, resquiggle, and level extract) for ATGCBS containing reads
            python xenomorph.py preprocess -w path/to/output/directory -f path/to/fast5/folder -r path/to/fasta.fa -o preprocess_output_file.csv -b -x 
        2) Perform raw signal extraction around region specified by level file 
            python xenomorph.py extract -w path/to/output/directory -f path/to/fast5/folder -r path/to/fasta.fa -l path/to/preprocess_output_file.csv -o path/to/raw_signal_output_file.csv



#### Model commands


        1) View all possible models accessible for the -m morph flag
            python xenomorph.py models -s all 
        2) View active models accessible for the -m morph flag
            python xenomorph.py models -s active 
        3) Activate the C-glycoside S model and inactivate the N-glycoside S model
            python xenomorph.py models -a Sc 
        4) Activate the N-glycoside S model and inactivate the C-glycoside S model
            python xenomorph.py models -a Sn 

#### Generate XNA model from raw reads and a reference FASTA


        1) Preprocess a raw nanopore run (convert, basecall, resquiggle, and level extract) 
            python xenomorph.py preprocess -w path/to/output/directory -f path/to/fast5/folder -r path/to/fasta.fa -o preprocess_output_file.csv -b -x 
        2) Extract all measurements of each kmer to generate a kmer output file 
            python lib/xm_extract_levels.py preprocess_output_file.csv 
        3) Convert kmer level file into a model summary file 
            python lib/parse_kmer.py preprocess_output_file_kmers model_output.csv
        4) To use model for basecalling, add entry of model file by linking in model/config_model.csv. Model abbreviation and letter name need to be specified. 
        5) Set model to active (True) in model/config_model.csv OR use: python xenomorph.py models -a [Model_Abbreviation] 



## Cite us or read more about this work 
    Title: Synthesis and Sequencing of a 12-Letter Supernumerary DNA

    By: H. Kawabe, C. Thomas, S. Hoshika, Myong-Jung Kim, Myong-Sang Kim, L. Miessner, J. M. Craig, 
    J. Gundlach, A. Laszlo,  S. A. Benner, J. A. Marchand



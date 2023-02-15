########################################################################
########################################################################
"""
xm_fasta2x.py

Description: This script is used for handling fasta conversion where 
XNA sequences are part of the fasta sequence. Positions and identity of 
XNA are placed in fasta header generating xfasta files.Input fasta files should 
contain XNAs with sequences that contain XNA bases in xm_params.py. 

Title: Synthesis and Sequencing of 12-Letter Supernumerary DNA

By: H. Kawabe, C. Thomas, S. Hoshika, Myong-Jung Kim, Myong-Sang Kim, L. Miessner, J. M. Craig, 
J. Gundlach, A. Laszlo,  S. A. Benner, J. A. Marchand

Updated: 2/14/23
"""
########################################################################
########################################################################


from xm_params import *
import sys
import numpy as np


input_fasta=sys.argv[1]
output_fasta=sys.argv[2]




#Create writable output file
f = open(output_fasta, "w")
detected_xfasta_header = False 
detected_xna = False


with open(input_fasta, "r") as fh:
    for line in fh:
    
        #Get header
        if line[0]=='>':
            header = line
            if 'XPOS[' in header: 
                detected_xfasta_header=True
            
        #Get sequence
        if line[0]!='>':
        
            #Make all upper case
            uline = line.upper()
            
            #Look for non-standard base
            diff = list(set(uline.replace('\n',''))-(set(standard_bases)))
            
            #This sequence contains an XNA
            if len(diff)>0:
            
                #Format header

                
                #Purge XNAs from sequence

                uline_gap = uline #Null position
                
                #Get position of each XNA
                for x in range(0,len(diff)):

                    #Location of every modification 
                    x_loc=[i for i in range(len(uline)) if uline.startswith(diff[x], i)]
                    

                    #For every substitution base
                    for b in substitution_bases: 
                    #Setup output header 
                        header_x = header.replace('\n','')+'+'+b+'+XPOS['
                        uline_clean = uline #Standard base set
                        #At every modification position 
                        for xi in x_loc:
                            posx=diff[x]+':'+str(xi)+'-'
                            header_x = header_x+posx
                            uline_clean = uline_clean.replace(diff[x],b)
                            uline_gap = uline_gap.replace(diff[x],'-')
                            detected_xna = True 

                        #Close and write 
                        header_xw=header_x[:-1]+']\n'
                        f.write(header_xw)
                        f.write(uline_clean)

                    if write_gaps==True: 
                        header_gap = header.replace('\n','')+'+-+_GAP[]\n'
                        uline_gap = uline #Standard base set
                        for xi in x_loc:
                            uline_gap = uline_gap.replace(diff[x],'-')
                        f.write(header_gap)
                        f.write(uline_gap.replace('-',''))
            elif len(diff)==0 and write_no_xna_seq==True: 
                f.write(header)
                f.write(uline)







if detected_xfasta_header == True: 
    print('Xenomorph Status - [Error] Fasta input file already in xfasta format')
    fasta_input_error=True 
else: 
    if detected_xna == False: 
        print('Xenomorph Status - [Error] No XNAs (BS/PZ/KX/JV/XY) detected in fasta input sequence.')
        fasta_input_error=True 


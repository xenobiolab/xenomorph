########################################################################
########################################################################
"""
xm_xfasta2bed.py 

Description: Converts a xfasta formatted file into bed file. This is a 
utility tool created for handling fasta with XNA bases. 

Usage: python xenomorph.py -h for help 

Title: XNA tailing enables nanopore sequencing of a 12-letter genetic code

By: H. Kawabe, C. Thomas, A. Laszlo, S. Hoshika, L. Miessner, J. M. Craig, 
J. Gundlach, Myong-Jung Kim, Myong-Sang Kim, S. A. Benner, J. A. Marchand

Updated: 1/15/23
"""
########################################################################
########################################################################



from xm_params import *
from xm_tools import *
import sys
import numpy as np
from Bio.Seq import Seq
import os


###Input fasta and XNA bases 
input_fasta = sys.argv[1]
xna_base = sys.argv[2]


def xfasta2bed(input_fasta, xna_base):
    output_fasta = os.path.expanduser(input_fasta).replace('.fa','.bed')
    fr = open(output_fasta,"w")
    with open(os.path.expanduser(input_fasta), "r") as fo:
        for line in fo: 

            if 'GAP' not in line: 
                if line[0]=='>':
                    header = line[1:].replace('\n','')
                    x_pos_base = fetch_xna_pos(header)
                    x_pos_to_rc =[]


                    for x in x_pos_base: 
                        x_base = x[0]

                        x_pos = int(''.join(filter(str.isdigit, x[1])))

                        if x_base == xna_base: 
                            strand = '+'
                            bed_line = header+'  '+str(x_pos)+'  '+str(int(x_pos)+1)+'  '+x_base+'  0  '+strand+'\n'

                            fr.write(bed_line)
                        elif x_base == xna_base_rc(xna_base,xna_base_pairs): 
                            x_base = xna_base
                            strand ='-'
                            bed_line = header+'  '+str(x_pos)+'  '+str(int(x_pos)+1)+'  '+x_base+'  0  '+strand+'\n'

                            fr.write(bed_line)
    fr.close()

xfasta2bed(input_fasta,xna_base)
                    
                    


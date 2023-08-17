########################################################################
########################################################################
"""
xr_xfasta2bed.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. A. Marchand

Updated: 3/2/23
"""
########################################################################
########################################################################


from xr_params import *
from xr_tools import *
import sys
import numpy as np
from Bio.Seq import Seq
import os


input_fasta = sys.argv[1]
output_file = sys.argv[2]
xna_base = sys.argv[3]
sub_base = sys.argv[4]



def xfasta2bed(input_fasta, xna_base, output_file, sub_base):

    fr = open(output_file,"w")
    with open(os.path.expanduser(input_fasta), "r") as fo:
        for line in fo: 

            if 'GAP' not in line.upper(): 
                if line[0]=='>':
                    header = line[1:].replace('\n','')
                    x_pos_base = fetch_xna_pos(header)
                    x_pos_to_rc =[]


                    for x in x_pos_base: 
                        x_base = x[0]

                        x_pos = int(''.join(filter(str.isdigit, x[1])))

                        if x_base == xna_base: 
                            strand = '+'
                            bed_line = header+'  '+str(x_pos)+'  '+str(int(x_pos)+1)+'  '+sub_base+'  0  '+strand+'\n'
                            fr.write(bed_line)


                        elif x_base == xna_base_rc(xna_base,xna_base_pairs): 
                            x_base = xna_base
                            strand ='-'
                            bed_line = header+'  '+str(x_pos)+'  '+str(int(x_pos)+1)+'  '+sub_base+'  0  '+strand+'\n'

                            fr.write(bed_line)
    fr.close()


xfasta2bed(input_fasta,xna_base, output_file, sub_base)
                    
                    


########################################################################
########################################################################
"""
xr_rescale.py 
This is a rescaling function that uses a selected scaling method to correct for remora-based scaling. 
Rescaling is generally required if Remora segmentation is chosen, and not required for Tombo segmentation. 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. A. Marchand

Updated: 8/16/23


"""

#######################################################################



#import everything
import re
import os
import sys
import subprocess

from alive_progress import alive_bar; import time
import pod5
import pysam
import numpy as np
import pandas as pd
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

from remora import io, refine_signal_map, util
from time import sleep
from tqdm import tqdm
from csv import writer
from Bio import SeqIO
from Bio.Seq import Seq
from xm_tools import *
from xr_tools import *
from xm_params import *
from xr_params import *




def global_rescale(level_file_input, kmer_model)


    return level_file_output


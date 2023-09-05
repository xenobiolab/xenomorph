########################################################################
########################################################################
"""
xm_tools.py 

Description: Contains various functions used for handling XNA sequences. 

Title: Synthesis and Sequencing of a 12-letter supernumerary DNA

By: H. Kawabe, C. Thomas, S. Hoshika, Myong-Jung Kim, Myong-Sang Kim, L. Miessner, J. M. Craig, 
J. Gundlach, A. Laszlo,  S. A. Benner, J. A. Marchand

Updated: 9/4/23
"""
########################################################################
########################################################################


import pandas as pd


def fetch_xna_pos(xm_header):
    pos=xm_header[xm_header.find('XPOS[')+5:-1].split('-')
    xpos = [x.split(':') for x in pos]
    return xpos


def xna_base_rc(xna_base, xna_bp): 
	for x in xna_bp: 
		if xna_base in x:
			if len(x)>1:
			    xx=list(x)
			    xx.remove(xna_base)
			    xrc=xx[0]
			else:
			   xrc=False
	return xrc


def check_xfasta_format(xfasta_file,standard_bases): 
    xfasta_header=False
    xna_in_sequence = False 
    with open(xfasta_file, "r") as fh:
        for line in fh:
            #Get header
            if line[0]=='>':
                header = line
                if 'XPOS[' in header: 
                    xfasta_header=True
                
            #Get sequence
            if line[0]!='>':
            
                #Make all upper case
                uline = line.upper()
                
                #Look for non-standard base
                diff = list(set(uline.replace('\n',''))-(set(standard_bases)))
                
                #This sequence contains an XNA
                if len(diff)>0:
                    xna_in_sequence = True 
    if xfasta_header==True and xna_in_sequence == False: 
        return True 
    else: 

        return False 


def compile_model(model_files, model_bases): 
    kmer_model =[] 
    for i in range(0,len(model_files)):
	    ifn = model_files[i]
	    model_part = pd.read_csv(ifn, sep=',', dtype={'kmer_xy': str, 'mean_level': object})
	    model_part = model_part[model_part['KXmer'].str.contains(model_bases[i])]
	    kmer_model.append(model_part)
    kmer_model = pd.concat(kmer_model).reset_index()


    kmer_model = kmer_model.drop_duplicates(subset='KXmer', keep="first")
    km = list(set(kmer_model['KXmer']))
    for i in range(0,len(km)): 
        kk = kmer_model[kmer_model['KXmer']==km[i]]
        if len(kk)>1:
            print(kk)
    return kmer_model 


def parse_model_files(model_code, active_status): 
    model_config = pd.read_csv('models/config_model.csv') 
    active_model = model_config[model_config['active']==active_status]
    model_files=[]
    for i in range(0,len(model_code)): 
        try: 
            model_files.append(active_model[active_model['letter_code'].str.contains(model_code[i])]['path '+flowcell_version].values[0])
        except: 
        ....print('Xenomorph Status - [Warning] Unsupported flowcell version detected. Model defaulting to 10.4.1')
            model_files.append(active_model[active_model['letter_code'].str.contains(model_code[i])]['path 10.4.1'].values[0])
    return model_files


def model_summary(active_status): 
    model_config = pd.read_csv('models/config_model.csv') 
    print('\n')
    if active_status.lower() == 'active': 
        model_config = model_config[model_config['active']==True]
        print(model_config[['letter_code','abbreviation','description', 'active']])
    elif active_status.lower() == 'inactive': 
        model_config = model_config[model_config['active']==False]
        print(model_config[['letter_code','abbreviation','description', 'active']])
    else: 
        print(model_config[['letter_code','abbreviation','description', 'active']])
    print('\n')
    print('***Kmer models can be configured by editing models/config_models.csv')
    print('***Use xenomorph models -a [ABBREVIATION] to set models as active')



def activate_model(base_abbreviation): 
    model_config = pd.read_csv('models/config_model.csv') 
    l_code = model_config[model_config['abbreviation']==base_abbreviation]['letter_code'].values[0]
    model_config.loc[model_config['abbreviation']==base_abbreviation,'active']=True 
    model_config.loc[(model_config['abbreviation']!=base_abbreviation) & (model_config['letter_code']==l_code), 'active']=False
    model_config.to_csv('models/config_model.csv', index = False) 



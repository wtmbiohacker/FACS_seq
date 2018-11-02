# this script is used as the main function for the HTsensor pipeline

import os
import sys
import pandas

config_file=sys.argv[1]
config={}
with open (config_file,'r') as f:
    for line in f:
        line=line.strip()
        if len(line)==0:
            continue
        if line[0]=='#':
            continue
        row=line.split('\t')
        config[row[0]]=[]
        if len(row)==1:
            config[row[0]]=''
        elif len(row)==2:
            config[row[0]]=row[1]
        else:
            config[row[0]]=row[1:]

#read the fastq file and corresponding label
for i in range(1,len(config['fastq'])):
    config['fastq'][i]=config['fastq'][0]+config['fastq'][i]
del config['fastq'][0]
fastq=','.join(config['fastq'])

#read the corresponding label of fastq file
if 'sample-label' not in config or config['sample-label']=='':
    sample_label=''
else:
    sample_label=','.join(config['sample-label'])

#the label of initial library
if 'iniLib' not in config or config['iniLib']=='':
    print 'incorrect setting about initial library'
    sys.exit()
else:
    iniLib=config['iniLib']

#Prefix nucleotide to identify the variable region
if 'prefix_nucl' not in config or config['prefix_nucl']=='' or 'True' in [(x not in 'ATCG') for x in config['prefix_nucl'].upper()]:
    print 'incorrect prefix nucleotide input!'
    sys.exit()
else:
    prefix_nucl=config['prefix_nucl']

#Suffix nucleotide to identify the variable region
if 'suffix_nucl' not in config or config['suffix_nucl']=='' or 'True' in [(x not in 'ATCG') for x in config['suffix_nucl'].upper()]:
    print 'incorrect suffix nucleotide input!'
    sys.exit()
else:
    suffix_nucl=config['suffix_nucl']

#specify the length of the sgRNAs (without PAM sequence)
if 'variable_region_len' not in config or config['variable_region_len']=='':
    print('Please tell the program the length of the variable region!')
    sys.exit()
else:
    variable_region_len=config['variable_region_len']

#Path to library design file (csv format, columns: id, sequence, gene)
list_seq=config['list-seq']

#Sub configure file to define the experimental setting
experiment_configure=config['experiment_configure']

#Sub configure file to define the cell number sorted into each bin in the original experiment
cell_count_configure=config['cell_count_configure']

#Sub configure file to define the Pij (occupation of each bin) for each sorting experiemnt
Bin_occupation_configure=config['bin_occupation_configure']

#Sub configure file to define the boundaries between bins for each sorting experiemnt
Bin_boundary_configure=config['bin_boundary_configure']

#number of reads in initial library that below this threshold won't be incoperated in to the following analysis pipleine 
if 'ReadsThreshold' not in config or config['ReadsThreshold']=='':
    ReadsThreshold=100
    print('the default value of ReadsThreshold is 100!')
else:
    ReadsThreshold=int(config['ReadsThreshold'])

# range of Log10u and its deviation during optimization, format:lower limit,upper limit
if 'log10u_sigma_range' not in config or config['log10u_sigma_range']=='':
    print('Please tell the program the search range of log10u and sigma!')
    sys.exit()
else:
    log10u_sigma_range=config['log10u_sigma_range']

# number of points in the grid space of each dimension during optimization
if 'search' not in config or config['search']=='':
    search=20
    print('the default value of ReadsThreshold is 20!')
else:
    search=int(config['search'])

# prefix for naming of all files
if 'prefix' not in config or config['prefix']=='':
    prefix='Test'
    print('the default value of prefix is Test!')
else:
    prefix=config['prefix']

# ////////////////////////////////////////////////////////////////////
# Process the raw .fq NGS data
outputname='%s_results/%s_rawcount/%s'%(prefix,prefix,prefix)
#os.system('python HTsensor_fqpre.py --unmapped-to-file --fastq %s --sample_label %s --prefix_nucl %s --suffix_nucl %s --variable_region_len %s --list_seq %s --output_prefix %s 2>error.log'%(fastq,sample_label,prefix_nucl,suffix_nucl,variable_region_len,list_seq,outputname))
print('NGS data pretreatment finished')

# ////////////////////////////////////////////////////////////////////
# Clean the preprocessed data and construct necessary Python objects for sensor response profile inferring
rawData='%s_results/%s_rawcount/%s.count.txt'%(prefix,prefix,prefix)
unmap='%s_results/%s_rawcount/%s.unmapped.txt'%(prefix,prefix,prefix)
outputname='%s_results/%s_cleandataset/%s'%(prefix,prefix,prefix)
#os.system('python HTsensor_predataset.py --ctab %s --unmap %s --readth %s --exp_con %s --ini_lib %s --output_prefix %s 2>>error.log'%(rawData,unmap,ReadsThreshold,experiment_configure,iniLib,outputname))
print('Clean dataset preparation finished')

# ////////////////////////////////////////////////////////////////////
# optimization process to infer the Log10u and sigma of each sensor mutant at each experiment condition
ini_ab='iniab_Dic.pickle'
ctab='sensor_read_Dic.pickle'
total='total_read_DF.pickle'
outputname='%s_results/%s_optimization/%s'%(prefix,prefix,prefix)
#os.system('python HTsensor_optimize.py --ini_ab %s --ctab %s --total %s --exp_con %s --cell_con %s --bin_occ_con %s --bin_bou_con %s --rrange %s --search %s --output_prefix %s 2>>error.log'%(ini_ab,ctab,total,experiment_configure,cell_count_configure,Bin_occupation_configure,Bin_boundary_configure,str(log10u_sigma_range),search,outputname))
print ('Bayesian optimization process finished')

# ////////////////////////////////////////////////////////////////////
# plot
Log10u_csv='sensor_Log10u.csv'
sigma_csv='sensor_sigma.csv'
replicate_con='replicae_configure.txt'
outputname='%s_results/%s_plot/%s'%(prefix,prefix,prefix)
os.system('python HTsensor_plot.py --replicate_con %s --rrange %s --log10u %s --sigma %s --output_prefix %s 2>error.log'%(replicate_con,str(log10u_sigma_range),Log10u_csv,sigma_csv,outputname))
print ('plotting finished')

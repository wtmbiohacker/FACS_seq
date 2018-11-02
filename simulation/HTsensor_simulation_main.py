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

#Sub configure file to define the experimental setting
experiment_configure=config['experiment_configure']

#Sub configure file to define the cell number sorted into each bin in the original experiment
cell_count_configure=config['cell_count_configure']

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

# range of Log10u and its deviation during optimization, format:lower limit,upper limit
if 'log10u_sigma_dis' not in config or config['log10u_sigma_dis']=='':
    print('Please tell the program the distribution of log10u and sigma!')
    sys.exit()
else:
    log10u_sigma_dis=config['log10u_sigma_dis']

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
# generate random library
outputname='%s_results/%s_library/%s'%(prefix,prefix,prefix)
print ('python HTsensor_simulation_generateLib.py --exp_con %s --bin_bou_con %s --distribution %s --readth %s --output_prefix %s 2>error.log'%(experiment_configure,    Bin_boundary_configure,log10u_sigma_dis,str(ReadsThreshold),outputname))
os.system('python HTsensor_simulation_generateLib.py --exp_con %s --bin_bou_con %s --distribution %s --readth %s --output_prefix %s 2>error.log'%(experiment_configure,Bin_boundary_configure,str(log10u_sigma_dis),str(ReadsThreshold),outputname))
print('Simulated dataset prepared')

# ////////////////////////////////////////////////////////////////////
# perform the optimization as for real data 
ini_ab='iniab_Dic.pickle'
ctab='sensor_read_Dic.pickle'
total='total_read_DF.pickle'
Bin_occupation_configure='bin_occupation_configure.txt'
cell_count_configure='cell_count_configure.txt'
outputname='%s_results/%s_optimization/%s'%(prefix,prefix,prefix)
print ('python HTsensor_optimize.py --ini_ab %s --ctab %s --total %s --exp_con %s --cell_con %s --bin_occ_con %s --bin_bou_con %s --rrange %s --search %s --output_prefix %s 2>>error.log'%(ini_ab,ctab,total,experiment_configure,cell_count_configure,Bin_occupation_configure,Bin_boundary_configure,str(log10u_sigma_range),search,outputname))
os.system('python HTsensor_optimize.py --ini_ab %s --ctab %s --total %s --exp_con %s --cell_con %s --bin_occ_con %s --bin_bou_con %s --rrange %s --search %s --output_prefix %s 2>>error.log'%(ini_ab,ctab,total,experiment_configure,cell_count_configure,Bin_occupation_configure,Bin_boundary_configure,str(log10u_sigma_range),search,outputname))
print ('Bayesian optimization process finished')

# ////////////////////////////////////////////////////////////////////
# Compare the calculation result with the inherent design of the mutant
my_lib='my_library.txt'
log10uResult='%s_results/%s_optimization/sensor_Log10u.csv'%(prefix,prefix)
sigmaResult='%s_results/%s_optimization/sensor_sigma.csv'%(prefix,prefix)
outputname='%s_results/%s_comparison/%s'%(prefix,prefix,prefix)
print ('python htsensor_simulation_compare.py --my_lib %s --log10u %s --sigma %s --rrange %s --output_prefix %s 2>>error.log'%(my_lib,log10uResult,sigmaResult,str(log10u_sigma_range),outputname))
os.system('python htsensor_simulation_compare.py --my_lib %s --log10u %s --sigma %s --rrange %s --output_prefix %s 2>>error.log'%(my_lib,log10uResult,sigmaResult,str(log10u_sigma_range),outputname))
print('Comparison finished!')


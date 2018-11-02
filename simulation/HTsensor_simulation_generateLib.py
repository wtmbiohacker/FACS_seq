# this script read in the readcount.csv file; and loop to find the best parameter set for each sequence

from __future__ import print_function
import os
import sys
import numpy as np
import pickle
import scipy.stats
import math
import pandas as pd
import logging
import argparse
import matplotlib.pyplot as plt
from random import shuffle

def htsensorsimu_generate_parseargs():
    """
    Parse arguments. Only used when htsensorsimu_generate.py is executed directly.
    """
    parser=argparse.ArgumentParser(description='Collecting read counts for multiple samples.')
    parser.add_argument('-n','--lib_size', type=int,default=1500,help='number of sensor mutants')
    parser.add_argument('-e','--exp_con',required=True,help='A file containing the design of the experiment, treatment of ligand concentration as well as the bins of cell sorting. Support file format: csv with tab delimiter.')
    parser.add_argument('-b','--bin_bou_con',required=True,help='A file containing the boundary of each bin in cell sorting experiment. Support file format: csv with tab delimiter.')
    parser.add_argument('-d','--distribution',required=True,help='string to define the distribution of Log10u and sigma. Format:(lower,upper(Log10u);lower,upper(sigma))')
    parser.add_argument('-t','--total_reads',type=int,default=10000000,help='total read number per sequencing library')
    parser.add_argument('-c','--total_cells',type=int,default=100000,help='total cell number in sorting experiment')
    parser.add_argument('-r', '--readth', type=int, default=20,help='read threshold in initial library to exclude sensor mutant for further analysis')
    parser.add_argument('-o','--output_prefix',default='sample1',help='The prefix of the output file(s). Default sample1.')
    args=parser.parse_args()
    return args

# ///////////////////////////////////////////////////////
def construct_library(lib_size,range_Log10u,range_sigma):
    """
    Construct a random library for simulation
    Parameters
    __________
    lib_size:positive int
    range_Log10u:[mean,sigma] of Log10u
    range_sigma:[mean.sigma] of sigma
    Return
    __________
    pandas DataFrame of my library
    columns: 'Log10u' 'sigma' 'Abundance'
    rows: mutants from mutant0 to mutantX
    """
    my_lib={}
    #Log10u_mean=range_Log10u[0]
    #Log10u_dev=range_Log10u[1]
    #sigma_mean=range_sigma[0]
    #sigma_dev=range_sigma[1]
    Log10u0=range_Log10u[0]
    Log10u1=range_Log10u[1]
    sigma0=range_sigma[0]
    sigma1=range_sigma[1]
    temp_array=(np.random.rand(lib_size)+0.2)
    # with a typical Gini index ~ 0.25
    ab_array=temp_array/np.sum(temp_array)
    for i in range(lib_size):
        mutant='mutant%d'%(i)
        #Log10u=np.random.normal(Log10u_mean,Log10u_dev,1)[0]
        #sigma=np.random.normal(sigma_mean,sigma_dev,1)[0]
        Log10u=np.random.uniform(Log10u0,Log10u1,1)[0]
        sigma=np.random.uniform(sigma0,sigma1,1)[0]
        my_lib[mutant]={}
        my_lib[mutant]['Log10u']=Log10u
        my_lib[mutant]['sigma']=sigma
        my_lib[mutant]['Abundance']=ab_array[i]
    my_lib_DF=pd.DataFrame(my_lib).T
    my_lib_DF.index.name='Sensor'
    return my_lib_DF


def Cal_area(Log10u,sigma,BinBoundary_Lst):
    """
    Calculate the area under the normal distribution curve in defined space: N(Log10u,sigma)[F0,F1]
    Using this funciton, we get pj,sensor, the possibility of a sensor seq into a particular bin under relevant ND assumption
    Parameters
    __________
    Log10u, mean of the normal distribution
    sigma, deviation of the normal distribution
    BinBoundary_Lst=[A0,A1], the defined X axis space of calculation
    #
    Return
    __________
    0=<area<=1
    """
    if len(BinBoundary_Lst)==2 and BinBoundary_Lst[1]>BinBoundary_Lst[0]:
        F0=BinBoundary_Lst[0]
        F1=BinBoundary_Lst[1]
        area=scipy.stats.norm(Log10u,sigma).cdf(F1)-scipy.stats.norm(Log10u,sigma).cdf(F0)
        return area
    else:
        logging.error('incorrect bin boundary!')
        sys.exit(-1)


def cytometry_sorting(pop_DF,exp_DF,Bin_bou_con,total_cells):
    """
    Funciton used to simulate a FACS experiment 
    Parameters
    __________
    pop_DF: pandas dataframe of the infomation about each sensor mutant
    columns = Log10u, sigma, relative abundance
    rows = sensor mutant
    values = real number
    #
    exp_DF: pandas DataFrame of experiment setting
    columns = ligand1R1 ligand1R2 ligand2R1 ligand2R2 ... (experiment condition)
    rows = bins of sorting
    values = library label
    #
    Bin_bou_con: sub configure file to define the boundaries between bins for each sorting experiemnt
    columns = ligand1R1 ligand1R2 ligand2R1 ligand2R2 ... (experiment condition)
    rows = bins of sorting
    values = [LogA0,LogA1]
    #
    total_cells: number of total cells sorted in this experiment
    #
    Return
    __________
    tuple:
    1. Dict of experiments targeting to relevant pandas DF of sorting result:
    {exp:DF}
    DF:
    columns = bins
    rows = sensor mutant
    values = integer
    2. DataFrame of bin occupations
    columns: experiment conditions
    rows: bins
    values: relative abundance of each bin against its relevant experiment, sum(every column)=1
    3. DataFrame of bin cell number
    columns: experiment conditions
    rows: bins
    values: absolute cell number sorted into each bin
    """
    sort_Dic={}
    exp_Lst=exp_DF.columns.tolist()
    bins_Lst=exp_DF.index.tolist()
    if Bin_bou_con.index.tolist()!=bins_Lst:
        logging.error('Contradictory bin setting!')
        sys.exit(-1)
    sensor_Lst=pop_DF.index.tolist()
    sort_Dic={}
    for exp in exp_Lst:
        sort_Dic[exp]={}
        for sensor in sensor_Lst:
            sort_Dic[exp][sensor]={}
            Log10u=pop_DF.loc[sensor,'Log10u']
            sigma=pop_DF.loc[sensor,'sigma']
            cells=pop_DF.loc[sensor,'Abundance']*total_cells
            sort_Dic[exp][sensor]=(pd.Series({bins:Cal_area(Log10u,sigma,Bin_bou_con.loc[bins,exp]) for bins in bins_Lst})*cells).to_dict()
        sort_Dic[exp]=pd.DataFrame(sort_Dic[exp]).T
    cell_con_Dic={}
    for exp in exp_Lst:
        cell_con_Dic[exp]={bins:np.sum(sort_Dic[exp][bins]) for bins in bins_Lst}
    cell_con_DF=pd.DataFrame(cell_con_Dic)
    cell_con_DF.index.name = 'Bin'
    bin_occ_DF=cell_con_DF.div(cell_con_DF.sum(),axis=1)
    return (sort_Dic,bin_occ_DF,cell_con_DF)

def sequencing(pop_Series, total_reads):
    """
    Funciton used to simulate a sequencing experiment 
    Parameters
    __________
    pop_DF: pandas Series of the abundance (number) about each sensor mutant
    sensor_mutant ab (number)
    #
    total_reads: total number of reads
    #
    Return
    __________
    pandas Series of seqeuncing result:
    sensor_mutant reads
    """
    popab_Series=pop_Series/float(np.sum(pop_Series))
    sensor_Lst=popab_Series.index.tolist()
    seq_Dic={sensor:np.random.poisson(popab_Series[sensor]*total_reads*10**(np.random.uniform(-0.5,0.5,1)[0]),size=1)[0] for sensor in sensor_Lst}
    return pd.Series(seq_Dic)

# //////////////////////////////////////////
def htsensorsimu_generate_main(args):
    """
    Main entry for mageck count module
    """
    """
    get input parameters
    """
    # Number of mutant in the library
    lib_size=args.lib_size
    #distribution of Log10u during optimization
    distribution=args.distribution.split(',')
    range_Log10u=[]
    range_sigma=[]
    if len(distribution)!=4:
        logging.error('incorrect Log10u setting!')
        sys.exit(-1)
    else:
        range_Log10u=[float(distribution[0]),float(distribution[1])]
        range_sigma=[float(distribution[2]),float(distribution[3])]
    # output directory
    list_files=args.output_prefix.split('/')
    output_dir=''
    for i in range(len(list_files)-1):
        output_dir=output_dir+list_files[i]+'/'
    os.system('mkdir -p %s' %(output_dir))
    # Sub configure file to define the experimental setting
    exp_DF=pd.read_csv(filepath_or_buffer=args.exp_con,sep='\t',index_col='Bin')
    exp_Lst=exp_DF.columns.tolist()
    exp_Lst.sort()
    bins_Lst=exp_DF.index.tolist()
    bins_Lst.sort()
    #
    # Sub configure file to define the boundaries between bins for each sorting experiemnt
    # columns = ligand1R1 ligand1R2 ligand2R1 ligand2R2 ... (experiment condition)
    # rows = bins of sorting
    # values = [LogA0,LogA1]
    Bin_bou_DF=pd.read_csv(filepath_or_buffer=args.bin_bou_con,sep='\t',index_col='Bin')
    Bin_bou_Dic=Bin_bou_DF.to_dict()
    for exp in exp_Lst:
        for bins in bins_Lst:
            temp=Bin_bou_DF.loc[bins,exp].split(',')
            Bin_bou_Dic[exp][bins]=[float(temp[0]),float(temp[1])]
    Bin_bou_DF=pd.DataFrame(Bin_bou_Dic)
    # total reads of every sequencing library
    total_reads=args.total_reads
    # total cells of sorting experiment
    total_cells=args.total_cells
    # read count threshold for sensor elimination in initial library
    readth=args.readth
    """
    Generate a random library
    """
    my_lib_DF=construct_library(lib_size,range_Log10u,range_sigma)
    """
    initial library abundance, with those mutants with reads above threshold
    """
    ini_ab=my_lib_DF['Abundance']
    all_iniab_Series=sequencing(ini_ab, total_reads)
    count=np.sum(all_iniab_Series)
    iniab_Dic=(all_iniab_Series[all_iniab_Series>=readth]/float(count)).to_dict()
    print ('initial library constructed!')
    """
    cytometry
    """
    (sorting_Dic,bin_occ_DF,cell_con_DF)=cytometry_sorting(my_lib_DF,exp_DF,Bin_bou_DF,total_cells)
    print ('cytometry experiment finalized!')
    """
    Sequencing for each cytometry library and construct ctab: dict {sensor:DF}
    DF:
    columns = exp
    rows = bins
    value = read count
    """
    total_read_Dic={exp:{bins:0 for bins in bins_Lst} for exp in exp_Lst}
    ctab={sensor:{exp:{bins:0 for bins in bins_Lst} for exp in exp_Lst} for sensor in iniab_Dic}
    for exp in exp_Lst:
        for bins in bins_Lst:
            thisLib=sequencing(sorting_Dic[exp][bins] , total_reads)
            allLibsensor=thisLib.index.tolist()
            for sensor in allLibsensor:
                if sensor in ctab:
                    ctab[sensor][exp][bins]+=thisLib[sensor]
            total_read_Dic[exp][bins]=np.sum(thisLib)
    for (sensor,dict0) in ctab.items():
        ctab[sensor]=pd.DataFrame(dict0)
    total_read_DF=pd.DataFrame(total_read_Dic)
    print ('NGS finished!')
    # save some Python objects as pickle files
    pickle.dump(iniab_Dic,open('iniab_Dic.pickle','wb'))
    pickle.dump(ctab,open('sensor_read_Dic.pickle','wb'))
    pickle.dump(total_read_DF,open('total_read_DF.pickle','wb'))
    my_lib_DF.to_csv('my_library.txt',sep='\t')
    bin_occ_DF.to_csv('bin_occupation_configure.txt',sep='\t')
    cell_con_DF.to_csv('cell_count_configure.txt',sep='\t')

if __name__ == '__main__':
    try:
        args=htsensorsimu_generate_parseargs()
        htsensorsimu_generate_main(args)
    except KeyboardInterrupt:
        sys.stderr.write("Interrupted.\n")
        sys.exit(0)

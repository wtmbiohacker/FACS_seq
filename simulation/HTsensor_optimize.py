# this script read in the readcount.csv file; and loop to find the best parameter set for each sequence

from __future__ import print_function
import os
import sys
import numpy as np
import pickle
import scipy.stats
import scipy.optimize as optimize
import math
import pandas as pd
import logging
import argparse
import matplotlib.pyplot as plt
from random import shuffle

def htsensoropt_parseargs():
    """
    Parse arguments. Only used when htsensoropt.py is executed directly.
    """
    parser=argparse.ArgumentParser(description='Collecting read counts for multiple samples.')
    parser.add_argument('-i', '--ini_ab', required=True,help='pickle file of a Python dictionary about the abundance of sensor mutants in the initial library')
    parser.add_argument('-c', '--ctab', required=True,help='pickle file of a Python dictionary about the read count of each sensor mutant across all experiments and bins')
    parser.add_argument('-t', '--total', required=True,help='pickle file of a Python pandas DataFrame about the total read count of all seqeuncing libraries across all experiments and bins')
    parser.add_argument('-e','--exp_con',required=True,help='A file containing the design of the experiment, treatment of ligand concentration as well as the bins of cell sorting. Support file format: csv with tab delimiter.')
    parser.add_argument('-m','--cell_con',required=True,help='A file containing the cell number sorted into each bin in the original experiment. Support file format: csv with tab delimiter.')
    parser.add_argument('-o','--bin_occ_con',required=True,help='A file containing the occupation of each bin in cell sorting experiment. Support file format: csv with tab delimiter.')
    parser.add_argument('-b','--bin_bou_con',required=True,help='A file containing the boundary of each bin in cell sorting experiment. Support file format: csv with tab delimiter.')
    parser.add_argument('-r','--rrange',required=True,help='string to define the search range of Log10u and sigma. Format:lower,upper(Log10u);lower,upper(sigma)')
    parser.add_argument('-s','--search', type=int,default=20,help='number of points in the grid space of each dimension during optimization')
    parser.add_argument('-k','--cell_th', type=int,default=10,help='Minimum total number of cell for one sensor mutant across all bins of one sorting experiment, the mutant below this threshold will be excluded from further analysis.')
    parser.add_argument('-n','--output_prefix',default='sample1',help='The prefix of the output file(s). Default sample1.')
    args=parser.parse_args()
    return args

# ///////////////////////////////////////////////////////
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


def Cal_negativelnP(Log10u,sigma,BinBoundary_Lst,Readjsensor,Readj,Pj,Psensor):
    """
    Funciton used to return negative log10P for a sensor at one particular concentration within a relevant bin (given theoretical read count and actual read count, assume Poisson distribution)
    note that j=bin
    Parameters
    __________
    Log10u, mean of the normal distribution
    sigma, deviation of the normal distribution
    BinBoundary_Lst=[A0,A1], the defined X axis space of calculation
    Based on these three parameters, we have pjsensor, the theoretical probability of cells carrying this sensor mutant sorted into this bin
    Readjsensor, observed read count at this bin
    Readj, total read count of library belonging to corresponding bin
    Pj, relative abundance of cells in this bin against the whole population in one sorting experiment
    Psensor, relative abundance of cells carrying this sensor mutant in the whole population
    In our model, we have:
    rjsensor=Readj*Psensor*pjsensor/Pj, the lambda value of the Poisson distribution, based on which we calculate the probability to observe actually Readjsensor count by Poisson(lambda=rjsensor, k=Readjsensor)
    Return
    __________
    -log(Poisson(lambda=rjsensor, k=Readjsensor),10)
    """
    pjsensor=Cal_area(Log10u, sigma, BinBoundary_Lst)
    poisson_lambda=float(Readj)*float(Psensor)*float(pjsensor)/float(Pj)
    if poisson_lambda==0.0:
        poisson_lambda=1.0/(10.0**300)
    P=scipy.stats.poisson(poisson_lambda).pmf(int(Readjsensor))
    if P==0.0:
        return 350.0
    else:
        return -math.log(P,10)


def ObjectiveF(tobeop, *paras):
    """
    Function to calculate the -sumlnP for one sensor mutant at one ligand concentration with read numbers observed at all bins, this function accepts Log10u0, NGS and Binflu DataFrame (below), Psensor, sigma, Bin_bou_DF (below)
    Parameters
    __________
    tobeop: to be optimized, tuple (Log10u, sigma)
    *paras: tupled parameters, including:
    Psensor, relative abundance of cells carrying this sensor mutant in the whole population
    DF: pandas dataframe of parameters in different bins
    columns = 'Readjsensor' 'Readj' 'Pj' 'BinBoundary_Lst', the meaning of which, see the annotation in function Cal_negativelnP
    rows = bins of sorting
    values
    Return
    __________
    Summary of sum(-log(P,10)) given the Log10u and sigma across all bins defined in DF
    
    """
    sum_negativelnP=0.0
    Log10u, sigma = tobeop
    DF, Psensor = paras
    for index, bin in DF.iterrows():
        Readjsensor=bin['Readjsensor']
        Readj=bin['Readj']
        Pj=bin['Pj']
        BinBoundary_Lst=bin['BinBoundary_Lst']
        sum_negativelnP+=Cal_negativelnP(Log10u,sigma,BinBoundary_Lst,Readjsensor,Readj,Pj,Psensor)
    return sum_negativelnP

# //////////////////////////////////////////
def htsensoropt_main(args):
    """
    Main entry for mageck count module
    """
    """
    get input parameters
    """
    # seq abundance in initial concentration
    # {'seq0':float,'seq1':float, ... ,'seqk':float}
    iniab_Dicpickle=args.ini_ab
    iniab_Dic=pickle.load(open(iniab_Dicpickle,'rb'))
    #
    # read count table in Dic format
    # each node in one dic point to a DataFramge
    # {'seq0':DF,'seq1':DF, ..., 'seqk':DF, ...}
    # DF structure:
    # columns = ligand1R1 ligand1R2 ligand2R1 ligand2R2 ... (experiment condition)
    # rows = bins of sorting
    # values = read count
    sensor_read_Dicpickle=args.ctab
    sensor_read_Dic=pickle.load(open(sensor_read_Dicpickle,'rb'))
    sensor_Lst=sensor_read_Dic.keys()
    sensor_Lst.sort()
    #shuffle(sensor_Lst)
    #
    # DataFrame of total read counts for each library (conibinj)
    # columns = ligand1R1 ligand1R2 ligand2R1 ligand2R2 ... (experiment condition)
    # rows = bins of sorting
    # values = read count
    total_read_DFpickle=args.total
    total_read_DF=pickle.load(open(total_read_DFpickle,'rb'))
    #
    # Sub configure file to define the experimental setting
    exp_DF=pd.read_csv(filepath_or_buffer=args.exp_con,sep='\t',index_col='Bin')
    exp_Lst=exp_DF.columns.tolist()
    bins_Lst=exp_DF.index.tolist()
    #
    # Sub configure file to define the cell number sorted into each bin in the original experiment
    cell_DF=pd.read_csv(filepath_or_buffer=args.cell_con,sep='\t',index_col='Bin')
    #
    # Sub configure file to define the Pij (occupation of each bin) for each sorting experiemnt
    # columns = ligand1R1 ligand1R2 ligand2R1 ligand2R2 ... (experiment condition)
    # rows = bins of sorting
    # values = occupation of each bin in its relevant experiment
    Bin_occ_DF=pd.read_csv(filepath_or_buffer=args.bin_occ_con,sep='\t',index_col='Bin')
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
    #
    #search range of Log10u and sigma during optimization
    rrange=args.rrange.split(',')
    range_Log10u=[]
    range_sigma=[]
    if len(rrange)!=4:
        logging.error('incorrect Log10u setting!')
        sys.exit(-1)
    else:
        range_Log10u=[float(rrange[0]),float(rrange[1])]
        range_sigma=[float(rrange[2]),float(rrange[3])]
    #
    # output directory
    list_files=args.output_prefix.split('/')
    output_dir=''
    for i in range(len(list_files)-1):
        output_dir=output_dir+list_files[i]+'/'
    os.system('mkdir -p %s' %(output_dir))
    os.system('mkdir -p %sheatmap/' %(output_dir))
    os.system('mkdir -p %scell_count/'%(output_dir))
    # search calculation intensity
    # number of trial for each parameter in a constructed grid space within the search range
    search_number=args.search
    #
    # construct the 2D grid space
    rrange = ((range_Log10u[0], range_Log10u[1]), (range_sigma[0], range_sigma[1]))
    #
    # cell count threshold to exclude sensor mutant with total cell count below in one sorting experiment
    cellth=args.cell_th
    os.system('cat /dev/null > %seliminated_cell_exp.txt'%(output_dir))
    eliminated_list=open('%seliminated_cell_exp.txt'%(output_dir),'w')
    """
    Search for the peak at the grid 2D space of Log10u and sigma
    """
    # dict to store the optimization process
    sensor_Log10uDic={}
    sensor_sigmaDic={}
    sensor_negLog10PDic={}
    for sensor in sensor_Lst:
        Psensor=iniab_Dic[sensor]
        sensor_Log10uDic[sensor]={}
        sensor_sigmaDic[sensor]={}
        sensor_negLog10PDic[sensor]={}
        cellcountDic={}
        for exp in exp_Lst:
            sensor_Log10uDic[sensor][exp]=0.0
            sensor_sigmaDic[sensor][exp]=0.0
            sensor_negLog10PDic[sensor][exp]=0.0
            cellcountDic[exp]={}
            DF_dic={}
            for bins in bins_Lst:
                cell_count=cell_DF.loc[bins,exp]
                Readj=total_read_DF.loc[bins,exp]
                # check whether Readj>>Cellj
                if Readj>10*cell_count:
                    correction_factor=float(cell_count)/float(Readj)
                else:
                    correction_factor=1.0
                DF_dic[bins]={}
                DF_dic[bins]['Readj']=int(total_read_DF.loc[bins,exp]*correction_factor)
                DF_dic[bins]['Readjsensor']=int(sensor_read_Dic[sensor].loc[bins,exp]*correction_factor)
                cellcountDic[exp][bins]=DF_dic[bins]['Readjsensor']
                DF_dic[bins]['Pj']=Bin_occ_DF.loc[bins,exp]
                DF_dic[bins]['BinBoundary_Lst']=Bin_bou_DF.loc[bins,exp]
            DF=pd.DataFrame(DF_dic).T
            # cell count across all bins > a given threshold
            if np.sum(DF['Readjsensor'])>cellth:
                paras=(DF, Psensor)
                res = optimize.brute(ObjectiveF, rrange, args = paras, Ns = search_number, full_output = True, finish = optimize.fmin)
                # reject those results beyond the search range we defined due to the finish = optimize.fmin step
                if (res[0][0]<range_Log10u[0]) or (res[0][0]>range_Log10u[1]) or (res[0][1]<range_sigma[0]) or (res[0][1]>range_sigma[1]):
                    res = optimize.brute(ObjectiveF, rrange, args = paras, Ns = search_number*3, full_output = True, finish = None)
                sensor_Log10uDic[sensor][exp]=res[0][0]
                sensor_sigmaDic[sensor][exp]=res[0][1]
                sensor_negLog10PDic[sensor][exp]=res[1]
                print ('%s: %s result: %s %s'%(sensor,exp,res[0],res[1]))
                # res[0]: best (u); res[1]: best -sumlnP, res[2]: grid space array; res[3]: -sumlnP landscape array
                X= res[2][0]
                Y= res[2][1]
                Z= res[3]
                # to imporve the resolution of minimum -Log10P region in heatmap, we restructure this matrix by resetting each value as min(10*globalMin,value)
                Zmin=Z.min()
                Z[Z >10*Zmin] = 10*Zmin
                # store the optimization process of each condition as a 2D heatmap
                #plt.contourf(X,Y,Z)
                #plt.colorbar()
                #plt.xlabel('Log10u')
                #plt.ylabel('sigma')
                #plt.savefig('%sheatmap/%s_%s_opt.png'%(output_dir,sensor,exp))
                #plt.clf()
            else:
                sensor_Log10uDic[sensor][exp]=np.NaN
                sensor_sigmaDic[sensor][exp]=np.NaN
                sensor_negLog10PDic[sensor][exp]=np.NaN
                print('%s in %s has total cell count below threshold'%(sensor,exp))
                print('%s_%s'%(sensor,exp),file=eliminated_list)
        # write in the cell count of each bin of one particular sensor
        cellcountDF=pd.DataFrame(cellcountDic)
        cellcountDF.to_csv('%scell_count/%s.csv'%(output_dir,sensor),sep='\t')
    # pandas sort the df automatically during construction
    sensor_Log10uDF=pd.DataFrame(sensor_Log10uDic).T
    sensor_sigmaDF=pd.DataFrame(sensor_sigmaDic).T
    sensor_negLog10PDF=pd.DataFrame(sensor_negLog10PDic).T
    sensor_Log10uDF.index.name='Sensor'
    sensor_sigmaDF.index.name='Sensor'
    sensor_negLog10PDF.index.name='Sensor'
    # write to file
    # save some Python objects as pickle files
    sensor_Log10uDF.to_csv('%ssensor_Log10u.csv'%(output_dir),sep='\t')
    sensor_sigmaDF.to_csv('%ssensor_sigma.csv'%(output_dir),sep='\t')
    sensor_negLog10PDF.to_csv('%ssensor_negLog10P.csv'%(output_dir),sep='\t')
    eliminated_list.close()

if __name__ == '__main__':
    try:
        args=htsensoropt_parseargs()
        htsensoropt_main(args)
    except KeyboardInterrupt:
        sys.stderr.write("Interrupted.\n")
        sys.exit(0)


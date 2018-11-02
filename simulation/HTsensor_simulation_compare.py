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
from scipy.stats import pearsonr

def htsensorsimu_com_parseargs():
    """
    Parse arguments. Only used when htsensorsimu_com.py is executed directly.
    """
    parser=argparse.ArgumentParser(description='Collecting read counts for multiple samples.')
    parser.add_argument('-l','--my_lib', required=True,help='.csv format library raw data')
    parser.add_argument('-u','--log10u',required=True,help='.csv format inferred sensor mutant log10u data, one or multiple replicates')
    parser.add_argument('-d','--sigma',required=True,help='.csv format inferred sensor mutant sigma data, one or multiple replicates')
    parser.add_argument('-r','--rrange',required=True,help='string to define the search range of Log10u and sigma. Format:lower,upper(Log10u);lower,upper(sigma)')
    parser.add_argument('-o','--output_prefix',default='sample1',help='The prefix of the output file(s). Default sample1.')
    args=parser.parse_args()
    return args

# ///////////////////////////////////////////////////////
def htsensorsimu_com_main(args):
    """
    Main entry
    """
    """
    get input parameters
    """
    # simulated library
    my_lib_DF=pd.read_csv(filepath_or_buffer=args.my_lib,sep='\t',index_col='Sensor')
    # calculated Log10u for sensors in multiple replicates
    Log10u_rep_DF=pd.read_csv(filepath_or_buffer=args.log10u,sep='\t',index_col='Sensor')
    # calculated sigma for sensors in multiple replicates
    sigma_rep_DF=pd.read_csv(filepath_or_buffer=args.sigma,sep='\t',index_col='Sensor')
    #search range of Log10u and sigma during optimization
    rrange=args.rrange.split(',')
    range_Log10u=[]
    range_sigma=[]
    if len(rrange)!=4:
        logging.error('incorrect Log10u setting!')
        sys.exit(-1)
    else:
        range_Log10u=[float(rrange[0])-0.1,float(rrange[1])+0.1]
        range_sigma=[float(rrange[2])-0.05,float(rrange[3])+0.05]
    # output directory
    list_files=args.output_prefix.split('/')
    output_dir=''
    for i in range(len(list_files)-1):
        output_dir=output_dir+list_files[i]+'/'
    os.system('mkdir -p %s' %(output_dir))
    # remove sensor mutants filtered by the cell count threshold (with NaN in dataframe)
    Log10u_repdropna_DF=Log10u_rep_DF.dropna(axis=0, how='any')
    sigma_repdropna_DF=sigma_rep_DF.dropna(axis=0, how='any')
    # calculate the Log10u and sigma for each sensor mutant by avaraging multiple replicates
    sensor_Lst=Log10u_repdropna_DF.index.tolist()
    sensor_mean_Dic={sensor:{} for sensor in sensor_Lst}
    for sensor in sensor_Lst:
        sensor_mean_Dic[sensor]['Calculated Log10u']=np.mean(Log10u_repdropna_DF.loc[sensor])
        sensor_mean_Dic[sensor]['Calculated sigma']=np.mean(sigma_repdropna_DF.loc[sensor])
    sensor_mean_DF=pd.DataFrame(sensor_mean_Dic).T
    # mapping with my library to extract those passing initial lib read and sorting experiment cell number threshold
    my_final_lib_DF=pd.concat([my_lib_DF, sensor_mean_DF], axis=1, join='inner')
    # draw scatter plot for theoretical and calculated Log10u
    Xaxis=np.array(my_final_lib_DF['Log10u'])
    Yaxis=np.array(my_final_lib_DF['Calculated Log10u'])
    pearsonCE,pValue=pearsonr(Xaxis,Yaxis)
    plt.scatter(Xaxis,Yaxis,s=30,color='#5DADE2')
    plt.xlabel('Log10u')
    plt.ylabel('Calculated Log10u')
    XYmin=min(min(Xaxis),min(Yaxis))-0.1*(range_Log10u[1]-range_Log10u[0])
    XYmax=max(max(Xaxis),max(Yaxis))+0.1*(range_Log10u[1]-range_Log10u[0])
    plt.xlim(XYmin,XYmax)
    plt.ylim(XYmin,XYmax)
    plt.text(XYmin+0.02*(XYmax-XYmin),XYmax-0.05*(XYmax-XYmin),'Pearson Correlation Coefficient=%s; n=%d'%(("{0:.3f}".format(pearsonCE)),len(Xaxis)),fontsize=12)
    plt.savefig('%sLog10u_comparison.png'%(output_dir))
    plt.clf()
    # draw scatter plot for theoretical and calculated sigma
    Xaxis=np.array(my_final_lib_DF['sigma'])
    Yaxis=np.array(my_final_lib_DF['Calculated sigma'])
    pearsonCE,pValue=pearsonr(Xaxis,Yaxis)
    plt.scatter(Xaxis,Yaxis,s=30,color='#5DADE2')
    plt.xlabel('sigma')
    plt.ylabel('Calculated sigma')
    XYmin=min(min(Xaxis),min(Yaxis))-0.1*(range_sigma[1]-range_sigma[0])
    XYmax=max(max(Xaxis),max(Yaxis))+0.1*(range_sigma[1]-range_sigma[0])
    plt.xlim(XYmin,XYmax)
    plt.ylim(XYmin,XYmax)
    plt.text(XYmin+0.02*(XYmax-XYmin),XYmax-0.05*(XYmax-XYmin),'Pearson Correlation Coefficient=%s; n=%d'%(("{0:.3f}".format(pearsonCE)),len(Xaxis)),fontsize=12)
    plt.savefig('%ssigma_comparison.png'%(output_dir))
    plt.clf()
    my_final_lib_DF.to_csv('%ssimulation_data.csv'%(output_dir),sep='\t')

if __name__ == '__main__':
    try:
        args=htsensorsimu_com_parseargs()
        htsensorsimu_com_main(args)
    except KeyboardInterrupt:
        sys.stderr.write("Interrupted.\n")
        sys.exit(0)

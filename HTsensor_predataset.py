# Prepare dataset for Bayesian calculation based on the tables derived from the .fq raw data

from __future__ import print_function
import os
import sys
import pandas as pd
import numpy as np
import logging
import argparse
import pickle

def htsensorpre_parseargs():
    """
    Parse arguments.
    """
    parser=argparse.ArgumentParser(description='Collecting read counts for multiple samples.')
    parser.add_argument('-c','--ctab',required=True,help='A file containing the read count information of all fastq files. Support file format: csv and txt.')
    parser.add_argument('-u','--unmap',required=True,help='A file containing the unmapped read count information of all fastq files. Support file format: csv and txt.')
    parser.add_argument('-r', '--readth', type=int, help='read threshold in initial library to exclude sensor mutant for further analysis')
    parser.add_argument('-e','--exp_con',required=True,help='A file containing the design of the experiment, treatment of ligand concentration as well as the bins of cell sorting. Support file format: csv and txt.')
    parser.add_argument('-i','--ini_lib',required=True,help='The sample label of initial library without soring, used to exclude sensor mutants without enough reads from further analysis')
    parser.add_argument('-n','--output_prefix',default='sample1',help='The prefix of the output file(s). Default sample1.')
    args=parser.parse_args()
    return args


def get_the_sum(ctab,unmap,iniLib,exp_configure):
    """
    calculate the sum of read counts for each library
    Parameters
    __________
    ctab, pandas DataFrame object
    columns = lib1 lib2 lib3 ...
    rows = sensor mutants
    values = read count
    ...
    #
    unmap, pandas DataFrame object,unmapped read count across all libraries
    columns = lib1 lib2 lib3 ...
    rows = unmapped read
    values = read count
    ...
    #
    iniLib defines the initial library derived from the unsorted population
    #
    exp_configure, pandas DataFrame of experiment setting
    columns = ligand1R1 ligand1R2 ligand2R1 ligand2R2 ... (experiment condition)
    rows = bins of sorting
    values = library label
    Return
    __________
    pandas DataFrame of total read counts for each library (conibinj) except for the initial library
    columns = ligand1R1 ligand1R2 ligand2R1 ligand2R2 ... (experiment condition)
    rows = bins of sorting
    values = total read count

    """
    libLst=ctab.columns.tolist()
    if libLst!=unmap.columns.tolist():
        logging.error('Different library setting between mapped and unmapped read files!')
        sys.exit(-1)
    libLst.remove(iniLib)
    expLst=exp_configure.columns.tolist()
    binsLst=exp_configure.index.tolist()
    total_count_Dic={lib:(np.sum(ctab[lib])+np.sum(unmap[lib])) for lib in libLst}
    tr_expbins_Dic={exp:{bins:total_count_Dic[exp_configure.loc[bins,exp]] for bins in binsLst} for exp in expLst}
    return pd.DataFrame(tr_expbins_Dic)


def initial_ab(ctab,unmap,iniLib,readth):
    """
    calculate the abundance of every sensor in initial library
    remove those below the read threshold
    Parameters
    __________
    ctab, pandas DataFrame object,mapped read count of each sensor mutant across all libraries
    columns = lib1 lib2 lib3 ...
    rows = sensor mutants
    values = read count
    ...
    #
    unmap, pandas DataFrame object,unmapped read count across all libraries
    columns = lib1 lib2 lib3 ...
    rows = unmapped read
    values = read count
    ...
    iniLib defines the initial library derived from the unsorted population
    #
    readth defined the threshold of read count below which the sensor mutant in the initial library will be excluded from further analysis
    Return
    __________
    Dictionary of sensor (those passing the cutoff) and abundance in initial library
    {sensor:ab}
    List of sensor (fail to pass)
    (Dic, Lst)
    """
    total_count=np.sum(ctab[iniLib])+np.sum(unmap[iniLib])
    iniab_Dic={}
    eliminateLst=[]
    for (sensor,read) in ctab[iniLib].to_dict().items():
        if read>=readth:
            iniab_Dic[sensor]=float(read)/float(total_count)
        else:
            eliminateLst.append(sensor)
    return (iniab_Dic, eliminateLst)


def include_sensor(ctab,iniLib,exp_configure,sensorLst):
    """
    convert the table of sensor mutant read count of each library to new table, associating library with relevant experiment condition and binss:
    Parameters
    __________
    ctab, pandas DataFrame object
    columns = lib1 lib2 lib3 ...
    rows = sensor mutants
    values = read count
    ...
    #
    iniLib defines the initial library derived from the unsorted population
    #
    exp_configure, pandas DataFrame of experiment setting
    columns = ligand1R1 ligand1R2 ligand2R1 ligand2R2 ... (experiment condition)
    rows = bins of sorting
    values = library label
    #
    sensorLst
    sensors passing the cutoff of initial library read count
    Return
    __________
    Dictionary of sensor (those passing the cutoff according to the sensorLst) and relevant read count pandas DataFrame
    {sensor:df}
    the strucuture of df
    columns = ligand1R1 ligand1R2 ligand2R1 ligand2R2 ... (experiment condition)
    rows = bins of sorting
    values = read count of sensor mutant here
    """
    sensor_read_Dic={}
    libLst=ctab.columns.tolist()
    libLst.remove(iniLib)
    expLst=exp_configure.columns.tolist()
    binsLst=exp_configure.index.tolist()
    for sensor in sensorLst:
        thisDic={exp:{bins:ctab.loc[sensor][exp_configure.loc[bins,exp]] for bins in binsLst} for exp in expLst}
        sensor_read_Dic[sensor]=pd.DataFrame(thisDic)
    return sensor_read_Dic


def htsensorpre_main(args):
    """
    Main entry
    """
    # Sub configure file to define the experimental setting
    exp_configure=pd.read_csv(filepath_or_buffer=args.exp_con,sep='\t',index_col='Bin')
    # Read count table of all .fq raw data
    ctab=pd.read_csv(filepath_or_buffer=args.ctab,sep='\t',index_col='sensor')
    # Unmapped read count table of all .fq raw data
    unmap=pd.read_csv(filepath_or_buffer=args.unmap,sep='\t',index_col='unmapped read')
    # read count threshold for sensor elimination in initial library
    readth=args.readth
    # label to indicate the initial library
    iniLib=args.ini_lib
    #
    # calculation
    total_read_DF=get_the_sum(ctab,unmap,iniLib,exp_configure)
    (iniab_Dic, eliminateLst)=initial_ab(ctab,unmap,iniLib,readth)
    sensor_read_Dic=include_sensor(ctab,iniLib,exp_configure,iniab_Dic.keys())
    iniab_ss=pd.Series(iniab_Dic)
    #
    # write to file
    # save some Python objects as pickle files
    list_files=args.output_prefix.split('/')
    output_dir=''
    for i in range(len(list_files)-1):
        output_dir=output_dir+list_files[i]+'/'
    os.system('mkdir -p %s' %(output_dir))
    fl=args.output_prefix+'.eliminate.txt'
    os.system('cat /dev/null > %s'%(fl))
    ofilel=open(fl,'w')
    ofilel.write('#Read threshold: %d\n'%(readth))
    for sensor in eliminateLst:
        ofilel.write('%s\n'%(sensor))
    iniab_ss.to_csv('%sinitial_ab.csv'%(output_dir),sep='\t')
    ofilel.close()
    # write in all sensor read count table as .csv
    os.system('mkdir -p %ssensor_ctab/'%(output_dir))
    for (sensor,ctab_df) in sensor_read_Dic.items():
        ctab_df.to_csv('%ssensor_ctab/%s.csv'%(output_dir,sensor),sep='\t')
    pickle.dump(total_read_DF,open('total_read_DF.pickle','wb'))
    pickle.dump(iniab_Dic,open('iniab_Dic.pickle','wb'))
    pickle.dump(sensor_read_Dic,open('sensor_read_Dic.pickle','wb'))
    os.system('cp *.pickle %s'%(output_dir))
    return 0


if __name__ == '__main__':
    try:
        args=htsensorpre_parseargs()
        htsensorpre_main(args)
    except KeyboardInterrupt:
        sys.stderr.write("Interrupted.\n")
        sys.exit(0)

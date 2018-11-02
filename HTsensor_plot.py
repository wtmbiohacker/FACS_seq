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
from scipy.stats import mannwhitneyu
from Bio.Data import CodonTable
standard_table = CodonTable.unambiguous_dna_by_id[1]
plt.rcParams["font.family"] = "Arial"

def htsensor_plot_parseargs():
    """
    Parse arguments. Only used when htsensor_plot.py is executed directly.
    """
    parser=argparse.ArgumentParser(description='Collecting read counts for multiple samples.')
    parser.add_argument('-w','--wt_seq', default='ATGAATATCTTACATATATGTGTGACCTCAAAATGGTTCAATATTGACAACAAAATTGTCGATCACCGCCCTTGA',help='wild type sequences')
    parser.add_argument('-a','--aa_range', default='2,24',help='start and end of variable region')
    parser.add_argument('-n','--wt_name', default='TnaC',help='Wild type name of the sequence')
    parser.add_argument('-c','--replicate_con', required=True,help='configure file of replicates')
    parser.add_argument('-r','--rrange',required=True,help='string to define the search range of Log10u and sigma. Format:lower,upper(Log10u);lower,upper(sigma)')
    parser.add_argument('-u','--log10u',required=True,help='.csv format inferred sensor mutant log10u data, one or multiple replicates')
    parser.add_argument('-d','--sigma',required=True,help='.csv format inferred sensor mutant sigma data, one or multiple replicates')
    parser.add_argument('-o','--output_prefix',default='sample1',help='The prefix of the output file(s). Default sample1.')
    args=parser.parse_args()
    return args

# ///////////////////////////////////////////////////////
def scatter_plot(Xaxis,Yaxis,xlabel,ylabel,output_dir,name,Xlen,Ylen):
    """
    funciton to draw scatter plot
    Parameters
    ________________
    numpy array: Xaixs, Yaxis
    output_dir, name
    Xlen, Ylen, X,Y lengths of figures
    """
    plt.figure(figsize=(Xlen,Ylen))
    pearsonCE,pValue=pearsonr(Xaxis,Yaxis)
    plt.scatter(Xaxis, Yaxis, s=30, alpha=0.5, edgecolors='#2E86C1', facecolors='#FDFEFE')
    #plt.scatter(Xaxis,Yaxis,s=30,color='#5DADE2')
    XYmin=min(min(Xaxis),min(Yaxis))
    XYmax=max(max(Xaxis),max(Yaxis))
    temp=0.1*(XYmax-XYmin)
    XYmin=XYmin-temp
    XYmax=XYmax+temp
    plt.xlim(XYmin,XYmax)
    plt.ylim(XYmin,XYmax)
    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95)
    #plt.text(XYmin+0.02*(XYmax-XYmin),XYmax-0.05*(XYmax-XYmin),'Pearson Correlation Coefficient=%s; n=%d'%(("{0:.3f}".format(pearsonCE)),len(Xaxis)),fontsize=12)
    plt.savefig('%s%s.png'%(output_dir,name),dpi=400)
    plt.clf()
    return 0

def hist_plot(data,where,bin_num,xlabel,output_dir,name):
    """
    funciton to draw histogram plot
    Parameters
    ________________
    numpy array: data
    where, tuple (lower, upper)
    """
    lower,upper=where
    bins=np.linspace(lower,upper,bin_num)
    plt.hist(data,bins,range=(np.nanmin(data),np.nanmax(data)))
    plt.xlabel(xlabel)
    plt.savefig('%s%s.png'%(output_dir,name),dpi=400)
    plt.clf()
    return 0

def heatmap_plot(data,cmap,vmin,vmax,xlabel,ylabel,xticks,yticks,tick_fontsize,output_dir,name,wt_coordinate,Xlen,Ylen,title,aspect,shrink):
    """
    function to draw heatmap
    Parameters
    ________________________
    data: np array as the raw data for heatmap
    vmin and vmax: minimun and maxium of value along Z axis in heatmap
    xlabel and ylabel: label of axis
    xticks and yticks: ticks of axis
    name: saved figure file
    wt_coordinate: wt codon or aa in heatmap to be highlighted
    Xlen, Ylen: X,Y lengths of figures
    title: colorbar legend
    aspect: length:width ratio of colorbar
    shrink: occupying ratio of colorbar
    """
    data_copy=data
    # mask the NaN cell
    data_copy=np.ma.masked_invalid(data_copy)
    fig, ax = plt.subplots(figsize=(Xlen,Ylen))
    heatmap = ax.matshow(data_copy, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.patch.set(hatch='///////', edgecolor='green')
    cbar = fig.colorbar(heatmap, aspect=aspect, shrink=shrink)
    cbar.ax.set_ylabel(title,fontsize=18)
    for (xi,yi), in zip(wt_coordinate[['Xaxis','Yaxis']].values):
        rec = plt.Rectangle((xi-0.5,yi-0.5),1.0,1.0, fill=False, edgecolor='black', lw=2 )
        ax.add_artist(rec)
    # add the ticks
    plt.xticks(range(0,len(xticks)),xticks,fontsize=tick_fontsize)
    plt.yticks(range(0,len(yticks)),yticks,fontsize=tick_fontsize)
    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95)
    plt.savefig('%s%s.png'%(output_dir,name),dpi=800)
    plt.close(fig)
    return 0

def codon_to_aa(codon_df,column_name,aa_order,aa_range):
    """
    core function of this script, used to convert a matrix (dataframe) with codon as index to matrix with aa as index (combine codon if degenerate exist)
    Parameters
    ____________________
    codon_df: pandas DataFrame
    columns: conditionXR1 conditionXR2 ... conditionXmean Position codon aa wt_aa whetherWT, etc
    rows: codon_level mutants
    column_name: column to be averaged for degenerate codon
    aa_order: aa ordered to sort the returned DF (Y axis, index)
    aa_range: position fo the returned DF (X axis, columns)
    Return
    ___________________
    DF with X axis: aa_range, Yaxis: aa_order, values: average value of column_name parameter of degenerate codons
    """
    aa_Dic={position:{aa:0.0 for aa in aa_order} for position in aa_range}
    for position in aa_range:
        for aa in aa_order:
            aa_Dic[position][aa]=np.mean(codon_df[(codon_df['aa']==aa)&(codon_df['Position']==position)][column_name])
    aa_DF=pd.DataFrame(aa_Dic)
    return aa_DF.reindex(aa_order)

def codon_to_aa_noise(mean_df,sigma_df,column_name,aa_order,aa_range):
    """
    core function of this script, used to convert a matrix (dataframe) with codon as index to matrix with aa as index (combine codon if degenerate exist, calculate the difference between degenerate codons)
    Parameters
    ____________________
    mean(sigma)_df: pandas DataFrame of mean(sigma)
    columns: conditionXR1 conditionXR2 ... conditionXmean Position codon aa wt_aa whetherWT, etc
    rows: codon_level mutants
    column_name: column of mean for this parameter
    aa_order: aa ordered to sort the returned DF (Y axis, index)
    aa_range: position fo the returned DF (X axis, columns)
    Return
    ___________________
    For aa with degenerate codons: return -log10P p from scipy.stats.f_oneway(np.random.normal(mean1,sigma1,20);..;..)
    DF with X axis: aa_range, Yaxis: aa_order, values: -log10P
    """
    anova_number=10
    aa_Dic={position:{aa:0.0 for aa in aa_order} for position in aa_range}
    for position in aa_range:
        for aa in aa_order:
            mean_Dic=mean_df[(mean_df['aa']==aa)&(mean_df['Position']==position)].T.to_dict()
            sigma_Dic=sigma_df[(sigma_df['aa']==aa)&(sigma_df['Position']==position)].T.to_dict()
            array_Dic={}
            for ((k1,v1),(k2,v2)) in zip(mean_Dic.items(),sigma_Dic.items()):
                if np.isnan(v1[column_name])==False and np.isnan(v2[column_name])==False:
                    array_Dic[k1]=np.random.normal(v1[column_name],v2[column_name],anova_number)
            if len(array_Dic)==0:
                aa_Dic[position][aa]=np.NaN
            elif len(array_Dic)==1:
                aa_Dic[position][aa]=0
            else:
                (F,p)=scipy.stats.f_oneway(*array_Dic.values())
                aa_Dic[position][aa]=-math.log(max(p,10**(-300)),10)
    aa_DF=pd.DataFrame(aa_Dic)
    return aa_DF.reindex(aa_order)

def codon_to_codon(codon_df,column_name,codon_order,codon_range):
    """
    core function of this script, used to convert a matrix (dataframe) with codon as index to matrix with codon as index (do not combine codon if degenerate exist)
    Parameters
    ____________________
    codon_df: pandas DataFrame
    columns: conditionXR1 conditionXR2 ... conditionXmean Position codon codon wt_aa whetherWT, etc
    rows: codon_level mutants
    column_name: column to be averaged for degenerate codon
    codon_order: codon ordered to sort the returned DF (Y axis, index)
    codon_range: position fo the returned DF (X axis, columns)
    Return
    ___________________
    DF with X axis: codon_range, Yaxis: codon_order, values: average value of column_name parameter of degenerate codons
    """
    codon_Dic={position:{codon:0.0 for codon in codon_order} for position in codon_range}
    for position in codon_range:
        for codon in codon_order:
            codon_Dic[position][codon]=np.mean(codon_df[(codon_df['codon']==codon)&(codon_df['Position']==position)][column_name])
    codon_DF=pd.DataFrame(codon_Dic)
    return codon_DF.reindex(codon_order)

def compare_syn_un(codon_df,column_name,aa_range,name,output_dir):
    """
    organize the dataset into two parts, mutant1 vs mutant2: classified as synonymous and unsynonymous mutations, and calculate the metrics
    Parameters
    ____________________
    codon_df: pandas DataFrame
    columns: conditionXR1 conditionXR2 ... conditionXmean Position codon aa wt_aa whetherWT, etc
    rows: codon_level mutants
    column_name: column of mean for this parameter
    aa_range: position fo the returned DF (X axis, columns)
    name: name of the saved .csv and .png files
    output_dir: name of the dierectory to save the output
    Return
    ___________________
    saved .png(histogram) files
    return 0
    """
    syn_Lst=[]
    unsyn_Lst=[]
    for position in aa_range:
        this_DF=codon_df[(codon_df['Position']==position)]
        sensor_Lst=this_DF.index.tolist()
        for (x,sensorX) in enumerate(sensor_Lst):
            for (y,sensorY) in enumerate(sensor_Lst):
                if y>x:
                    aaX=this_DF.loc[sensorX,'aa']
                    aaY=this_DF.loc[sensorY,'aa']
                    wt_codonX=this_DF.loc[sensorX,'whetherWTcodon']
                    wt_codonY=this_DF.loc[sensorY,'whetherWTcodon']
                    diff=this_DF.loc[sensorX,column_name]-this_DF.loc[sensorY,column_name]
                    if wt_codonX!='WT' or wt_codonY!='WT':
                        if aaX==aaY:
                            syn_Lst.append(diff)
                        else:
                            unsyn_Lst.append(diff)
    vmin=min(np.nanmin(syn_Lst),np.nanmin(unsyn_Lst))
    vmax=max(np.nanmax(syn_Lst),np.nanmax(unsyn_Lst))
    bins=np.linspace(vmin,vmax,50)
    U,p=mannwhitneyu(np.array(syn_Lst),np.array(unsyn_Lst),True,'two-sided')
    plt.hist(unsyn_Lst,bins,range=(vmin,vmax),label='Unsynonymous codon',alpha=0.3)
    plt.hist(syn_Lst,bins,range=(vmin,vmax),label='Synonymous codon',alpha=0.3)
    plt.xlim([vmin,vmax])
    plt.legend(loc='upper right')
    plt.xlabel('Difference of Log10u between mutant pairs at %s'%(column_name.lower()))
    plt.ylabel('Number of mutant pairs')
    plt.annotate('Mann Whitney U test\nLog10(P) value=%s'%(("{0:.3f}".format(math.log(p,10)))),xy=(0.02, 0.85), xycoords='axes fraction',fontsize=12)
    plt.savefig('%s%s.png'%(output_dir,name),dpi=400)
    plt.clf()
    return 0

def htsensor_plot_main(args):
    """
    Main entry
    """
    """
    input parameters
    """
    # order of amino acid according to their chemistry
    aa_order=['H','K','R','D','E','C','M','N','Q','S','T','A','G','I','L','P','V','F','W','Y','*']
    # order of codon according to their amino acid chemistry
    codon_order=['CAT_H','CAC_H','AAA_K','AAG_K','CGC_R','AGA_R','AGG_R','CGA_R','CGG_R','CGT_R','GAC_D','GAT_D','GAA_E','GAG_E','TGT_C','TGC_C','ATG_M','AAT_N','AAC_N','CAA_Q','CAG_Q','TCA_S','AGC_S','AGT_S','TCC_S','TCG_S','TCT_S','ACC_T','ACA_T','ACG_T','ACT_T','GCA_A','GCC_A','GCG_A','GCT_A','GGA_G','GGC_G','GGG_G','GGT_G','ATA_I','ATC_I','ATT_I','TTA_L','CTA_L','CTC_L','CTG_L','CTT_L','TTG_L','CCT_P','CCA_P','CCC_P','CCG_P','GTG_V','GTC_V','GTA_V','GTT_V','TTC_F','TTT_F','TGG_W','TAC_Y','TAT_Y','TAA_*','TAG_*','TGA_*']
    codon_order_Lst=['CAT','CAC','AAA','AAG','CGC','AGA','AGG','CGA','CGG','CGT','GAC','GAT','GAA','GAG','TGT','TGC','ATG','AAT','AAC','CAA','CAG','TCA','AGC','AGT','TCC','TCG','TCT','ACC','ACA','ACG','ACT','GCA','GCC','GCG','GCT','GGA','GGC','GGG','GGT','ATA','ATC','ATT','TTA','CTA','CTC','CTG','CTT','TTG','CCT','CCA','CCC','CCG','GTG','GTC','GTA','GTT','TTC','TTT','TGG','TAC','TAT','TAA','TAG','TGA']
    # wild type sequence
    wt_seq=args.wt_seq
    # start and end aa of variable region in wild type sequence
    aa_range=[int(item) for item in args.aa_range.split(',')]
    aa_range_Lst=[i for i in range(aa_range[0],aa_range[1]+1)]
    # Wild type name of the sequence
    wt_name=args.wt_name
    # replicate configure
    replicate_con=pd.read_csv(filepath_or_buffer=args.replicate_con,sep='\t',index_col='Condition')
    condition_Lst=replicate_con.columns.tolist()
    replicate_Lst=replicate_con.index.tolist()
    condition_Dic={condition:replicate_con[replicate_con[condition]==1].index.tolist() for condition in condition_Lst}
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
    # calculated Log10u for sensors in multiple replicates
    Log10u_DF=pd.read_csv(filepath_or_buffer=args.log10u,sep='\t',index_col='Sensor')
    # calculated sigma for sensors in multiple replicates
    sigma_DF=pd.read_csv(filepath_or_buffer=args.sigma,sep='\t',index_col='Sensor')
    # output directory
    list_files=args.output_prefix.split('/')
    output_dir=''
    for i in range(len(list_files)-1):
        output_dir=output_dir+list_files[i]+'/'
    os.system('mkdir -p %s' %(output_dir))
    """
    Calculate the consistency of replicates and means of biological replicates
    """
    # remove sensor mutants filtered by the cell count threshold (with NaN in dataframe)
    Log10u_dropna_DF=Log10u_DF.dropna(axis=0, how='any')
    sigma_dropna_DF=sigma_DF.dropna(axis=0, how='any')
    # draw scatter plot of biological repliates
    for (condition,v) in condition_Dic.items():
        scatter_plot(np.array(Log10u_dropna_DF[v[0]]),np.array(Log10u_dropna_DF[v[1]]),'Response, replicate 1','Response, replicate 2',output_dir,condition+'_Log10u_replicate',6,6)
        scatter_plot(np.array(sigma_dropna_DF[v[0]]),np.array(sigma_dropna_DF[v[1]]),'Response, replicate 1','Response, replicate 2',output_dir,condition+'_sigma_replicate',6,6)
        # calculate the means
        Log10u_DF[condition]=Log10u_DF[v].mean(axis=1)
        sigma_DF[condition]=sigma_DF[v].mean(axis=1)
        # draw the histogram
        hist_plot(np.array(Log10u_DF[condition]),(range_Log10u[0],range_Log10u[1]),50,'Log10u_%s'%(condition),output_dir,condition+'Log10u_hist')
        hist_plot(np.array(sigma_DF[condition]),(range_sigma[0],range_sigma[1]),50,'sigma_%s'%(condition),output_dir,condition+'sigma_hist')
    """
    add codon wt sequence to the df, as new rows
    """
    wtcodon_Dic={}
    for p in [i for i,x in enumerate(wt_seq) if i%3==0]:
        aa_p=p/3+1
        if aa_p>=aa_range[0] and aa_p<=aa_range[1]:
            codon=wt_seq[p:p+3]
            wtcodon_Dic[aa_p]=codon
            wt_aa=standard_table.forward_table[codon]
            name='TnaC_%s%d%s_%s'%(wt_aa,aa_p,wt_aa,codon)
            Log10u_DF.loc[name]=Log10u_DF.loc[wt_name]
            sigma_DF.loc[name]=sigma_DF.loc[wt_name]
    # remove row with 'TnaC' index
    Log10u_DF=Log10u_DF[Log10u_DF.index.str.rstrip()!= wt_name]
    sigma_DF=sigma_DF[sigma_DF.index.str.rstrip()!= wt_name]
    # dict with index as keys for wt codon
    name_Lst=Log10u_DF.index.tolist()
    name_codon_Dic={x.replace('prematureSC','*'):wtcodon_Dic[int(x.split('_')[1].replace('prematureSC','*')[1:-1])] for x in name_Lst}
    """
    extract information from the name of each mutant, including aa, codon, position and whether aa == wt or not
    """
    Log10u_DF['Name']=Log10u_DF.index
    sigma_DF['Name']=sigma_DF.index
    Log10u_DF['Name']=Log10u_DF['Name'].str.replace('prematureSC','*')
    sigma_DF['Name']=sigma_DF['Name'].str.replace('prematureSC','*')
    Log10u_DF['wt_aa']=Log10u_DF['Name'].str.split('_').str[1].str[0]
    sigma_DF['wt_aa']=sigma_DF['Name'].str.split('_').str[1].str[0]
    Log10u_DF['wt_codon']=pd.Series(name_codon_Dic)
    sigma_DF['wt_codon']=pd.Series(name_codon_Dic)
    Log10u_DF['aa']=Log10u_DF['Name'].str.split('_').str[1].str[-1]
    sigma_DF['aa']=sigma_DF['Name'].str.split('_').str[1].str[-1]
    Log10u_DF['codon']=Log10u_DF['Name'].str.split('_').str[2]
    sigma_DF['codon']=sigma_DF['Name'].str.split('_').str[2]
    Log10u_DF['Position']=Log10u_DF['Name'].str.split('_').str[1].str[1:-1].apply(int)
    sigma_DF['Position']=sigma_DF['Name'].str.split('_').str[1].str[1:-1].apply(int)
    Log10u_DF['whetherWTaa']=np.where((Log10u_DF['aa']==Log10u_DF['wt_aa']),'WT','No')
    sigma_DF['whetherWTaa']=np.where((sigma_DF['aa']==sigma_DF['wt_aa']),'WT','No')
    Log10u_DF['whetherWTcodon']=np.where((Log10u_DF['codon']==Log10u_DF['wt_codon']),'WT','No')
    sigma_DF['whetherWTcodon']=np.where((sigma_DF['codon']==sigma_DF['wt_codon']),'WT','No')
    Log10u_DF['100uM-fold']=Log10u_DF['Ligand=100uM']-Log10u_DF['Ligand=0uM']
    Log10u_DF['500uM-fold']=Log10u_DF['Ligand=500uM']-Log10u_DF['Ligand=0uM']
    Log10u_DF['500-100uM-fold']=Log10u_DF['Ligand=500uM']-Log10u_DF['Ligand=100uM']
    Log10u_dropna_DF=Log10u_DF.dropna(axis=0, how='any')
    scatter_plot(np.array(Log10u_dropna_DF['Ligand=0uM']),np.array(Log10u_dropna_DF['500uM-fold']),'0uM Expression','0-500uM fold change',output_dir,'Fold_change_vs.Leakage',6,6)
    Log10u_DF.to_csv('%salldata.csv'%(output_dir),sep='\t')
    wt_matrix=Log10u_DF[(Log10u_DF['whetherWTcodon']=='WT')]
    wt_aa_Dic={'Xaxis':{},'Yaxis':{}}
    wt_codon_Dic={'Xaxis':{},'Yaxis':{}}
    for index, row in wt_matrix.iterrows():
        wt_aa_Dic['Xaxis'][index]=row['Position']-aa_range[0]
        wt_aa_Dic['Yaxis'][index]=aa_order.index(row['aa'])
        wt_codon_Dic['Xaxis'][index]=row['Position']-aa_range[0]
        wt_codon_Dic['Yaxis'][index]=codon_order.index(row['codon']+'_'+row['aa'])
    wt_aa_matrix=pd.DataFrame(wt_aa_Dic)
    wt_codon_matrix=pd.DataFrame(wt_codon_Dic)
    """
    Category the sensor mutant pair based on synonymous and unsynonymous mutations, test the difference between them by MWU test
    """
    compare_syn_un(Log10u_DF,'Ligand=0uM',aa_range_Lst,'Log10u_0uM_syn_vs_unsyn',output_dir)
    compare_syn_un(Log10u_DF,'Ligand=100uM',aa_range_Lst,'Log10u_100uM_syn_vs_unsyn',output_dir)
    compare_syn_un(Log10u_DF,'Ligand=500uM',aa_range_Lst,'Log10u_500uM_syn_vs_unsyn',output_dir)
    """
    convert to Xaxis:position, Yaxis:aa matrix
    """
    Log10u_0uM_aaDF=codon_to_aa(Log10u_DF,'Ligand=0uM',aa_order,aa_range_Lst)
    Log10u_100uM_aaDF=codon_to_aa(Log10u_DF,'Ligand=100uM',aa_order,aa_range_Lst)
    Log10u_500uM_aaDF=codon_to_aa(Log10u_DF,'Ligand=500uM',aa_order,aa_range_Lst)
    Log10u_100uM_foldDF=codon_to_aa(Log10u_DF,'100uM-fold',aa_order,aa_range_Lst)
    Log10u_500uM_foldDF=codon_to_aa(Log10u_DF,'500uM-fold',aa_order,aa_range_Lst)
    Log10u_500_100uM_foldDF=codon_to_aa(Log10u_DF,'500-100uM-fold',aa_order,aa_range_Lst)
    sigma_0uM_aaDF=codon_to_aa(sigma_DF,'Ligand=0uM',aa_order,aa_range_Lst)
    sigma_100uM_aaDF=codon_to_aa(sigma_DF,'Ligand=100uM',aa_order,aa_range_Lst)
    sigma_500uM_aaDF=codon_to_aa(sigma_DF,'Ligand=500uM',aa_order,aa_range_Lst)
    Log10u_0uM_aanoiseDF=codon_to_aa_noise(Log10u_DF,sigma_DF,'Ligand=0uM',aa_order,aa_range_Lst)
    Log10u_100uM_aanoiseDF=codon_to_aa_noise(Log10u_DF,sigma_DF,'Ligand=100uM',aa_order,aa_range_Lst)
    Log10u_500uM_aanoiseDF=codon_to_aa_noise(Log10u_DF,sigma_DF,'Ligand=500uM',aa_order,aa_range_Lst)
    """
    convert to Xaxis:position, Yaxis:codon matrix
    """
    Log10u_0uM_codonDF=codon_to_codon(Log10u_DF,'Ligand=0uM',codon_order_Lst,aa_range_Lst)
    Log10u_100uM_codonDF=codon_to_codon(Log10u_DF,'Ligand=100uM',codon_order_Lst,aa_range_Lst)
    Log10u_500uM_codonDF=codon_to_codon(Log10u_DF,'Ligand=500uM',codon_order_Lst,aa_range_Lst)
    Log10u_100uM_fold_codonDF=codon_to_codon(Log10u_DF,'100uM-fold',codon_order_Lst,aa_range_Lst)
    Log10u_500uM_fold_codonDF=codon_to_codon(Log10u_DF,'500uM-fold',codon_order_Lst,aa_range_Lst)
    Log10u_500_100uM_fold_codonDF=codon_to_codon(Log10u_DF,'500-100uM-fold',codon_order_Lst,aa_range_Lst)
    sigma_0uM_codonDF=codon_to_codon(sigma_DF,'Ligand=0uM',codon_order_Lst,aa_range_Lst)
    sigma_100uM_codonDF=codon_to_codon(sigma_DF,'Ligand=100uM',codon_order_Lst,aa_range_Lst)
    sigma_500uM_codonDF=codon_to_codon(sigma_DF,'Ligand=500uM',codon_order_Lst,aa_range_Lst)
    """
    draw the heatmap
    """
    vmin=range_Log10u[0]
    vmax=range_Log10u[1]
    temp=max(abs(vmin),abs(vmax))
    vmin=-temp
    vmax=temp
    aspect=20
    shrink=0.8
    tick_size=16
    Xlen=8
    Ylen=6
    heatmap_plot(np.array(Log10u_0uM_aaDF),plt.cm.seismic,vmin,vmax,'Position','Amino acid',aa_range_Lst,aa_order,tick_size,output_dir,'Log10u_0uM_matrix',wt_aa_matrix,Xlen,Ylen,'Response',aspect,shrink)
    heatmap_plot(np.array(Log10u_100uM_aaDF),plt.cm.seismic,vmin,vmax,'Position','Amino acid',aa_range_Lst,aa_order,tick_size,output_dir,'Log10u_100uM_matrix',wt_aa_matrix,Xlen,Ylen,'Response',aspect,shrink)
    heatmap_plot(np.array(Log10u_500uM_aaDF),plt.cm.seismic,vmin,vmax,'Position','Amino acid',aa_range_Lst,aa_order,tick_size,output_dir,'Log10u_500uM_matrix',wt_aa_matrix,Xlen,Ylen,'Response',aspect,shrink)
    aspect=20
    shrink=0.6
    tick_size=8
    Xlen=5
    Ylen=8
    heatmap_plot(np.array(Log10u_0uM_codonDF),plt.cm.seismic,vmin,vmax,'Position','Codon',aa_range_Lst,codon_order,tick_size,output_dir,'Log10u_0uM_codon_matrix',wt_codon_matrix,Xlen,Ylen,'Response',aspect,shrink)
    heatmap_plot(np.array(Log10u_100uM_codonDF),plt.cm.seismic,vmin,vmax,'Position','Codon',aa_range_Lst,codon_order,tick_size,output_dir,'Log10u_100uM_codon_matrix',wt_codon_matrix,Xlen,Ylen,'Response',aspect,shrink)
    heatmap_plot(np.array(Log10u_500uM_codonDF),plt.cm.seismic,vmin,vmax,'Position','Codon',aa_range_Lst,codon_order,tick_size,output_dir,'Log10u_500uM_codon_matrix',wt_codon_matrix,Xlen,Ylen,'Response',aspect,shrink)
    vmin=range_sigma[0]
    vmax=range_sigma[1]
    aspect=20
    shrink=0.8
    tick_size=16
    Xlen=8
    Ylen=6
    heatmap_plot(np.array(sigma_0uM_aaDF),plt.cm.Greys,vmin,vmax,'Position','Amino acid',aa_range_Lst,aa_order,tick_size,output_dir,'sigma_0uM_matrix',wt_aa_matrix,Xlen,Ylen,'Sigma of response',aspect,shrink)
    heatmap_plot(np.array(sigma_100uM_aaDF),plt.cm.Greys,vmin,vmax,'Position','Amino acid',aa_range_Lst,aa_order,tick_size,output_dir,'sigma_100uM_matrix',wt_aa_matrix,Xlen,Ylen,'Sigma of response',aspect,shrink)
    heatmap_plot(np.array(sigma_500uM_aaDF),plt.cm.Greys,vmin,vmax,'Position','Amino acid',aa_range_Lst,aa_order,tick_size,output_dir,'sigma_500uM_matrix',wt_aa_matrix,Xlen,Ylen,'Sigma of response',aspect,shrink)
    aspect=20
    shrink=0.6
    tick_size=8
    Xlen=5
    Ylen=8
    heatmap_plot(np.array(sigma_0uM_codonDF),plt.cm.Greys,vmin,vmax,'Position','Codon',aa_range_Lst,codon_order,tick_size,output_dir,'sigma_0uM_codon_matrix',wt_codon_matrix,Xlen,Ylen,'Sigma of response',aspect,shrink)
    heatmap_plot(np.array(sigma_100uM_codonDF),plt.cm.Greys,vmin,vmax,'Position','Codon',aa_range_Lst,codon_order,tick_size,output_dir,'sigma_100uM_codon_matrix',wt_codon_matrix,Xlen,Ylen,'Sigma of response',aspect,shrink)
    heatmap_plot(np.array(sigma_500uM_codonDF),plt.cm.Greys,vmin,vmax,'Position','Codon',aa_range_Lst,codon_order,tick_size,output_dir,'sigma_500uM_codon_matrix',wt_codon_matrix,Xlen,Ylen,'RSigma of response',aspect,shrink)
    temp=np.append(np.array(Log10u_100uM_foldDF),np.array(Log10u_500uM_foldDF))
    vmin=np.nanmin(temp)
    vmax=np.nanmax(temp)
    temp=max(abs(vmin),abs(vmax))
    vmin=-temp
    vmax=temp
    aspect=20
    shrink=0.8
    tick_size=16
    Xlen=8
    Ylen=6
    heatmap_plot(np.array(Log10u_100uM_foldDF),plt.cm.seismic,vmin,vmax,'Position','Amino acid',aa_range_Lst,aa_order,tick_size,output_dir,'Log10u_100uM_fold_matrix',wt_aa_matrix,Xlen,Ylen,'Induction (100 uM vs. 0 uM)',aspect,shrink)
    heatmap_plot(np.array(Log10u_500uM_foldDF),plt.cm.seismic,vmin,vmax,'Position','Amino acid',aa_range_Lst,aa_order,tick_size,output_dir,'Log10u_500uM_fold_matrix',wt_aa_matrix,Xlen,Ylen,'Induction (500 uM vs. 0 uM)',aspect,shrink)
    heatmap_plot(np.array(Log10u_500_100uM_foldDF),plt.cm.seismic,vmin,vmax,'Position','Amino acid',aa_range_Lst,aa_order,tick_size,output_dir,'Log10u_500vs100uM_fold_matrix',wt_aa_matrix,Xlen,Ylen,'Induction (500 uM vs. 100 uM)',aspect,shrink)
    aspect=20
    shrink=0.6
    tick_size=8
    Xlen=5
    Ylen=8
    heatmap_plot(np.array(Log10u_100uM_fold_codonDF),plt.cm.seismic,vmin,vmax,'Position','Codon',aa_range_Lst,codon_order,tick_size,output_dir,'Log10u_100uM_fold_codon_matrix',wt_codon_matrix,Xlen,Ylen,'Induction (100 uM vs. 0 uM)',aspect,shrink)
    heatmap_plot(np.array(Log10u_500uM_fold_codonDF),plt.cm.seismic,vmin,vmax,'Position','Codon',aa_range_Lst,codon_order,tick_size,output_dir,'Log10u_500uM_fold_codon_matrix',wt_codon_matrix,Xlen,Ylen,'Induction (500 uM vs. 0 uM)',aspect,shrink)
    heatmap_plot(np.array(Log10u_500_100uM_fold_codonDF),plt.cm.seismic,vmin,vmax,'Position','Codon',aa_range_Lst,codon_order,tick_size,output_dir,'Log10u_500vs100uM_fold_codon_matrix',wt_codon_matrix,Xlen,Ylen,'Induction (500 uM vs. 100 uM)',aspect,shrink)
    #temp=np.append(np.array(Log10u_0uM_aanoiseDF),np.array(Log10u_100uM_aanoiseDF),np.array(Log10u_500uM_aanoiseDF))
    vmin=0.0
    vmax=30.0
    aspect=20
    shrink=0.8
    tick_size=16
    Xlen=8
    Ylen=6
    heatmap_plot(np.array(Log10u_0uM_aanoiseDF),plt.cm.Greens,vmin,vmax,'Position','Amino acid',aa_range_Lst,aa_order,tick_size,output_dir,'Log10u_0uM_dege_noise_matrix',wt_aa_matrix,Xlen,Ylen,'Response deviation between\nsynonymous codons, -log$_{10}P$',aspect,shrink)
    heatmap_plot(np.array(Log10u_100uM_aanoiseDF),plt.cm.Greens,vmin,vmax,'Position','Amino acid',aa_range_Lst,aa_order,tick_size,output_dir,'Log10u_100uM_dege_noise_matrix',wt_aa_matrix,Xlen,Ylen,'Response deviation between\nsynonymous codons, -log$_{10}P$)',aspect,shrink)
    heatmap_plot(np.array(Log10u_500uM_aanoiseDF),plt.cm.Greens,vmin,vmax,'Position','Amino acid',aa_range_Lst,aa_order,tick_size,output_dir,'Log10u_500uM_dege_noise_matrix',wt_aa_matrix,Xlen,Ylen,'Response deviation between\nsynonymous codons, -log$_{10}P$',aspect,shrink)


if __name__ == '__main__':
    try:
        args=htsensor_plot_parseargs()
        htsensor_plot_main(args)
    except KeyboardInterrupt:
        sys.stderr.write("Interrupted.\n")
        sys.exit(0)

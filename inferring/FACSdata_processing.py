# this script is used to process the raw .fcs data for FACS analysis

import os
import sys
import FlowCytometryTools
from FlowCytometryTools import FCMeasurement
from FlowCytometryTools import ThresholdGate, PolyGate
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

from matplotlib import cm
import warnings

dataFile=sys.argv[1]
print dataFile
prefix=sys.argv[2]
sample = FCMeasurement(ID=prefix, datafile=dataFile)
print sample.data.describe()
print ''

# Log10 treatment for selected channels of FACS raw data
def custom_log(original_sample,channelLst):
    # Copy the original sample
    new_sample = original_sample.copy()
    new_data = new_sample.data
    original_data = original_sample.data
    
    # Our transformation goes here
    for channel in channelLst:
        new_data[channel] = np.log10(original_data[channel])
    new_data = new_data.dropna() # Removes all NaN entries
    new_sample.data = new_data
    return new_sample

# //////////////////////////////////////////////
# threshold gating to remove all negative FSC and SSC events
preP1th_value1=10.0**1.6699667767116002
preP1th_channel1='FSC-A'
preP1th_value2=10.0**2.8012542209507507
preP1th_channel2='SSC-A'
preP1gate1=ThresholdGate(preP1th_value1, preP1th_channel1, region='above')
preP1gate2=ThresholdGate(preP1th_value2, preP1th_channel2, region='above')
preP1gate=preP1gate1 & preP1gate2

preP1_sample = sample.gate(preP1gate)
print 'profile after preP1 gating'
print preP1_sample.data.describe()
print ''

preP1log = custom_log (preP1_sample, ['FSC-A', 'SSC-A'])
print 'profile after preP1 gating and log10 for FSC and SSC'
print preP1log.data.describe()
print ''

# //////////////////////////////////////////////
# P1 gating based on FSC and SSC (Log10 space)
P1_gatevalue=[(4.481395264353071,4.727180089479611),(1.761395345415388,3.4679208677492026),(1.6699667767116002,2.8259455782395824),(3.7042524303708757,2.8012542209507507),(4.869966681344168,4.529649231168959)]
P1gate_channel=('FSC-A','SSC-A')
P1gate = PolyGate(P1_gatevalue, P1gate_channel, region='in', name='P1 gate')
P1_sample=preP1log.gate(P1gate)
print 'profile after P1 gating'
print P1_sample.data.describe()
print ''

# ////////////////////////////////////
# threshold gate to remove all negative mCherry (PE-Texas Red) and GFP (FITC) events
prePxth_value1=10.0**2.0065723832855875
prePxth_channel1='PE-Texas Red-A'
prePxth_value2=10.0**1.426888943216241
prePxth_channel2='FITC-A'
prePxgate1=ThresholdGate(prePxth_value1, prePxth_channel1, region='above')
prePxgate2=ThresholdGate(prePxth_value2, prePxth_channel2, region='above')
prePxgate=prePxgate1 & prePxgate2

prePx_sample = P1_sample.gate(prePxgate)
print 'profile after prePx gating'
print prePx_sample.data.describe()
print ''

prePxlog = custom_log (prePx_sample, ['PE-Texas Red-A', 'FITC-A'])
prePxlog.data['Log10u'] = prePxlog.data['FITC-A']-prePxlog.data['PE-Texas Red-A']
print 'profile after prePx gating and log10 for PE-Texas Red and FITC'
print prePxlog.data.describe()
print ''

# //////////////////////////////////////////////
# Px gating based on mCherry and GFP (Log10 space)
P2_gatevalue=[(4.898879989395794,5.426888824006959),(2.0065723832855875,2.854863848467208),(2.0065723832855875,5.418538093566893)]
P2gate_channel=('PE-Texas Red-A','FITC-A')
P2gate = PolyGate(P2_gatevalue, P2gate_channel, region='in', name='P2 gate')

P3_gatevalue=[(2.0065723832855875,2.8465131180271332),(2.0065723832855875,2.462379517784184),(5.04930733534006,5.142963989044773),(5.056144941973891,5.4101873631268305),(4.905717596029625,5.426888824006957)]
P3gate_channel=('PE-Texas Red-A','FITC-A')
P3gate = PolyGate(P3_gatevalue, P3gate_channel, region='in', name='P3 gate')

P4_gatevalue=[(2.0065723832855875,2.4540287873441176),(2.0065723832855875,2.0698951871011677),(5.04930733534006,4.75047965836176),(5.056144941973891,5.017703032443813),(5.04930733534006,5.134613258604714)]
P4gate_channel=('PE-Texas Red-A','FITC-A')
P4gate = PolyGate(P4_gatevalue, P4gate_channel, region='in', name='P4 gate')

P5_gatevalue=[(2.0065723832855875,2.061544456661103),(2.0065723832855875,1.6774108564181534),(5.04930733534006,4.357995327678745),(5.056144941973891,4.625218701760799),(5.04930733534006,4.7421289279217005)]
P5gate_channel=('PE-Texas Red-A','FITC-A')
P5gate = PolyGate(P5_gatevalue, P5gate_channel, region='in', name='P5 gate')

P6_gatevalue=[(2.0065723832855875,1.6774108564181365),(2.0065723832855875,1.435239673656306),(2.184350155765174,1.4352396736563082),(5.04930733534006,3.982212457875856),(5.04930733534006,4.357995327678738)]
P6gate_channel=('PE-Texas Red-A','FITC-A')
P6gate = PolyGate(P6_gatevalue, P6gate_channel, region='in', name='P6 gate')

P7_gatevalue=[(5.04930733534006,3.9738617274358004),(2.184350155765174,1.426888943216241),(5.042469728706229,1.4268889432162477)]
P7gate_channel=('PE-Texas Red-A','FITC-A')
P7gate = PolyGate(P7_gatevalue, P7gate_channel, region='in', name='P7 gate')

# 'Or' Px gate to construct an overall fluorescnece gate
P234567_gate=P2gate | P3gate | P4gate | P5gate | P6gate | P7gate
P234567_sample = prePxlog.gate(P234567_gate)
print 'profile of P1234567 population'
print P234567_sample.data.describe()
print ''
if prefix=='wt':
    meanZLog10u, stdZLog10u = norm.fit(P234567_sample.data['Log10u'])
    meanZLog10GFP, stdZLog10GFP = norm.fit(P234567_sample.data['FITC-A'])
    print 'mean = %s, std = %s for fitting ND of Log10u'%(str(meanZLog10u), str(stdZLog10u))
    print 'mean = %s, std = %s for fitting ND of Log10GFP'%(str(meanZLog10GFP), str(stdZLog10GFP))
    print ''
    os.system('cat /dev/null > %s_NDpara.csv'%(prefix))
    g=open('%s_NDpara.csv'%(prefix),'r+')
    g.write('Log10u_mean\tLog10u_std\n')
    g.write('%s\t%s\n'%(str(meanZLog10u), str(stdZLog10u)))
    g.write('Log10GFP_mean\tLog10GFP_std\n')
    g.write('%s\t%s\n'%(str(meanZLog10GFP), str(stdZLog10GFP)))
    g.close()
Log10uMax=P234567_sample.data.describe()['Log10u']['max']
Log10uMin=P234567_sample.data.describe()['Log10u']['min']
Log10uGFP=P234567_sample.data.describe()['FITC-A']['max']
Log10uGFP=P234567_sample.data.describe()['FITC-A']['min']
#bins=np.linspace(Log10uMin, Log10uMax, 200)
bins=np.linspace(-2.0, 5.0, 500)
plt.hist(np.array(P234567_sample.data['Log10u']),bins,alpha=0.5, color='r', linewidth=0.0, label='Log10 (GFP/mCherry)')
plt.hist(np.array(P234567_sample.data['FITC-A']),bins,alpha=0.5, color='g', linewidth=0.0, label='Log10 GFP')
plt.ylabel('cytometry event')
plt.xlim([-2.0,5.0])
plt.legend(loc='upper left')
plt.savefig('%s.all.png'%(prefix))
plt.clf()

P1P2_sample = prePxlog.gate(P2gate)
P1P3_sample = prePxlog.gate(P3gate)
P1P4_sample = prePxlog.gate(P4gate)
P1P5_sample = prePxlog.gate(P5gate)
P1P6_sample = prePxlog.gate(P6gate)
P1P7_sample = prePxlog.gate(P7gate)

P12Logu=np.array(P1P2_sample.data['Log10u'])
P12Logu.sort()
P13Logu=np.array(P1P3_sample.data['Log10u'])
P13Logu.sort()
P14Logu=np.array(P1P4_sample.data['Log10u'])
P14Logu.sort()
P15Logu=np.array(P1P5_sample.data['Log10u'])
P15Logu.sort()
P16Logu=np.array(P1P6_sample.data['Log10u'])
P16Logu.sort()
P17Logu=np.array(P1P7_sample.data['Log10u'])
P17Logu.sort()

bins=np.linspace(-2.0, 2.0, 200)
plt.hist(P12Logu,bins,color='b',label='P1&P2',alpha=0.5,linewidth=0.0)
plt.hist(P13Logu,bins,color='r',label='P1&P3',alpha=0.5,linewidth=0.0)
plt.hist(P14Logu,bins,color='g',label='P1&P4',alpha=0.5,linewidth=0.0)
plt.hist(P15Logu,bins,color='y',label='P1&P5',alpha=0.5,linewidth=0.0)
plt.hist(P16Logu,bins,color='c',label='P1&P6',alpha=0.5,linewidth=0.0)
plt.hist(P17Logu,bins,color='m',label='P1&P7',alpha=0.5,linewidth=0.0)
plt.legend(loc='upper left')
plt.xlabel('Log10 (GFP/mCherry)')
plt.ylabel('cytometry event')
plt.xlim([-2.0,1.5])
plt.savefig('%s.gate.png'%(prefix))
plt.clf()

print 'profile of P1&P2 population'
print P1P2_sample.data.describe()
print ''
print 'profile of P1&P3 population'
print P1P3_sample.data.describe()
print ''
print 'profile of P1&P4 population'
print P1P4_sample.data.describe()
print ''
print 'profile of P1&P5 population'
print P1P5_sample.data.describe()
print ''
print 'profile of P1&P6 population'
print P1P6_sample.data.describe()
print ''
print 'profile of P1&P7 population'
print P1P7_sample.data.describe()
print ''

P234567_sample.data.describe().to_csv('P1P234567.stats.csv')
P1P2_sample.data.describe().to_csv('P1P2.stats.csv')
P1P3_sample.data.describe().to_csv('P1P3.stats.csv')
P1P4_sample.data.describe().to_csv('P1P4.stats.csv')
P1P5_sample.data.describe().to_csv('P1P5.stats.csv')
P1P6_sample.data.describe().to_csv('P1P6.stats.csv')
P1P7_sample.data.describe().to_csv('P1P7.stats.csv')
P1P2ratio=float(P1P2_sample.data.describe().loc['count','FITC-A'])/float(P234567_sample.data.describe().loc['count','FITC-A'])
P1P3ratio=float(P1P3_sample.data.describe().loc['count','FITC-A'])/float(P234567_sample.data.describe().loc['count','FITC-A'])
P1P4ratio=float(P1P4_sample.data.describe().loc['count','FITC-A'])/float(P234567_sample.data.describe().loc['count','FITC-A'])
P1P5ratio=float(P1P5_sample.data.describe().loc['count','FITC-A'])/float(P234567_sample.data.describe().loc['count','FITC-A'])
P1P6ratio=float(P1P6_sample.data.describe().loc['count','FITC-A'])/float(P234567_sample.data.describe().loc['count','FITC-A'])
P1P7ratio=float(P1P7_sample.data.describe().loc['count','FITC-A'])/float(P234567_sample.data.describe().loc['count','FITC-A'])
os.system('cat /dev/null > %s_ratio.csv'%(prefix))
g=open('%s_ratio.csv'%(prefix),'r+')
g.write('%s bin occupation against P1P234567\n'%(prefix))
g.write('P1P2\tP1P3\tP1P4\tP1P5\tP1P6\tP1P7\n')
g.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(str(P1P2ratio), str(P1P3ratio), str(P1P4ratio), str(P1P5ratio), str(P1P6ratio), str(P1P7ratio)))
g.close()

# /////////////////////////////////////////////////////////////
# get the overlap between each two neighbour bins on the Log10u space, using which to calculate the Log10u boundary of two neighbour bins
# The function to define the overlap between each two neighbour bins on the Log10u space
# Briefly, feed two arrays (Log10u of two bins) and a grid bin space; search the grid space one bin by another, extract subarray within this small bin from two arrays and compare with each other. When both arrays have members within this bin and the comtamination (ratio of member belonging to two arrays) > 1%, extract the minority append it into a new array. calculate the mean of the new array.
def find_ovarlap(array1,array2,bins):
    array1copy=array1
    array2copy=array2
    array1copy.sort()
    array2copy.sort()
    bins.sort()
    myarray=np.array([])
    for i,item in enumerate(bins):
        if i<len(bins)-1:
            thisbin=[bins[i],bins[i+1]]
            array1here=array1copy[(array1copy>thisbin[0]) & (array1copy<thisbin[1])]
            array2here=array2copy[(array2copy>thisbin[0]) & (array2copy<thisbin[1])]
            if len(array1here)>0 and len(array2here)>0:
                if len(array1here)<len(array2here) and float(len(array1here))/float(len(array2here))>0.01:
                    myarray=np.append(array1here,myarray)
                elif len(array2here)<len(array1here) and float(len(array2here))/float(len(array1here))>0.01:
                    myarray=np.append(array2here,myarray)
    return myarray

P2P3overlap=find_ovarlap(P12Logu, P13Logu, bins)
P3P4overlap=find_ovarlap(P13Logu, P14Logu, bins)
P4P5overlap=find_ovarlap(P15Logu, P14Logu, bins)
P5P6overlap=find_ovarlap(P15Logu, P16Logu, bins)
P6P7overlap=find_ovarlap(P16Logu, P17Logu, bins)

os.system('cat /dev/null > %s_boundary.csv'%(prefix))
g=open('%s_boundary.csv'%(prefix),'r+')
g.write('%s overlap boundary\n'%(prefix))
g.write('P2P3\tP3P4\tP4P5\tP5P6\tP6P7\n')
g.write('%s\t%s\t%s\t%s\t%s\n'%(str(np.mean(P2P3overlap)),str(np.mean(P3P4overlap)),str(np.mean(P4P5overlap)),str(np.mean(P5P6overlap)),str(np.mean(P6P7overlap)),))
g.write('Overlap event number, threshold = 0.01\n')
g.write('P2P3\tP3P4\tP4P5\tP5P6\tP6P7\n')
g.write('%d\t%d\t%d\t%d\t%d\n'%(len(P2P3overlap),len(P3P4overlap),len(P4P5overlap),len(P5P6overlap),len(P6P7overlap)))
g.close()
os.system('mkdir %s'%(prefix))
os.system('mv *.csv *.png %s'%(prefix))

print 'P2P3 overlap: '
print len(P2P3overlap)
print 'P2P3 mean'
print np.mean(P2P3overlap)
print norm.fit(P2P3overlap)
print ''
print 'P3P4 overlap: '
print len(P3P4overlap)
print 'P3P4 mean'
print np.mean(P3P4overlap)
print norm.fit(P3P4overlap)
print ''
print 'P4P5 overlap: '
print len(P4P5overlap)
print 'P4P5 mean'
print np.mean(P4P5overlap)
print norm.fit(P4P5overlap)
print ''
print 'P5P6 overlap: '
print len(P5P6overlap)
print 'P5P6 mean'
print np.mean(P5P6overlap)
print norm.fit(P5P6overlap)
print ''
print 'P6P7 overlap: '
print len(P6P7overlap)
print 'P6P7 mean'
print np.mean(P6P7overlap)
print norm.fit(P6P7overlap)
print ''

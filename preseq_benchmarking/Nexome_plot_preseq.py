import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
from scipy import interpolate
import os 
#use .matplotlib-1.3.1-python-2.7.1-sqlite3-rtrees
#use .scipy-0.13.0-python-2.7.1-sqlite3-rtrees

sample_root = '/seq/picard_aggregation/D5227'
out_base = '/humgen/gsa-hpprojects/dev/hogstrom/depth_by_read_group/Nexome'
# picard_path = '/seq/software/picard/current/bin/picard.jar'
picard_path='/home/unix/hogstrom/picard_dup_by_insert.jar' #modified jar

# Nexome samples
samples=['359781',
 '361433',
 '359877',
 '360361',
 '360457',
 '362428',
 '361337',
 '379543',
 '388072',
 '379639',
 '381153',
 '384114',
 '364464',
 '375034',
 '377582',
 '388264',
 '384210',
 '379214',
 '384498',
 '380335',
 '372754',
 '387332',
 '386674',
 '391478',
 '393370',
 '363907',
 '385972',
 '388168',
 '373158',
 '387976',
 '379118',
 '376110',
 '385390',
 '384594',
 '391382',
 '371834',
 '383609',
 '381057',
 '386068',
 '384694',
 '383320',
 '384790',
 '383416',
 '387428',
 '382684',
 '382191',
 '392292',
 '376734']

############################
### plot preseq results ####
############################

wkdir = out_base + '/preseq_curve_estimates' 
if (not os.path.isdir(wkdir)):
    os.mkdir(wkdir)
for sample in samples:
    print sample
    sample = 'NexPond-' + sample 
    outdir = out_base + '/' + sample
    # ccurve  = wkdir +'/' + sample+'_rgSet' + str(16) + '.c_curve'
    ccurve  = wkdir +'/' + sample+'_chrm21' + '.c_curve'
    iextrap = 3
    lcextrap = wkdir +'/' + sample+'_rgSet' + str(iextrap) + '.lc_extrap'
    # outlog  = wkdir +'/' + sample+'_rgSet' + str(i) + '_c_curve.log'
    ye = pd.read_csv(lcextrap,sep='\t',index_col=0) #load in yield estimates
    cc = pd.read_csv(ccurve,sep='\t',index_col=0) #load in yield estimates
    ### load picard dup counts
    fin2 = outdir + '/' + sample + '_chrm21_read_grp_dup_metrics.txt'
    dm = pd.read_csv(fin2,index_col=0) #load in yield estimates
    #dm.PERCENT_DUPLICATION.convert_objects('convert_numeric')
    dm['UNIQUE_READS'] = dm.READ_PAIRS_EXAMINED - dm.READ_PAIR_DUPLICATES
    ### make plot
    ye.plot(ye.index)
    plt.plot(cc.index.values,cc.distinct_reads.values,'or',label='preseq observed')
    plt.plot(dm.READ_PAIRS_EXAMINED.values,dm.UNIQUE_READS.values,'og',label='picard observed')
    plt.legend(loc='lower right',numpoints=1)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.title('preseq yield estimate: based on ' + str(iextrap) + ' read \n groups for Nexome sample 376014')
    plt.ylabel('expected distinct')
    plt.xlabel('reads observed')
    outF = wkdir  +'/' + sample+'_rgSet' + str(iextrap) + '_c_curve.png'
    plt.savefig(outF, bbox_inches='tight')
    plt.close()

### plot of read counts
for sample in samples:
    print sample
    sample = 'NexPond-' + sample 
    outdir = out_base + '/' + sample
    # ccurve  = wkdir +'/' + sample+'_rgSet' + str(16) + '.c_curve'
    outlog  = wkdir +'/' + sample+'_chrm21' + '_c_curve.log'    
    df = pd.read_csv(outlog,sep='\t',index_col=0,skiprows=8,names=['read occurences','count']) #load in yield estimates
    isStr = [not ("sample size" in x) for x in df.index.values]
    df = df.ix[isStr]
    df_normed = df['count']/df['count'].sum()
    iInts = [int(x) for x in df_normed.index]
    plt.plot(iInts,df_normed)
#plt.xlim((0, 10))
plt.title('occurences of unique reads in Nexome samples')
plt.ylabel('normed count')
plt.xlabel('read occurences')
outF = wkdir  +'/Nexome_read_counts.png'
plt.savefig(outF, bbox_inches='tight')
plt.close()
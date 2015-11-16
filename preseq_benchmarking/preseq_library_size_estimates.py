import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv

### parse log files, retrieve duplicate metrics, save table to disk ###

# sample NexPond-376014
outdir='/humgen/gsa-hpprojects/dev/hogstrom/depth_by_read_group/376014'

### read in duplication metric details
metric_dict = {}
for i in range(1,10):
    rg = 'rg'+str(i)
    fname = outdir + "/NexPond-376014_chrm21_" + rg + "_dup_Metrics.txt"
    lines=list(csv.reader(open(fname))) 
    nMetrics = lines[6][0].split('\t')
    mvals = lines[7][0].split('\t')
    metric_dict[rg] = mvals

# put into dataframe
df_dup = pd.DataFrame.from_dict(metric_dict,orient='index')
df_dup.columns = nMetrics
# df_dup['platform'] = 'WGS'

#load in picard yield estimates
fname = outdir + '/NexPond-376014_chrm21_rg2_yield_estimates.txt'
g = pd.read_csv(fname,sep='\t')
#lines=list(csv.reader(open(fname))) 

fout = outdir + '/NexPond-376014_chrm21_read_grp_dup_metrics.txt'
df_dup.to_csv(fout)

#scp hogstrom@gsa3:/humgen/gsa-hpprojects/dev/hogstrom/depth_by_read_group/376014/NexPond-376014_chrm21_rg2_yield_estimates.txt "/Users/hogstrom/Dropbox (MIT)/genomic_library_complexity/NexPond-376014_chrm21_rg2_yield_estimates.txt"
#scp hogstrom@gsa3:/humgen/gsa-hpprojects/dev/hogstrom/depth_by_read_group/376014/NexPond-376014_chrm21_read_grp_dup_metrics.txt "/Users/hogstrom/Dropbox (MIT)/genomic_library_complexity/NexPond-376014_chrm21_read_grp_dup_metrics.txt"
#scp hogstrom@gsa3:/humgen/gsa-hpprojects/dev/hogstrom/depth_by_read_group/376014/NexPond-376014_chrm21_rg3_yield_estimates_withPairedEndFlag.txt "/Users/hogstrom/Dropbox (MIT)/genomic_library_complexity/NexPond-376014_chrm21_rg3_yield_estimates_withPairedEndFlag.txt"

### read in locally
import matplotlib.pyplot as plt
import pandas as pd
import os 

wkdir='/Users/hogstrom/Documents/code/palantir/Analysis/library_complexity/figures'
indir='/Users/hogstrom/Dropbox (MIT)/genomic_library_complexity/'
### load preseq yield etimates
#fin1=os.path.join(indir,"NexPond-376014_chrm21_rg3_yield_estimates.txt")
fin1=os.path.join(indir,"NexPond-376014_chrm21_rg3_yield_estimates_withoutPairedEndFlag.txt")
ye = pd.read_csv(fin1,sep='\t',index_col=0) #load in yield estimates

yeshort = ye[:100]
yeshort.plot(yeshort.index)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.title('preseq yield estimate: based on 2 read \n groups for Nexome sample 376014')
plt.ylabel('expected distinct')
plt.xlabel('reads observed')
#plt.show()
outF = os.path.join(wkdir,'376014_chrm21_yield_estimate_curve.png')
plt.savefig(outF, bbox_inches='tight')

### load duplication metrics for observed read groups
fin2="/Users/hogstrom/Dropbox (MIT)/genomic_library_complexity/NexPond-376014_chrm21_read_grp_dup_metrics.txt"
dm = pd.read_csv(fin2,index_col=0) #load in yield estimates
dm['UNIQUE_READS'] = dm.READ_PAIRS_EXAMINED - dm.READ_PAIR_DUPLICATES
dm.plot('READ_PAIRS_EXAMINED','UNIQUE_READS')

### plot observed with predicted 
# yeshort = ye[:15]
yeshort.plot(yeshort.index)
plt.plot(dm.READ_PAIRS_EXAMINED,dm.UNIQUE_READS,'o',label='observed')
plt.legend(loc='lower right')
plt.title('library complexity: based on 3 read \n groups for Nexome sample 376014')
plt.xlabel('read pairs examined')
plt.ylabel('unique reads')
#plt.show()
outF = os.path.join(wkdir,'376014_chrm21_observed_yield_curve_rg3_withoutPairedEndFlag.png')
plt.savefig(outF, bbox_inches='tight')
plt.close()

# check percent duplicates 
# 1 - np.divide(dm.UNIQUE_READS,dm.READ_PAIRS_EXAMINED.astype('float'))
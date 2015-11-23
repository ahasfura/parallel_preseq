import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
from scipy import interpolate
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

# library size estimate equations 
def f(x,c,n):
    return c/float(x) - 1 + np.exp(-n/float(x))

def est_lib_size(n, c):
    # n = number of observed reads
    # c = number of unique reads
    if (c > n):
        print("error in unique vs. total reads")
    m=1
    M=100
    while (f(M*c, c, n) >= 0):
        M *= 10.0
    for i in range(1,41):
        r = (m+M)/float(2)
        u = f(r*c,c,n)
        if (u == 0):
            break
        else:
            if (u>0):
                m = r
            elif (u<0):
                M = r
    return (c * (m+M)/2.0)

def estimateRoi(EstLibSize, x, n, c):
    # EstLibSize = estimated library size
    # x = multiple of sequencing to be simulated (how many X sequencing)
    # n = read pairs observed
    # c = unique read pairs 
    return EstLibSize * (1 - np.exp(-(x*n)/EstLibSize))/c


### plot ROI curves
### calculate extrapolated estimate errors across many samples & multiple models
wkdir= out_base + '/roi_estimates_minus_OptDups'
colSet = ['ROI_percent_error','poly_perc_err_deg1','poly_perc_err_deg2','spline_perc_err_deg1','spline_perc_err_deg2']
dfError = pd.DataFrame(columns=colSet)
for sample in samples:
    sample = 'NexPond-' + sample 
    print sample 
    outdir = out_base + '/' + sample
    fin2 = outdir + '/' + sample + '_chrm21_read_grp_dup_metrics.txt'
    dm = pd.read_csv(fin2,index_col=0) #load in yield estimates
    #dm.PERCENT_DUPLICATION.convert_objects('convert_numeric')
    dm['UNIQUE_READS'] = dm.READ_PAIRS_EXAMINED - dm.READ_PAIR_DUPLICATES
    dm['read_pairs_minus_OptDups'] = dm.READ_PAIRS_EXAMINED - dm.READ_PAIR_OPTICAL_DUPLICATES # adapt read pair count to match picard 
    dm = dm.applymap(lambda x: np.nan if isinstance(x,basestring) else x)
    dm['est_lib_size'] = dm.apply(lambda row: est_lib_size(row['read_pairs_minus_OptDups'], row['UNIQUE_READS']), axis=1)
    #dm[['ESTIMATED_LIBRARY_SIZE','est_lib_size']] #check calculation to markdup output
    ### picard-based ROI estimate
    iRg = 'rgSet3'
    # read_pairs_for_estimation = 'READ_PAIRS_EXAMINED' # same as appears in picard
    read_pairs_for_estimation = 'read_pairs_minus_OptDups' # proposed set for roi prediction
    dm['inMultiple'] = dm[read_pairs_for_estimation]/dm.ix[iRg,read_pairs_for_estimation]
    dm['ROI_' + iRg + '_mult'] = dm.apply(lambda row: estimateRoi(dm.ix[iRg,'est_lib_size'],row['inMultiple'],dm.ix[iRg,read_pairs_for_estimation], dm.ix[iRg,'UNIQUE_READS']), axis=1)
    dm['ROI_' + iRg + '_est'] = dm['ROI_' + iRg + '_mult']*dm.ix[iRg,'UNIQUE_READS']
    rpExamined = float(dm.ix[-1,read_pairs_for_estimation])
    urObserved = float(dm.ix[-1,'UNIQUE_READS'])
    dfError.ix[sample,'ROI_percent_error'] = abs((dm.ix[-1,'ROI_' + iRg + '_est'] - urObserved) / urObserved)
    ### test multiple degree estimates for polynomail & spline
    for sdegree in range(1,3):
        print sdegree
        #polynomial fit
        z = np.polyfit(dm.ix[:iRg,read_pairs_for_estimation].values, dm.ix[:iRg,'UNIQUE_READS'].values, sdegree)
        p = np.poly1d(z)
        dfError.ix[sample,'poly_perc_err_deg' +str(sdegree)] = abs((p(rpExamined) - urObserved) / urObserved)
        #spline fit
        tck = interpolate.splrep(dm.ix[:iRg,read_pairs_for_estimation].values, dm.ix[:iRg,'UNIQUE_READS'].values, s=0,k=sdegree)
        yspline = interpolate.splev(dm.ix[:,read_pairs_for_estimation].values, tck, der=0)
        dfError.ix[sample,'spline_perc_err_deg' +str(sdegree)] = abs((yspline[-1] - urObserved) / urObserved)
        ### make plot of polynomail fit
        xp = np.linspace(0, dm[read_pairs_for_estimation].max(), 100)
        polylabel= 'polynomial fit, deg= ' + str(sdegree)
        plt.plot(xp,p(xp),'g-',label=polylabel)
        plt.plot(dm[read_pairs_for_estimation],dm['UNIQUE_READS'],'-ob',label='observed data')
        plt.plot(dm[read_pairs_for_estimation],dm['ROI_' + iRg + '_est'],'-.r',label='roi-based estimate')
        plt.plot([dm.ix[iRg,read_pairs_for_estimation]],[dm.ix[iRg,'UNIQUE_READS']],'ro',label='extrapolation point')
        plt.grid()
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.legend(loc='lower right',numpoints=1)
        plt.title('library complexity: WGS sample ' + sample)
        plt.ylabel('unique reads')
        outF = wkdir + '/'+ sample +'_' + iRg + '_poly_curve_deg' + str(sdegree) +'.png'
        plt.savefig(outF, bbox_inches='tight')
        plt.close()
        ### make plot of spline fit
        slabel= 'spline fit, deg= ' + str(sdegree)
        plt.plot(dm.ix[:,read_pairs_for_estimation].values,yspline,'g-',label=slabel)
        plt.plot(dm[read_pairs_for_estimation],dm['UNIQUE_READS'],'-ob',label='observed data')
        plt.plot(dm[read_pairs_for_estimation],dm['ROI_' + iRg + '_est'],'-.r',label='roi-based estimate')
        plt.plot([dm.ix[iRg,read_pairs_for_estimation]],[dm.ix[iRg,'UNIQUE_READS']],'ro',label='extrapolation point')
        plt.grid()
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.legend(loc='lower right',numpoints=1)
        plt.title('library complexity: WGS sample ' + sample)
        plt.ylabel('unique reads')
        plt.xlabel(read_pairs_for_estimation)
        outF = wkdir + '/'+ sample +'_' + iRg + '_spline_curve_deg' + str(sdegree) +'.png'
        plt.savefig(outF, bbox_inches='tight')
        plt.close()
### make plot of error rates for multiple models
outF = wkdir + '/polynomial_spline_percent_error.png'
plt.boxplot(dfError.values)
plt.xticks(np.arange(1,len(colSet)+2), colSet, rotation=90)
plt.ylabel('error')
plt.title('unique read estimates \n extrapolated from ' + iRg + ' read groups')
plt.savefig(outF, bbox_inches='tight')
plt.close()


### calculate fraction of optical duplicates
dataDict = {}
for sample in samples:
    sample = 'NexPond-' + sample 
    print sample 
    outdir = out_base + '/' + sample
    fin2 = outdir + '/' + sample + '_chrm21_read_grp_dup_metrics.txt'
    dm = pd.read_csv(fin2,index_col=0) #load in yield estimates
    dm['fraction_of_dups_that_are_optical'] = dm.READ_PAIR_OPTICAL_DUPLICATES/ dm.READ_PAIR_DUPLICATES
    dm['fraction_optical_dup'] = dm.READ_PAIR_OPTICAL_DUPLICATES/ dm.READ_PAIRS_EXAMINED
    dm['fraction_PCR_dup'] = (dm.READ_PAIR_DUPLICATES - dm.READ_PAIR_OPTICAL_DUPLICATES)/ dm.READ_PAIRS_EXAMINED
    #dataDict[sample] = dm
    plt.plot(dm.READ_PAIRS_EXAMINED,dm.fraction_PCR_dup)
    # plt.plot(dm.READ_PAIRS_EXAMINED,dm.fraction_of_dups_that_are_optical)
    # plt.plot(dm.READ_PAIRS_EXAMINED,dm.fraction_optical_dup)
# dmMtrx = pd.Panel(dataDict)
# dmMtrx.ix[:,:,'fraction_optical_dup']
# plt.plot(dm.READ_PAIRS_EXAMINED,dm.fraction_of_dups_that_are_optical)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
# plt.title('fraction of duplicate events that are optical')
# plt.title('fraction of read pairs that are optical duplicates')
plt.title('fraction of read pairs that are PCR duplicates')
plt.xlabel('reads observed')
# plt.ylabel('optical fraction')
plt.ylabel('PCR fraction')
outF = out_base + '/nexome_fraction_PCR_dup.png'
# outF = out_base + '/nexome_fraction_of_dups_that_are_optical.png'
# outF = out_base + '/nexome_fraction_optical_dup.png'
plt.savefig(outF, bbox_inches='tight')
plt.close()


#######################
### run preseq call ###
#######################

wkdir = out_base + '/preseq_curve_estimates' 
if (not os.path.isdir(wkdir)):
    os.mkdir(wkdir)
for sample in samples:
    print sample
    sample = 'NexPond-' + sample 
    outdir = out_base + '/' + sample
    ccurve  = wkdir +'/' + sample+'_rgSet' + str(16) + '.c_curve'
    lcextrap = wkdir +'/' + sample+'_rgSet' + str(3) + '.lc_extrap'
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
    plt.title('preseq yield estimate: based on 2 read \n groups for Nexome sample 376014')
    plt.ylabel('expected distinct')
    plt.xlabel('reads observed')
    outF = wkdir  +'/' + sample+'_rgSet' + str(i) + '_c_curve.png'
    plt.savefig(outF, bbox_inches='tight')
    plt.close()

# ### scratch - fit general exponential 
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.optimize import curve_fit

# def func(x, a, b, c):
#     return a * np.exp(-b * x) + c

# def func2(x, a, b, c):
#     return a * (1-np.exp(-(x*b)/a))/c

# def estimateRoi(EstLibSize, x, n, c):
#     # EstLibSize = estimated library size
#     # x = multiple of sequencing to be simulated (how many X sequencing)
#     # n = read pairs observed
#     # c = unique read pairs 
#     return EstLibSize * (1 - np.exp(-(x*n)/EstLibSize))/c

# #c/float(x) - 1 + np.exp(-n/float(x))

# x = np.linspace(0,4,50)
# # y = func(x, 1, 1, 1)
# y = estimateRoi(1, x, 5, 1)
# yn = y + 0.01*np.random.normal(size=len(x))

# # popt, pcov = curve_fit(func, x, yn)
# popt, pcov = curve_fit(estimateRoi, x, yn)

# plt.figure()
# plt.plot(x, yn, 'ko', label="Original Noised Data")
# # plt.plot(x, func(x, *popt), 'r-', label="Fitted Curve")
# plt.plot(x, estimateRoi(x, *popt), 'r-', label="Fitted Curve")
# plt.legend()
# plt.show()


# ### scratch 
# x = np.linspace(0,4,50)
# # y = func(x, 1, 1, 1)
# y1 = estimateRoi(1000, x, 2000, 200)
# y2 = estimateRoi(1000, x, 2000, 400)

# plt.plot(x, y1, 'ro-', label="1")
# plt.plot(x, y2, 'bo-', label="2")
# plt.legend()
# plt.show()



# ### interpolate 3 points
# x = dm.ix[:iRg,read_pairs_for_estimation].values
# y = dm.ix[:iRg,'UNIQUE_READS'].values
# popt, pcov = curve_fit(estimateRoi, x, y)
# plt.plot(dm.ix[:,read_pairs_for_estimation].values, 
#         dm.ix[:,'UNIQUE_READS'].values, 
#         'ko', label="Original Noised Data")
# # plt.plot(x, func(x, *popt), 'r-', label="Fitted Curve")
# plt.plot(dm.ix[:,read_pairs_for_estimation].values, 
#         estimateRoi(dm.ix[:,read_pairs_for_estimation].values, 
#             *popt), 'r-', label="Fitted Curve")
# outF = out_base + '/nexome_scratch.png'
# plt.savefig(outF, bbox_inches='tight')
# plt.close()



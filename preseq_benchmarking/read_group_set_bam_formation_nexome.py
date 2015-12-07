import numpy as np
import pandas as pd
import os
import subprocess
import time
#use .matplotlib-1.3.1-python-2.7.1-sqlite3-rtrees
#use .scipy-1.3.1-python-2.7.1-sqlite3-rtrees

sample_root = '/seq/picard_aggregation/D5227'
out_base = '/humgen/gsa-hpprojects/dev/hogstrom/depth_by_read_group/Nexome'
chrm_bed_file = '/humgen/gsa-hpprojects/dev/hogstrom/depth_by_read_group/chrm21_full.bed'
picard_path = '/seq/software/picard/current/bin/picard.jar'
#jl_path = '/humgen/gsa-hpprojects/dev/hogstrom/code/parallel_preseq/read_counts/preseq_lh.jl'
jl_path = '/humgen/gsa-hpprojects/dev/hogstrom/code/parallel_preseq/read_counts/preseq_lh_parallel.jl'

processes = set()
max_processes = 27

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

#/seq/picard_aggregation/D5227/NexPond-376014/current
######################################################
### write chromosome specific file for each sample ###
######################################################

for sample in samples:
    sample = 'NexPond-' + sample 
    outdir = out_base + '/' + sample 
    if (not os.path.isdir(outdir)):
        os.mkdir(outdir)
    # bam_path = os.path.join(sample_root,sample,sample+'.bam')
    bam_path = sample_root + '/' + sample + '/current/' + sample+ '.bam'
    outChrmbam  = outdir +'/' + sample+'_chrm21.bam'
    cmd_bam = ' '.join(['samtools view',
                         '-h \"' + bam_path+ '\"',
                         '-L ' + chrm_bed_file,
                         '-o ' + outChrmbam])
    processes.add(subprocess.Popen(cmd_bam,shell=True))
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update(
            p for p in processes if p.poll() is not None)


#samtools view -h "/seq/picard_aggregation/D5227/NexPond-359781/current/NexPond-359781.bam" -L /humgen/gsa-hpprojects/dev/hogstrom/depth_by_read_group/chrm18_to_21_full.bed -o /humgen/gsa-hpprojects/dev/hogstrom/depth_by_read_group/Nexome/NexPond-359781/NexPond-359781_chrm18_to_21.bam

##############################################
### define read group sets for each sample ###
##############################################
processes = set()
for sample in samples:
    print sample
    sample = 'NexPond-' + sample 
    outdir = out_base + '/' + sample 
    if (not os.path.isdir(outdir)):
        os.mkdir(outdir)
    # bam_path = sample_root +'/' + sample +'/' + sample+ '.bam'
    bam_path = outdir +'/' + sample+'_chrm21.bam'
    # make comand line call for sametools
    cmd = ' '.join(['samtools view -H',
                         '\"' + bam_path + '\"',
                         '| grep \'^@RG\''])
    gout = os.popen(cmd).read()
    rgLines = gout.split('@RG')
    rgLines.pop(0)
    rgs = [rg[4:11] for rg in rgLines]
    ### make read group sets + write new bam file for each
    rgSer = pd.Series(rgs)
    i=rgSer.shape[0]
    while (i >= 1):
        print i
        ind = range(i)
        rgSubset = rgSer.ix[ind]
        # outRGtxt = os.path(outdir, sample+'_rgSet' + str(i)) + '.txt')
        outRGtxt = outdir +'/' + sample+'_rgSet' + str(i) + '.txt'
        rgSubset.to_csv(outRGtxt,index=False)
        # i = i-1
        # create bam file for individual chromsome
        # outRGbam = os.path(outdir, sample+'_rgSet' + str(i)) + '.bam')
        outRGbam  = outdir +'/' + sample+'_chrm21_rgSet' + str(i) + '.bam'
        cmd_bam = ' '.join(['samtools view',
                             '-h \"' + bam_path+ '\"',
                             '-L ' + chrm_bed_file,
                             '-R ' + outRGtxt,
                             '-o ' + outRGbam])
        processes.add(subprocess.Popen(cmd_bam,shell=True))
        i = i-1
        if len(processes) >= max_processes:
            # time.sleep(2)
            os.wait()
            processes_temp = processes.copy()
            processes.difference_update(
                p for p in processes_temp if p.poll() is not None)

############################################
### run mark Dup for each read group set ###
############################################
processes = set()
for sample in samples:
    print sample
    sample = 'NexPond-' + sample 
    outdir = out_base + '/' + sample
    bam_path = outdir +'/' + sample+'_chrm21.bam'
    cmd = ' '.join(['samtools view -H',
                         '\"' + bam_path + '\"',
                         '| grep \'^@RG\''])
    gout = os.popen(cmd).read()
    rgLines = gout.split('@RG')
    rgLines.pop(0)
    rgs = [rg[4:11] for rg in rgLines]
    ### make read group sets + write new bam file for each
    rgSer = pd.Series(rgs)
    i=rgSer.shape[0]
    while (i >= 1):
        print i
        inRGbam  = outdir +'/' + sample+'_chrm21_rgSet' + str(i) + '.bam'
        outRGbam  = outdir +'/' + sample+'_chrm21_rgSet' + str(i) + '_dup_marked.bam'
        # mark duplicates
        outRGmetrics  = outdir +'/' + sample+'_rgSet' + str(i) + '.duplication_metrics'
        cmd2 = ' '.join(['java -jar',
                             picard_path,
                             'MarkDuplicates',
                             'INPUT='+inRGbam,
                             'OUTPUT='+outRGbam,
                             'METRICS_FILE='+outRGmetrics])
        processes.add(subprocess.Popen(cmd2,shell=True))
        if len(processes) >= max_processes:
            os.wait()
            processes_temp = processes.copy()
            processes.difference_update(
                p for p in processes_temp if p.poll() is not None)
        i = i-1

### run DuplicationByInsertLength
processes = set()
max_processes = 5
picard_path='/home/unix/hogstrom/picard_dup_by_insert.jar' #modified jar
for sample in samples:
    print sample
    sample = 'NexPond-' + sample 
    outdir = out_base + '/' + sample
    bamName = sample+'_chrm21.bam'
    inRGbam = outdir +'/' + bamName
    fPrefix = bamName.split('.')[0]
    outMetric = outdir +'/' + fPrefix +'.insertSize_by_dup'
    cmd2 = ' '.join(['java -jar',
                         picard_path,
                         'DuplicationByInsertLength',
                         'INPUT='+inRGbam,
                         'OUTPUT='+outMetric,
                         'CHART=/home/unix/hogstrom/null.txt'])
    processes.add(subprocess.Popen(cmd2,shell=True))
    if len(processes) >= max_processes:
        os.wait()
        processes_temp = processes.copy()
        processes.difference_update(
            p for p in processes_temp if p.poll() is not None)

#######################
### run preseq call ###
#######################

def get_rgs(fpath):
    cmd = ' '.join(['samtools view -H',
                     '\"' + fpath + '\"',
                     '| grep \'^@RG\''])
    gout = os.popen(cmd).read()
    rgLines = gout.split('@RG')
    rgLines.pop(0)
    return [rg[4:11] for rg in rgLines]

processes = set()
max_processes = 13
# wkdir = out_base + '/preseq_curve_estimates' 
wkdir = out_base + '/parallel_preseq_27Nov'
if (not os.path.isdir(wkdir)):
    os.mkdir(wkdir)
for sample in samples:
    print sample
    sample = 'NexPond-' + sample 
    outdir = out_base + '/' + sample
    #RGbam  = outdir +'/' + sample+'_chrm21_rgSet' + str(i) + '.bam'
    ### run c_curve on whole file
    preseq_mode = 'c_curve'
    RGfull  = outdir +'/' + sample+'_chrm21.bam' # 
    i = len(get_rgs(RGfull))
    bamName = sample+'_chrm21_rgSet' + str(i) + '_dup_marked.bam' # largest read group set
    RGbam  = outdir +'/' + bamName
    # outRGmetrics  = wkdir +'/' + sample+'_chrm21' + '.' + preseq_mode
    # outlog  = wkdir +'/' + sample+'_chrm21' + '_' + preseq_mode + '.log'
    # cmd2 = ' '.join(['/humgen/gsa-hpprojects/dev/hogstrom/code/preseq/preseq',
    #                      preseq_mode + ' -v', #-P
    #                      '-B ' + RGbam,
    #                      '-s 1e+04',
    #                      '-o '+ outRGmetrics,
    #                      '&> ' + outlog])
    # # gout = os.popen(cmd2).read()
    # processes.add(subprocess.Popen(cmd2,shell=True))
    # if len(processes) >= max_processes:
    #     os.wait()
    #     processes_temp = processes.copy()
    #     processes.difference_update(
    #         p for p in processes_temp if p.poll() is not None)
    ### run preseq_jl
    fPrefix = bamName.split('.')[0]
    outMetric = wkdir +'/' + fPrefix +'.jl_read_counts'
    if not os.path.exists(outMetric):
        cmd2 = ' '.join(['julia ' + jl_path,
                             '--bam '+RGbam,
                             '--output '+outMetric])
        processes.add(subprocess.Popen(cmd2,shell=True))
        if len(processes) >= max_processes:
            os.wait()
            processes_temp = processes.copy()
            processes.difference_update(
                p for p in processes_temp if p.poll() is not None)
    ### run lc_extrap on smaller set of read groups
    i = 3 # use i number of read groups
    preseq_mode = 'lc_extrap'
    RGbam  = outdir +'/' + sample+'_chrm21_rgSet' + str(i) + '_dup_marked.bam' # mark duplicates
    outRGmetrics  = wkdir +'/' + sample+'_rgSet' + str(i) + '.' + preseq_mode
    outlog  = wkdir +'/' + sample+'_rgSet' + str(i) + '_' + preseq_mode + '.log'
    cmd3 = ' '.join(['/humgen/gsa-hpprojects/dev/hogstrom/code/preseq/preseq',
                         preseq_mode + ' -v', #-P
                         '-B ' + RGbam,
                         '-s 1e+04',
                         '-e 1e+06',
                         '-o '+ outRGmetrics,
                         '&> ' + outlog])
    # gout = os.popen(cmd3).read()
    processes.add(subprocess.Popen(cmd3,shell=True))
    if len(processes) >= max_processes:
        os.wait()
        processes_temp = processes.copy()
        processes.difference_update(
            p for p in processes_temp if p.poll() is not None)


################################
### run parallel preseq call ###
################################

def get_rgs(fpath):
    cmd = ' '.join(['samtools view -H',
                     '\"' + fpath + '\"',
                     '| grep \'^@RG\''])
    gout = os.popen(cmd).read()
    rgLines = gout.split('@RG')
    rgLines.pop(0)
    return [rg[4:11] for rg in rgLines]

wkdir = out_base + '/parallel_preseq_27Nov'
sample = 'NexPond-359781'
nprocs = 16
for nproc in range(1,17):
    print nproc
    #sample = 'NexPond-' + sample 
    outdir = out_base + '/' + sample
    #RGbam  = outdir +'/' + sample+'_chrm21_rgSet' + str(i) + '.bam'
    ### run c_curve on whole file
    preseq_mode = 'c_curve'
    RGfull  = outdir +'/' + sample+'_chrm21.bam' # 
    i = len(get_rgs(RGfull))
    #bamName = sample+'_chrm21_rgSet' + str(i) + '_dup_marked.bam' # largest read group set
    bamName = 'NexPond-359781_chrm21_rgSet8_dup_marked.bam'
    RGbam  = outdir +'/' + bamName
    fPrefix = bamName.split('.')[0]
    outMetric = wkdir +'/' + fPrefix +'.jl_read_counts'
    log_file =  wkdir +'/' + fPrefix +'.processor_timing'
    cmd2 = ' '.join(['julia -p ' + str(nproc),
                         jl_path,
                         '--bam '+RGbam,
                         '--output '+outMetric,
                         '>> ' + log_file])
    os.popen(cmd2)

### run preseq to estimate libary yield
# /humgen/gsa-hpprojects/dev/hogstrom/code/preseq/preseq lc_extrap -B \
# -s 1e+04 \
# -e 1e+06 \
# -o $outdir/NexPond-376014_chrm21_rg3_yield_estimates_withoutPairedEndFlag2.txt \
# $outdir/NexPond-376014_chrm21_rg3_dupMarked.bam \
# -v

### preseq complexity curve
# /humgen/gsa-hpprojects/dev/hogstrom/code/preseq/preseq c_curve -B \
# -s 1e+04 \
# -o $outdir/NexPond-376014_chrm21_rg2_complexity_curve.txt \
#  $outdir/NexPond-376014_chrm21_rg2_dupMarked.bam \
# -p \
# -v

import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import os 
import subprocess
#use .matplotlib-1.3.1-python-2.7.1-sqlite3-rtrees
#use Samtools
#use .julia-0.4.0

sample_root = '/seq/picard_aggregation/D5227'
out_base = '/humgen/gsa-hpprojects/dev/hogstrom/depth_by_read_group/Nexome'
# picard_path='/home/unix/hogstrom/picard_dup_by_insert.jar' #modified jar
jl_path = '/humgen/gsa-hpprojects/dev/hogstrom/code/parallel_preseq/read_counts/preseq.jl'

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

### run Julia 
processes = set()
max_processes = 13
for sample in samples:
    print sample
    sample = 'NexPond-' + sample 
    outdir = out_base + '/' + sample
    wkdir = out_base + '/parallel_preseq_27Nov'
    # bamName = sample+'_chrm21.bam' 
    bamName = sample+'_chrm21_rgSet3_dup_marked.bam' 
    inRGbam = outdir +'/' + bamName
    fPrefix = bamName.split('.')[0]
    outMetric = wkdir +'/' + fPrefix +'.jl_read_counts'
    if not os.path.exists(outMetric):
        cmd2 = ' '.join(['julia ' + jl_path,
                             '--bam '+inRGbam,
                             '--output '+outMetric])
        processes.add(subprocess.Popen(cmd2,shell=True))
        if len(processes) >= max_processes:
            os.wait()
            processes_temp = processes.copy()
            processes.difference_update(
                p for p in processes_temp if p.poll() is not None)

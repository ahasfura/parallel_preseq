## Read count calculation

Run on Ubuntu 14.10. Install pre-requisites (like GSL), compile samtools and preseq from source. 

### Possibility of speedup from pre-calculating read counts

Using test_chrm21_rg2.bam, a BAM file containing entries from chromosome 21, the complexity curve can be calculated using the preseq tool (along with execution time) by running:

```bash
time ./preseq c_curve -o complexity_from_bam.txt -B test_chrm21_rg2.bam -s 10000
```

The execution time for calculating the complexity curve using the .bam file is 1.4 s. We can test how much the execution time is reduced if we provide the preseq tool with pre-calculated read counts. To start off, we can just dump these read counts while they are being calculated by preseq in the above command by adding 

```c
std::cout << current_count << "\n";
```

towards the end of the second if-block of the update_se_duplicate_counts_hist function of load_data_for_complexity.cpp, and directing the output to a file (counts_preseq.txt). If we re-run preseq on this,

```bash
time ./preseq c_curve -o complexity_from_preseq_counts.txt -V counts_preseq.txt -s 10000
```

we get an execution time of 0.04 s, a factor of 30 speedup. So calculation of read counts really is a significant bottleneck in the complexity curve calculation.

### How to go about calculating read counts

The preseq manual has an example bash command for calculating read counts from a FASTQ file. We can try it out after converting our .bam file to a .fastq file using samtools:

```bash
samtools bam2fq test_chrm21_rg2.bam > test_chrm21_rg2.fastq
awk '{if (NR%4==2) print substr($0,1,20)}' test_chrm21_rg2.fastq | sort | uniq -c | awk '{print $1}' > counts_fastq.txt
```

This does not really give the correct complexity curve we pass it to preseq, 

```bash
time ./preseq c_curve -o complexity_from_fastq_counts.txt -V counts_fastq.txt -s 10000
```

so we need to look into what is going on here. In the end, this is probably not the way to go for calculating read counts, since conversion from .bam to .fastq may end up being quite expensive.

A better option is to re-implement the load_counts_BAM_se (and it's helper update_se_duplicate_counts_hist) from load_data_for_complexity.cpp in Julia, and parallelizing them to introduce speedup. Since these methods use helper methods from Smith lab C++ codes like GenomicRegion.cpp and MappedRead.cpp in order to parse BAM file entries, we need to figure out how to call C++ code from Julia. There is discussion on this here: docs.julialang.org/en/release-0.4/manual/calling-c-and-fortran-code

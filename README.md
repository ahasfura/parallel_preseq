# parallel_preseq
parallelizing the preseq calculation for DNA sequencing.

# library size estimates & predicted yield curve from preseq:

![](figures/376014_chrm21_yield_estimate_curve.png)

![](figures/376014_chrm21_observed_yield_curve.png)

yield curves for 3 read groups with Paired end read flag on:
![](figures/376014_chrm21_observed_yield_curve_rg3_withPairedEndFlag.png)

and paired end read flag off:
![](figures/376014_chrm21_observed_yield_curve_rg3_withoutPairedEndFlag.png)
*It seems strange that the predicted yeild curve got even the point for the input file wrong. Perhaps this command isn't being run correclty or on the right file.
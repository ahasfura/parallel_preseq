
using PyPlot

#################
#  Create Data  #
#################
x = randn(1000) # Values
nbins = 50 # Number of bins

println(typeof(x))

##########
#  Plot  #
##########
fig = figure("pyplot_histogram",figsize=(10,10)) # Not strictly required
ax = axes() # Not strictly required
h = plt[:hist](x,nbins) # Histogram

grid("on")
xlabel("X")
ylabel("Y")
title("Histogram")

savefig("pyplot_historgram.svg")

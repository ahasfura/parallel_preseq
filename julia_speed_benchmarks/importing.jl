using PyPlot

#=
modules = readdir(string(homedir(), "/.julia/lib/v0.4/"))
for mod in modules
  println(string(mod))
  tic()
  using mod[1:length(mod)-3]
  push!(times, toc())
end
=#
times = []
nbins = 2

println("DistributedArrays")
tic()
using DistributedArrays
push!(times, toc())

println("Interact")
tic()
using Interact
push!(times, toc())

println("Colors")
tic()
using Colors
push!(times, toc())

println("Clustering")
tic()
using Clustering
push!(times, toc())

println("Distributions")
tic()
using Distributions
push!(times, toc())

println("IJulia")
tic()
using IJulia
push!(times, toc())

println("JSON")
tic()
using JSON
push!(times, toc())

println("SortingAlgos")
tic()
using SortingAlgorithms
push!(times, toc())

println("DataFrames")
tic()
using DataFrames
push!(times, toc())

fig = figure("importing_hist", figsize=(10, 10))
ax = axes()

println(times)
h = plt[:hist](times, nbins) # Histogram

grid("on")
xlabel("time")
title("Importing Times")

savefig("hist_of_import_times.svg")

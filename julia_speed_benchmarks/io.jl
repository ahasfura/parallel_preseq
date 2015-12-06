using PyPlot


sizes = [10^i for i in range(1,9)]
read_times = []
write_times = []

for size in sizes
  println(size)

  #write portion
  tic()
  f = open("test.txt", "w")
  for i in range(0,size)
    write(f, "speed test")
  end
  close(f)
  push!(write_times, toc())

  #read portion
  tic()
  f = open("test.txt", "r")
  lines = readall(f)
  #=
  for line in lines
    i += 1
  end
  =#
  close(f)
  push!(read_times, toc())

end

println(read_times)
println(write_times)

log_sizes = [log10(x) for x in sizes]
log_read_times = [log10(x) for x in read_times]
log_write_times = [log10(x) for x in write_times]

plot(log_sizes, log_write_times)

title("Write Times vs File Size")
xlabel("Log File Size")
ylabel("Log Write Times")

savefig("write_times_vs_file_size.svg")

using PyPlot

x = [.372, 1.88, 3.78, 18.24, 36.04]
y = [.432, 1.07, 2.10, 6.05, 26.21]

plot(x,y)
title("Run Time Difference Vs. Execultable Run Time")
xlabel("Executable Run Time")
ylabel("Run Time Difference")

savefig("exec_run_time_diff.svg")

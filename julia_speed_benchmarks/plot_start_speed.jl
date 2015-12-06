using PyPlot

x = [1, 2, 3, 4, 5, 6, 7, 8]
y = [3.284, 4.217, 5.127, 5.726, 6.453, 7.738, 8.757, 9.183]

plot(x,y)
title("Start Time Vs. Num Cores")
xlabel("Num Cores")
ylabel("Start Time")

savefig("start_speed_vs_num_cores.svg")


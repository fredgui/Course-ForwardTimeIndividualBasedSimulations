

rm(list=ls())


# PLOT TEMPORAL DYNAMICS ................................................... ----
# info                    ----
# -> Lets first explore the temporal dynamics of evolution at the example of a single simulation with no plasticity g1=0
# -> to do so:
#    1) adapt the path in line 19 to your computer
#    2) mark all lines between 19 and 65, run them, and look at the resulting three graphs 
#    3) ProTip: you can collapse single sections when clicking on the arrows beside the line numbers
#    4) now you can switch to the simulation results with another resource distribution (e.g., uncomment lines 23 and 24), 
#       and repeat the plotting
#    5) add new resource distributions to line 45 of "Exercise_Resource_competition.ini", and rerun simulations
#       note that you need to adapt the file path in line 20-21 and 87 of this R script when plotting a different plasticity scenario

# load simulation data    ####
setwd("/path/to/7_Competition")
data          <- read.table("Competition_q_{{1,2,3,4,5}{1,2,3,4,5}}_bygen.txt", header=T)
data_genotype <- read.table("Competition_q_{{1,2,3,4,5}{1,2,3,4,5}}_2000_1.quanti", header=T)          # t=2000, replicate 1

# data          <- read.table("resource_q_{{0,1.5,3,4.5,6}{0,1.5,3,4.5,6}}_bygen.txt", header=T)
# data_genotype <- read.table("resource_q_{{0,1.5,3,4.5,6}{0,1.5,3,4.5,6}}_2000_1.quanti", header=T)          # t=2000, replicate 1

# resource variation      ----
# 1. res distr ...............
res_freq   <- c(0.2, 0.2, 0.2, 0.2, 0.2)    # the frequencies, line XY of ini-file
res_prop   <- c(1,2,3,4,5)              # resource properties, line XY of ini-file
res_mean   <- mean(res_prop)
res_var_w1 <- sum( res_freq*((res_prop - res_mean)^2 ) )
res_var_w1

par(mfrow=c(1,1))
plot(x=res_prop, y=res_freq,pch=16, lwd=8, xlim=c(-2,7), las=1, type="h", main="resource distribution", xlab="resource property", ylab="resource frequency")

# 2. res distr ...............
res_freq   <- c(0.2, 0.2, 0.2, 0.2, 0.2)    # the frequencies, line XY of ini-file
res_prop   <- c(0,1.5,3,4.5,6)              # resource properties, line XY of ini-file
res_mean   <- mean(res_prop)
res_var_w2 <- sum( res_freq*((res_prop - res_mean)^2 ) )
res_var_w2

par(mfrow=c(1,1))
plot(x=res_prop, y=res_freq,pch=16, lwd=8, xlim=c(-2,7), las=1, type="h", main="resource distribution", xlab="resource property", ylab="resource frequency")


# plot avg. z             ----
par(mfrow=c(1,3))
plot(x=data$generation, y=data$adlt.q1.p1, lwd=2, col="transparent", xlab="time", ylab="z", main="avg. consumer trait, z", las=1, ylim=c(0,5), type="l")
grid()
lines(x=data$generation, y=data$adlt.q1.p1, lwd=2, col="skyblue")
lines(x=data$generation, y=data$adlt.q1.p2, lwd=2, col="coral")

# plot var. z             ----
plot(x=data$generation, y=data$adlt.Va.q1.p1, lwd=2, col="transparent", xlab="time", ylab="VA", main="additive genetic variance, VA", las=1, ylim=c(0,6), type="l")
grid()
lines(x=data$generation, y=data$adlt.Va.q1.p1, lwd=2, col="skyblue")
lines(x=data$generation, y=data$adlt.Va.q1.p2, lwd=2, col="coral")

# plot hist z             ----
hist(data_genotype$P1, col="skyblue", xlab="trait value, z", las=1, xlim=c(-3,9), breaks=seq(-10,10,0.2), main="trait distribution")
par(mfrow=c(1,1))

#                         ----
# PLOT EQUILIBRIUM ......................................................... ----
# info                    ----
# -> In a second step, you can plot genetic variance as a function of resource variation
# -> it makes sense to use VA at evolutionary equilibrium (e.g., after 2,000 years), when evolution reached a plateau
# -> to this end:
#    1) adapt the path in line 19 to your computer
#    2) run lines 86-105, and make a first overview plot
#    3) then extend your ini-file for more resource distributions (adding additional matrices to line 45 of the ini-file), and run these simulations
#    4) add these distributions to the Rscript (in sections "a) compute o_w^2", "b) VA from simulations", and "c) create summary table")
#    6) again, you now can change a parameter in the ini-file (e.g., competition being absent or present), rerun the simulations, and repeat this kind of plot

# a) compute o_w^2        ----
res_freq <- c(0.2, 0.2, 0.2, 0.2, 0.2)    # frequencies

# 1. distribution .............
res_prop1   <- c(1,2,3,4,5)                                 # resource properties
res_mean1   <- mean(res_prop1)                              # mean properties 
res_var_w_1 <- sum( res_freq*((res_prop1 - res_mean1)^2 ) ) # resource variance

# 2. distribution .............
res_prop2   <- c(0,1.5,3,4.5,6)                             # resource properties
res_mean2   <- mean(res_prop2)                              # mean properties 
res_var_w_2 <- sum( res_freq*((res_prop2 - res_mean2)^2 ) ) # resource variance

# b) VA from simulations  ----
# 1. distribution .............
data1 <- read.table("Competition_q_{{1,2,3,4,5}{1,2,3,4,5}}_bygen.txt", header=T)[2000,"adlt.q1.Va"]

# 2. distribution .............
data2 <- read.table("Competition_q_{{0,1.5,3,4.5,6}{0,1.5,3,4.5,6}}_bygen.txt", header=T)[2000,"adlt.q1.Va"]

# c) create summary table ----
summaryData           <- c(data1, data2)
summaryData           <- cbind(summaryData,c(res_var_w_1, res_var_w_2))
colnames(summaryData) <- c("VA", "ow")

# d) plot                 ----
par(mfrow=c(1,1))
plot(x=summaryData[,"ow"], y=summaryData[,"VA"], ylim=c(0,7), type="o", lwd=2, xlab=expression(paste("within-patch resource variation, ",sigma[w]^2,sep="")), main="genetic variance, VA", ylab="VA", las=1)
par(mfrow=c(1,1))

#                         ----




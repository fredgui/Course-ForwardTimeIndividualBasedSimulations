

rm(list=ls())


# PLOT TEMPORAL DYNAMICS ................................................... ----
# info         ----
# -> Lets first explore the temporal dynamics of evolution at the example of a single simulation with no plasticity g1=0
# -> to do so:
#    1) adapt the path in line 19 to your computer
#    2) mark all lines between 19 and 75, run them, and look at the resulting four graphs 
#    3) ProTip: you can collapse single sections when clicking on the arrows beside the line numbers on the left
#    4) now you can change a parameter in the ini-file (e.g., dispersal rate or selection strength), 
#       rerun the simulations, and repeat the plotting
#    5) add new plasticity levels to line 47 of "Exercise_Phenotypic_plasticity.ini"
#       note that you need to add/adapt the file path in line 21 and line 89 of this R script when plotting a different plasticity scenario

# load data    ####
setwd("/path/to/directory/6_Plasticity")

data <- read.table("Plasticity_g1_{{0}}_bygen.txt", header=T)
# data <- read.table("Plasticity_g1_{{0.5}}_bygen.txt", header=T)

# plot avg. z  ----
par(mfrow=c(2,2))
plot(x=data$generation, y=data$adlt.z1.p1, lwd=2, col="transparent", xlab="time", ylab="phenotype, z", main="trait value, z", las=1, ylim=c(9,13), type="l")
grid()
lines(x=data$generation, y=data$adlt.z1.p1, lwd=2, col="skyblue")
lines(x=data$generation, y=data$adlt.z1.p2, lwd=2, col="coral")
legend(x=400, y=13, legend=c("patch 1", "patch 2"), lty=1, col=c("skyblue", "coral"), lwd=2)

# plot avg. g0 ----
plot(x=data$generation, y=data$adlt.g0_1.p1, lwd=2, col="transparent", xlab="time", ylab="intercept, g0", main="breeding values, g0", las=1, ylim=c(9,13), type="l")
grid()
lines(x=data$generation, y=data$adlt.g0_1.p1, lwd=2, col="skyblue")
lines(x=data$generation, y=data$adlt.g0_1.p2, lwd=2, col="coral")

# plot avg. g1 ----
plot(x=data$generation, y=data$adlt.g1_1.p1, lwd=2, col="skyblue", xlab="time", ylab="slope, g1", main="plasticity, g1", las=1, ylim=c(-1,1), type="l")
grid()
lines(x=data$generation, y=data$adlt.g1_1.p1, lwd=2, col="skyblue")
lines(x=data$generation, y=data$adlt.g1_1.p2, lwd=2, col="coral")

# plot NoR     ----
data_end <- data[1000,]

# g0-values
g0_p1 <- data_end$adlt.g0_1.p1 # patch 1
g0_p2 <- data_end$adlt.g0_1.p2 # patch 2

# g1-values
g1_p1 <- data_end$adlt.g1_1.p1 # patch 1
g1_p2 <- data_end$adlt.g1_1.p2 # patch 2

# e-values
e_p1 <- data_end$adlt.e1.p1 # patch 1
e_p2 <- data_end$adlt.e1.p2 # patch 2

# z-values 
z_p1 <- data_end$adlt.z1.p1 # patch 1
z_p2 <- data_end$adlt.z1.p2 # patch 2

# plot reaction norms
evalues   <- seq(-1,1,0.01)
z1_p1_nor <- g0_p1 + g1_p1 * evalues
z1_p2_nor <- g0_p2 + g1_p2 * evalues

plot( x=evalues, y=z1_p1_nor, lwd=3, col="skyblue", ylim=c(8,14), xlim=c(-1.2,1.2), las=1, ylab="phenotype", xlab="environment", main="reaction norm at t=1000", type="l")
lines(x=evalues, y=z1_p2_nor, lwd=3, col="coral")
points(x=c(-1,1), y=c(z_p1, z_p2), col=c("skyblue", "coral"), lwd=3)
points(x=c(-1,1), y=c(10,12), col=c("skyblue", "coral"), lwd=3, pch=4)
par(mfrow=c(1,1))


#              ----
# PLOT DIVERGENCE AT EQUILIBRIUM ........................................... ----
# info         ----
# -> In a second step, you can plot genetic and phenotypic divergence as a function of plasticity
# -> it makes sense to use Dg and DP at evolutionary equilibrium (e.g., after 1,000 year), when dispersal and selection balanced out (at least to a great extent)
# -> to this end:
#    1) adapt the path in line 15 to your computer
#    2) extend your simulations for more plasticity values (adding additional values to line 47 of the ini-file)
#    3) run these simulations
#    4) add these g1-values to line 87 of this script
#    5) run lines 86-105, and make an overview plot
#    6) again, you now can change a parameter in the ini-file (e.g., dispersal rate), rerun the simulations, and repeat this kind of plot

# summary table----
g1_value    <- c(0, 0.5)
summaryData <- NULL

for(g1 in g1_value){ # g1 <- -1
  
  data        <- read.table(paste("Plasticity_g1_{{",g1,"}}_bygen.txt",sep=""), header=T)[1000,]
  Dg          <- data$adlt.g0_1.p2 - data$adlt.g0_1.p1
  Dp          <- data$adlt.z1.p2  - data$adlt.z1.p1
  summaryData <- rbind(summaryData, c(g1, Dg, Dp))
  
}
colnames(summaryData) <- c("g1", "Dg", "Dp")

# plot         ----
par(mfrow=c(1,2))
plot(x=summaryData[,"g1"], y=summaryData[,"Dg"], ylim=c(0,4), type="o", lwd=2, xlab="plasticity, g1", main="genetic divergence, Dg", ylab="Dg", las=1)
abline(h=2, col="gray", lty=3)
plot(x=summaryData[,"g1"], y=summaryData[,"Dp"], ylim=c(0,4), type="o", lwd=2, xlab="plasticity, g1", main="phenotypic divergence, Dp", ylab="Dp", las=1)
abline(h=2, col="gray", lty=3)
par(mfrow=c(1,1))

#              ----




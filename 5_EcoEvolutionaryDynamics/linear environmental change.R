rm(list = ls(all.names = TRUE))
setwd("~/Desktop/Winterschool-Finlande/NEMOAGE-schoolversion/Evolution")

# check burn-in: genetic variance must have reached equilibrium

b=read.table("evolution-burnin_bygen.txt",header=T)
colnames(b)
genmax=max(b$generation)
par(mfrow=c(2,2))
plot(b$adlt.q1.Va,xlab="time",ylab="additive genetic variance",type="l")
eq_popsize=b[b$generation==genmax,]$pop.tot
plot(b$pop.tot,ylab="population size",xlab="time",type="l",ylim=c(0,1.1*eq_popsize))
# notice that fecundity does not include the loads. 
# if variance/selection too strong, 
# average fecundity is too low to permit population growth

# build change in optimum with time
d=read.table("evolution-envchange_bygen.txt",header=T)
colnames(d)
begin_change=250
k=0.05 # linear speed of environmental change
par(mfrow=c(2,2))
#trait value
plot(d$generation,d$off.q1,ylab="trait value",xlab="time",type="l")
abline(a=-(begin_change*k),b=k,col="red")
#additive genetic variance
plot(d$generation,d$off.q1.Va,ylab="additive genetic variance",xlab="time",type="l")
# notice that selection depletes genetic variance
#population size
plot(d$generation,d$pop.tot,ylab="population size",xlab="time",type="l")
abline(h=eq_popsize,col="blue") # problem for students: why does the population size is larger
# after environmental change than before ?

## Problems :
# change the pace of life history : maturity is reached early
# include clonal reproduction
# compare different adult survival rates




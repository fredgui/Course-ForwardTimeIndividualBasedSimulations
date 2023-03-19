# analysis of nemoage output for killerwhale demography

# load simulation data for killerwhale demography

setwd("~/Desktop/Winterschool-Finlande/NEMOAGE-schoolversion/Demography/")
dd=read.table("killerwhale.txt",header=T)

# expectation

m=matrix(0,nrow=4,ncol=4)
m[1,]=c(0, 0.0043, 0.1132, 0)
m[2,]=c(0.9775, 0.9111, 0, 0)
m[3,]=c(0, 0.0736, 0.9534, 0)
m[4,]=c(0, 0, 0.0452, 0.9804)

eigen(m)

lambda=max(eigen(m)$values)

w=eigen(m)$vectors[,1]/sum(eigen(m)$vectors[,1])

# total population size #

nrep=5
xmax=300
clist=c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e")

nini=20

d=dd[dd$generation<(xmax+1),]
plot(d[d$replicate==1,"generation"],log10(d[d$replicate==1,"pop.tot"]),
     col=clist[1],type="l",ylab="population size",xlab="time",xlim=c(0,xmax))
for(i in 2:nrep) lines(d[d$replicate==i,"generation"],log10(d[d$replicate==i,"pop.tot"]),col=clist[i])
abline(a=log10(nini),b=log10(lambda),col="red",lwd=2)

# to observe: population size initially fluctuates around the expected. 
# ultimately reach asymptotic growth rate, but population size is shifted from expectation
# due to initial conditions. 
# Demographic stochasticity at small population size

# age structure # 

clist2=c("#feebe2", "#fbb4b9", "#f768a1", "#c51b8a" ,"#7a0177")

dd1=dd[dd$replicate==1,]
plot(dd1$generation,log10(dd1$a1.tot),type="l",ylab="log number of individuals",xlab="time", col=clist2[1])
lines(dd1$generation,log10(dd1$a2.tot),col=clist2[2])
lines(dd1$generation,log10(dd1$a3.tot),col=clist2[3])
lines(dd1$generation,log10(dd1$a4.tot),col=clist2[4])

# to observe: the number of individuals in each age class asymptotically grows at rate lambda

dd1=dd[dd$replicate==1,]
plot(dd1$generation,dd1$a1.tot/(dd1$a1.tot+dd1$a2.tot+dd1$a3.tot+dd1$a4.tot),type="l",ylab="log number of individuals",xlab="time", col=clist2[1],ylim=c(0,1))
abline(h=w[1],lty=3)
lines(dd1$generation,dd1$a2.tot/(dd1$a1.tot+dd1$a2.tot+dd1$a3.tot+dd1$a4.tot),col=clist2[2])
abline(h=w[2],lty=2)
lines(dd1$generation,dd1$a3.tot/(dd1$a1.tot+dd1$a2.tot+dd1$a3.tot+dd1$a4.tot),col=clist2[3])
abline(h=w[3],lty=4)
lines(dd1$generation,dd1$a4.tot/(dd1$a1.tot+dd1$a2.tot+dd1$a3.tot+dd1$a4.tot),col=clist2[4])
abline(h=w[4],lty=1)

# to observe : the frequency of each age class asymptotically meet expectation

##### population regulation ######

# offspring survival rate is dependent on adult density. BH: N0={0,100,0}, k=0.001. RK: N0={0,1000,0}, k=0.0001

# get simulation data for Beverton-Holt regulation function

flist=c("0050","0500","1800")
flistn=c(50,500,1800)
for(i in 1:length(flist)){
  assign(paste("dbh",flistn[i],sep=""),read.table(paste("density-dependent-regulation-BH-f",flist[i],"_bygen.txt",sep=""),header=T))
}

ymax=max(c(dbh50$pop.tot,dbh500$pop.tot,dbh1800$pop.tot))

plot(dbh50$generation,dbh50$pop.tot,type="l",ylim=c(0,ymax),ylab="number of individuals",xlab="time")
lines(dbh500$generation,dbh500$pop.tot,col="blue")
lines(dbh1800$generation,dbh1800$pop.tot,col="dark green")

# get simulation data for Ricker regulation function 

flist=c("0050","0500","1800")
flistn=c(50,500,1800)
for(i in 1:length(flist)){
  assign(paste("drk",flistn[i],sep=""),read.table(paste("density-dependent-regulation-RK-f",flist[i],"_bygen.txt",sep=""),header=T))
}

ymax=max(c(drk50$pop.tot,drk500$pop.tot,drk1800$pop.tot))

plot(drk50$generation,drk50$pop.tot,type="l",ylim=c(0,ymax),ylab="number of individuals",xlab="time")
lines(drk500$generation,drk500$pop.tot,col="blue")
lines(drk1800$generation,drk1800$pop.tot,col="dark green")

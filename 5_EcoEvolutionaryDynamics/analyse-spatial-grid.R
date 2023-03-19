#### analyzing spatial eco-evolutionary dynamics ####
## clean repertory and load libraries ####

rm(list=ls(all=TRUE)) 
library(raster)
library(st)

## get simulation results ####

setwd("~/Desktop/Winterschool-Finlande/NEMOAGE-schoolversion/Grid")

# generation pooled # TO DO by replicate 

# burnin
d=read.table("grid-envchange-burnin_bygen.txt",header=T)
# projection
p=read.table("grid-envchange_bygen.txt",header=T)

# define grid size (remember to get information from SDM simulation)=
nl=10
ncol=30

## spatial demography ####

MyDemoGrid=function(dat,gen){
  
  L=dat$generation==gen
  
  mat=rep(0,(nl*ncol))
  
  for(i in 1:(nl*ncol)){
    
    a = paste("a1.fem.p",i,sep="")
    b = paste("a2.fem.p",i,sep="")
    c = paste("a3.fem.p",i,sep="")
    
    mat[i]=dat[L,a]+dat[L,b]+dat[L,c] # matrix is rotated
    
  }
  
  dim(mat)=c(nl,ncol)
  
  return(mat)
}

par(mfrow=c(2,2))

matg500=MyDemoGrid(d,500)
plot(raster(matg1),main="burn-in")

matg1proj=MyDemoGrid(p,1)
plot(raster(matg1proj),main="t0")

matg10proj=MyDemoGrid(p,10)
plot(raster(matg10proj),main="t10")

## evolutionary stat ####
par(mfrow=c(1,2))
plot(p$generation,p$adlt.q1, main="Environment variable 1",xlab="Time")
plot(p$generation,p$adlt.q2, main="Environment variable 2",xlab="Time")

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
dmut=read.table("grid-envchange-burnin-mut_bygen.txt",header=T)

# projection
p=read.table("grid-envchange_bygen.txt",header=T)
pmut=read.table("grid-envchange-mut_bygen.txt",header=T)


# define grid size (remember to get information from SDM simulation)
nl=10
ncol=30

## global demography ####

# no mutation
par(mfrow=c(1,2))
plot(d$pop.tot,type="l",xlab="Time",ylim=c(mean(d$pop.tot)*0.5,mean(d$pop.tot)*1.5))
plot(p[,"pop.tot"],type="l",xlab="Time",main="Prediction",col="red",ylim=c(mean(d$pop.tot)*0.5,mean(d$pop.tot)*1.5))

# mutation
par(mfrow=c(1,2))
plot(dmut$pop.tot,type="l",xlab="Time",ylim=c(mean(dmut$pop.tot)*0.5,mean(dmut$pop.tot)*1.5))
plot(pmut[,"pop.tot"],type="l",xlab="Time",main="Prediction",col="red",ylim=c(mean(dmut$pop.tot)*0.5,mean(dmut$pop.tot)*1.5))

## evolutionary dynamics

#no mutation
plot(d$adlt.q1,type="l",xlab="Time")
plot(p$adlt.q1,type="l",xlab="Time",col="red")

#mutation
plot(dmut$adlt.q1,type="l",xlab="Time")
plot(pmut$adlt.q1,type="l",xlab="Time",col="red")

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
plot(raster(matg500),main="burn-in")

matg1proj=MyDemoGrid(p,100)
plot(raster(matg1proj),main="before climate change")

matg1proj=MyDemoGrid(p,110)
plot(raster(matg1proj),main="10 years climate change")

matg10proj=MyDemoGrid(p,115)
plot(raster(matg10proj),main="15 years climate change")

##

par(mfrow=c(2,2))

matg500=MyDemoGrid(dmut,500)
plot(raster(matg500),main="burn-in")

matg1proj=MyDemoGrid(pmut,100)
plot(raster(matg1proj),main="before climate change")

matg1proj=MyDemoGrid(pmut,110)
plot(raster(matg1proj),main="10 years climate change")

matg10proj=MyDemoGrid(pmut,115)
plot(raster(matg10proj),main="15 years climate change")

## animation

PlotMyDemoGrid=function(dat,g){
  m=MyDemoGrid(dat,g)
  plot(raster(m),main=g)
}

PlotMyDemoGrid(dmut,100)

library("animation")
saveHTML({timeseq=c(1,seq(0,120,1)[-1]);
  for(i in 1:length(timeseq)) PlotMyDemoGrid(p,timeseq[i])},
  interval=0.1, img.name = "climchange", imgdir="gif_files",
  autoplay = FALSE, 
  autobrowse = FALSE,
  ani.height = 600,
  ani.width = 600,
  htmlfile = "rangeclimchg.html",
  title="range during climchange", 
  description="range during climchange")
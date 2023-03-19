## Clean repertories, load data and libraries ####

rm(list=ls(all=TRUE)) 

setwd("~/Desktop/Winterschool-Finlande/simulate-SDM")

library(raster, quietly = T)
library(virtualspecies)
library(rworldmap,quiet=T)
library(terra)
library(st)

# load current climate

worldclim <- getData("worldclim", var = "bio", res = 2.5)

# load futur climate (here rcp 4.5 for 2061-2080 ; see help for getData for details)

climateproj<-getData(name = 'CMIP5', var = 'bio', res = 2.5,
                            rcp = 45, model = 'IP', year = 70)

## Function to write matrices in nemo ####

write.matrix.nemo = function(mat, outfile) 
{
  rows = dim(mat)[1]
  cols = dim(mat)[2]
  
  cat("{",file=outfile)
  
  for(i in 1:rows) {
    cat("{",file=outfile, append=TRUE)
    cat(mat[i,],sep=",", file=outfile, append=TRUE)
    cat("}\n",file=outfile, append=TRUE)
  }
  
  cat("}\n",file=outfile, append=TRUE)
}



## Generate random species ####

# We use the package virtual species to generate random SDM

# get the spatial coordinate for a given area - here Finland 
finland <- getData("GADM", country = "FIN", level = 0)

# restrain data to the chosen area
var1=crop(worldclim[["bio1"]],finland)
var2=crop(worldclim[["bio12"]],finland)
my.stack=stack(var1,var2)

# simulate SDM (see help for virtual species ; http://borisleroy.com/files/virtualspecies-tutorial.html)
random.sp <- generateRandomSp(my.stack,relations="gaussian",realistic.sp = TRUE)

random.sp

# change projection system to have resolution in meters (more straightforward to think on dispersal rate in distance units)
suitmap=projectRaster(random.sp$suitab.raster,crs="+init=epsg:3035")
presmap=projectRaster(random.sp$pa.raster,crs="+init=epsg:3035")

# plot
par(mfrow=c(1,1))
plot(suitmap)

# and select a subgrid directly by drawing it on plot
myext=drawExtent()

# restrain the remaining of the analysis to the subgrid above defined
suitmap=crop(suitmap,myext)
presmap=crop(presmap,myext)

# visualize
par(mfrow=c(1,2))
plot(suitmap)
plot(presmap)

suitmap


## Generate matrices for i) suitability ii) carrying-capacity ####

datsuit=as.matrix(suitmap*100) # suitability is a probability
write.matrix.nemo(datsuit,"current_suitability.txt")

#threshold=0.25 # above this suitability threshold we assume the species to be present
#datpres=(datsuit>threshold)+0
datpres=as.matrix(presmap)
write.matrix.nemo(datpres*1000,"current_occurrence.txt") # 1000 as maximum carrying capacity - can be modified

## Generate matrix for environment ####

e1="bio1"
e2="bio12"

envmatrix=function(environmental.variable,area1,area2){
  
  # function to generate environment raster for the subgrid chosen
  
  d=worldclim[[environmental.variable]]
  
  d=crop(d,area1)
  
  d=projectRaster(d,crs="+init=epsg:3035")
  
  d=crop(d,area2)
  
  return(d)
  
}

de1 <- extract(envmatrix(e1,finland,myext),myext, cellnumbers=TRUE)
de2<-extract(envmatrix(e2,finland,myext),myext, cellnumbers=TRUE)

# environmental values are paired for nemo input - see nemoage manual
de=cbind(de1[,2],de2[,2])
write.matrix.nemo(de,"current_env.txt")

## Generate dispersal kernel ####

# we assume that the dispersal kernel is gaussian (other options possible)
# matrix with spatialized patch identity
# see nemoge manual to understand the outputs

BuildDispersalKernel=function(my_rast){

l=dim(my_rast)[1]

w=dim(my_rast)[2]

disp1=c(1:(l*w))

dim(disp1)=c(l,w); 

disp1=t(disp1);

#matrix miminum for interconnected patches: number of patches * size disp kernel*4

patch_connected={} 

probaD={}

final_matrix_pD=matrix(0,nrow=(l*w),ncol=l*w)

final_matrix_p_connected=matrix(0,nrow=(l*w),ncol=l*w)

for(i in 1:(l*w)){
  
  a=which(disp1==i,arr.ind=T) 
  
  for(j in 1:(l*w)){  #l*w  1:(l*w)
    
    b=which(disp1==j,arr.ind=T)
    
    c=(as.numeric((a[1]-b[1]))^2+as.numeric((b[2]-a[2]))^2)^(1/2)
    
    patch_connected=c(patch_connected,j);
  
    probaD=c(probaD,dnorm(c,mean=0,sd=1));
    
    #note that other distribution could be chosen for kernel
  
  }
  
  
  #order highest proba first
  
  ord=rbind(patch_connected,probaD)
  
  ord=ord[,order(ord[2,],decreasing=TRUE)]
  
  #final matrices
  
  final_matrix_pD[i,1:length(probaD)]=ord[2,]
  
  final_matrix_p_connected[i,1:length(patch_connected)]=ord[1,]
  
  patch_connected={} 
  
  probaD={}
  
}

#remove very small values in pdisp

for(i in 1:dim(final_matrix_pD)[1]){for(j in 1:dim(final_matrix_pD)[2]){if(final_matrix_pD[i,j]<10^(-6)){final_matrix_pD[i,j]=0;final_matrix_p_connected[i,j]=0}}}

#rescale to 1 final matrix p dispersion

for(i in 1:dim(final_matrix_pD)[1]){
  
  s=1/sum(final_matrix_pD[i,]); 
  
  final_matrix_pD[i,]=final_matrix_pD[i,]*s
  
}


return(list(final_matrix_pD,final_matrix_p_connected))

}

dispKernel=BuildDispersalKernel(suitmap)

write.matrix.nemo(dispKernel[[1]],"p_Disp.txt")

write.matrix.nemo(dispKernel[[2]],"p_Connected.txt")

## RCP scenarios and changing rate ####

# get future environment ; notice that names for environmental values are different from those of current climate
f1=paste("ip45bi70",1,sep="") # for bio1
f2=paste("ip45bi70",12,sep="") # for bio12

Projenvmatrix=function(environmental.variable,area1,area2){
  
  d=climateproj[[environmental.variable]]
  
  d=crop(d,area1)
  
  d=projectRaster(d,crs="+init=epsg:3035")
  
  d=crop(d,area2)
  
  return(d)
  
}

# projected temperature are given in 10*C

f70b1=extract(Projenvmatrix(f1,finland,myext)/10,myext,cellnumbers=TRUE)
f70b2=extract(Projenvmatrix(f2,finland,myext),myext,cellnumbers=TRUE)

# from here we need a rate of environmental change by cell
# we are currently in 2023 and we got the annual temperature for the decade 2061-2080

# the maximum rate of change is then (Tpt[2061]-Tpt[2023])/38

rateb1=(f70b1[,2]-de[,1])/(2061-2023)
rateb2=(f70b2[,2]-de[,2])/(2061-2023)

futurenv=cbind(rateb1,rateb2)
write.matrix.nemo(futurenv,"rate_chgt_env.txt")

## Define the realized niche ####

niche1=c(((datpres>0)+0)*as.matrix(envmatrix(e1,finland,myext)))
niche1=niche1[niche1!=0]
mean(niche1)
sd(niche1)

niche2=c(((datpres>0)+0)*as.matrix(envmatrix(e2,finland,myext)))
niche2=niche2[niche2!=0]
mean(niche2)
sd(niche2)

## Save information ####
# First, it is wise to save the R-environment if later analyses are required
save.image(file='mySpeciesEnvironment.RData')
# In addition, we may want to note of a number of information
notefile=file("note_myspecies.txt")
  # size of the grid
sink(notefile)

print(paste("The grid has ",dim(suitmap)[1]," lines and ",dim(suitmap)[2]," colomuns"))
  # properties of the niche

print(paste("the mean value of first environmental variable in occupied sites is", mean(niche1), "with sd ",sd(niche1)))
print(paste("the mean value of second environmental variable in occupied sites is", mean(niche2), "with sd ",sd(niche2)))

print(random.sp)

sink() #close external connexion




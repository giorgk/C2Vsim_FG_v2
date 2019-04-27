library(h5)
library(pracma)

# Set as current directory the directory where this script is run
setwd('f:/UCDAVIS/C2VSIM_FG_OR/C2Vsim_FG_v2/Rwrkspc/')
setwd(getwd())
load("FaceData.RData")

# Set where the C2Vsim model is
c2vsim_path <- "f:/UCDAVIS/C2VSIM_FG_OR/C2Vsim_FG_v2/C2VSimFG-BETA_PublicRelease/"
results_path <- paste0(c2vsim_path, "Results/")
preproc_path <- paste0(c2vsim_path, "Preprocessor/")

# Read Budget file
GW_BDGinfo <- h5file(name =  paste0(results_path, "C2VSimFG_GW_ZBudget.hdf"), mode = "r")
hfileGroups <- list.groups(GW_BDGinfo)
hfileDataSets <- list.datasets(GW_BDGinfo)

# Deep percolation
# Assume that the units are AC.FT /month
# 1 ACFT -> 1233.48 m^3
cnvrt <- 1233.48
DeepPerc <- GW_BDGinfo[hfileDataSets[24]]
Vflow <- array(dim = c(dim(DeepPerc)[2], dim(DeepPerc)[1], 4))
for (i in 1:dim(DeepPerc)[1]) {
  Vflow[,i,1] <- DeepPerc[i,]*cnvrt
}

# The vertical flows are written per node.
# Read the vertical flows
ids <- c(46, 73, 100)
vertflowNodes <- array(dim = c(dim(XY)[1], 504, length(ids)))
for (i in 1:length(ids)) {
  Vertflow <- GW_BDGinfo[hfileDataSets[ids[i]]]
  for (j in 1:dim(Vertflow)[1]) {
    vertflowNodes[,j,i] <- Vertflow[j,]
  }
}

# First we have to find out how many element share each node
NsharedElem <- vector(mode = "integer", length = dim(XY)[1])
for(i in 1:dim(XY)[1]){
  elemlist = which(MSH[c(-1,-6)] == i, arr.ind = TRUE)
  NsharedElem[i] <- dim(elemlist)[1]
}

# For each element we will add the vertical flows of the nodes divided by the number of elements each node is connected
for (i in 1:dim(MSH)[1]) {
  velemflow <- matrix(data = 0, nrow = 504, ncol = 3)
  for(j in 2:5){
    if (MSH[i,j] == 0)
      break
    velemflow <- velemflow + vertflowNodes[MSH[i,j],,]/as.numeric(NsharedElem[MSH[i,j]])
  }
  
  for(j in 1:3){
    Vflow[i,,j+1] <- velemflow[,j]*cnvrt
  }
}

## Read the faceflows
hflow <- array(dim = c(dim(FcLm)[1], 504, 4))
ids <- c(28, 55, 82, 109)
for (i in 1:length(ids)) {
  tempfflow <- GW_BDGinfo[hfileDataSets[ids[i]]]
  for (j in 1:dim(tempfflow)[1]) {
    hflow[,j,i] <- tempfflow[j,]*cnvrt
  }
}

# Divide the vertical flows with the element areas
for (i in 1:dim(Vflow)[2]) {
  for (j in 1:dim(Vflow)[3]) {
    Vflow[,i,j] <- Vflow[,i,j]/elemArea
  }
}

# Divide the horizontal flows with the face areas
for (i in 1:dim(hflow)[2]) {
  for (j in 1:dim(hflow)[3]) {
    hflow[,i,j] <- hflow[,i,j]/faceArea[,j]
  }
}

## Convert dataframes to matrices
XYmat <- matrix(nrow = dim(XY)[1], ncol = 2)
XYmat[,1] <- XY$X
XYmat[,2] <- XY$Y

MSHmat <- matrix(nrow = dim(MSH)[1], ncol = 5)
MSHmat[,1] <- MSH$nd1
MSHmat[,2] <- MSH$nd2
MSHmat[,3] <- MSH$nd3
MSHmat[,4] <- MSH$nd4
MSHmat[,5] <- MSH$S

STRATmat <- matrix(nrow = dim(strat)[1], ncol = 5)
STRATmat[,1] <- strat$GSE
STRATmat[,2] <- strat$L1
STRATmat[,3] <- strat$L2
STRATmat[,4] <- strat$L3
STRATmat[,5] <- strat$L4



hf <- h5file(name = "PartTrackData.h5", mode = "w")
hf["flowdata/VFLOW"] <- Vflow
hf["flowdata/HFLOW"] <- hflow
hf["geodata/XY"] <- XYmat
hf["geodata/MSH"] <- MSHmat
hf["geodata/STRAT"] <- STRATmat
hf["geodata/FI"] <- faceIndex
hf["geodata/FCEL"] <- FcLm
h5close(hf)

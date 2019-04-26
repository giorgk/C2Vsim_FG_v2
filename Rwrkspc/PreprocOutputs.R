# Set as current directory the directory where this script is run
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
DeepPerc <- GW_BDGinfo[hfileDataSets[24]]
Vflow <- array(dim = c(dim(DeepPerc)[2], dim(DeepPerc)[1], 4))
for (i in 1:dim(DeepPerc)[1]) {
  Vflow[,i,1] <- DeepPerc[i,]
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
    Vflow[i,,j] <- velemflow[,j]
  }
}

## Read the faceflows
hflow <- array(dim = c(dim(FcLm)[1], 504, 4))
ids <- c(28, 55, 82, 109)
for (i in 1:length(ids)) {
  tempfflow <- GW_BDGinfo[hfileDataSets[ids[i]]]
  for (j in 1:dim(tempfflow)[1]) {
    hflow[,j,i] <- tempfflow[j,]
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

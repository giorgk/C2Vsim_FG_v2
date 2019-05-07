# load the required library
library(h5)
library(pracma)

# Define the paths
c2vsim_path <- "f:/UCDAVIS/C2VSIM_FG_OR/C2Vsim_FG_v2/C2VSimFG-BETA_PublicRelease/"
results_path <- paste0(c2vsim_path, "Results/")
preproc_path <- paste0(c2vsim_path, "Preprocessor/")

### Read Ascii files
## Read coordinate file
# Set the number of header lines and the number of XY nodes
NlinesSkip <- 90
nNodes <- 30179
XY <- read.table(file = paste0(preproc_path, "C2VSimFG_Nodes.dat"),
                 header = FALSE, sep = "", skip = NlinesSkip, nrows = nNodes,
                 quote = "",fill = TRUE,
                 col.names = c("ID", "X", "Y"))
## plot the nodes
# plot(XY$X,XY$Y)

## Read stratigraphy file
# GSE is the Ground Surface Elevation
# A columns contain the thicknes of aquicludes
# L columns contain the thicknes of main layers 
NlinesSkip <- 105
strat <- read.table(file = paste0(preproc_path, "C2VSimFG_Stratigraphy.dat"),
                    header = FALSE, sep = "", skip = NlinesSkip, nrows = nNodes,
                    quote = "",fill = TRUE,
                    col.names = c("ID", "GSE", "A1", "L1", "A2", "L2", "A3", "L3", "A4", "L4"))

# Convert feet to meter
strat[-1] <- strat[-1] * 0.3048
# For each of the 4 layers add the thickness of aquiclude to the layer
for(i in 1:4){
  strat[[ paste0("L", as.character(i))]] <- strat[[ paste0("L", as.character(i))]] + strat[[ paste0("A", as.character(i))]]
}
#Delete the aquicludes and Keep only the layer thicknesses
strat$A1 <- NULL
strat$A2 <- NULL
strat$A3 <- NULL
strat$A4 <- NULL
# Subtract the thickness from the surface to obtain the elevation of the bottom of the first layer
strat$L1 = strat$GSE - strat$L1 
strat$L2 = strat$L1 - strat$L2 
strat$L3 = strat$L2 - strat$L3
strat$L4 = strat$L3 - strat$L4


## Read Mesh file
NlinesSkip <- 142
nNodes <- 32537
MSH <- read.table(file = paste0(preproc_path, "C2VSimFG_Elements.dat"),
                  header = FALSE, sep = "", skip = NlinesSkip, nrows = nNodes,
                  quote = "",fill = TRUE,
                  col.names = c("ID", "nd1", "nd2", "nd3", "nd4", "S"))

elemArea <- vector(mode = "numeric", length = length(MSH[[1]]))

# Calculate the element areas
print("Calculate element areas...")
for (i in 1:length(MSH[[1]])){
  x <- XY$X[as.integer(c(MSH[i,2:5]))]
  y <- XY$Y[as.integer(c(MSH[i,2:5]))]
  elemArea[i] <- polyarea(x,y)
}

# Calculate the face areas
# initialize a list of matrices to hold the face areas for each layer

faceArea <- array( dim = c(length(MSH[[1]]), 4, 4 )) # elements, face per element, layers

### Read the HDF5 budget file
GW_BDGinfo <- h5file(name =  paste0(results_path, "C2VSimFG_GW_ZBudget.hdf"), mode = "r")
hfileGroups <- list.groups(GW_BDGinfo)
hfileDataSets <- list.datasets(GW_BDGinfo)
faceElem <- GW_BDGinfo[hfileDataSets[17]]

#Gface <- make_empty_graph(n = dim(MSH)[1], directed = FALSE)
#for (i in 1:dim(faceElem)[1]){
#  add_edges(Gface, faceElem[i])
#}

faceArea = matrix(data = 0, nrow = dim(faceElem)[1], ncol = 4)
faceIndex = matrix(data = 0L, nrow = dim(MSH)[1], ncol = 4)

## Calculate the face areas and indices for the inner faces of the Mesh
print("Calculate inner element face areas and indices...")
for (i in 1:dim(faceElem)[1]){
  ela <- faceElem[i][1]
  if (ela != 0 ){
    na <- 4
    if (MSH[ela,5] == 0)
      na <- 3
  }
  
  elb <- faceElem[i][2]
  if (elb != 0 ){
    nb <- 4
    if (MSH[elb,5] == 0)
      nb <- 3
  }
  
  if (ela != 0 && elb != 0){
    # find the vertex indices for the common face between the two elements
    breatThis <- FALSE
    for (ii in 1:na){
      if (faceIndex[ela,ii] != 0)
        next
      ia <- MSH[ela,ii+1]
      iii <- ii + 1
      if (ii == na)
        iii = 1
      ib <- MSH[ela, iii+1]
      for (jj in 1:nb){
        if (faceIndex[elb,jj] != 0)
          next
        ja <- MSH[elb,jj+1]
        jjj <- jj +1
        if (jj == nb)
          jjj <- 1
        jb <- MSH[elb,jjj+1]
        
        if ((ia == ja & ib == jb) | (ia == jb & ib == ja)){
          faceIndex[ela,ii] = i
          faceIndex[elb,jj] = -i
          # calculate the area of the face for each layer
          L <- sqrt((XY$X[ia] - XY$X[ib])^2 + (XY$Y[ia] - XY$Y[ib])^2)

          xv <- c(0, L, L, 0)
          for (k in 1:4){
            yv <- c(strat[ia,k+2], strat[ib,k+2], strat[ib,k+1], strat[ia,k+1])
            faceArea[i,k] <- polyarea(xv,yv)
          }
          breatThis <- TRUE
          break
        }
      }
      if (breatThis)
        break
    }
  } 
}

FcLm = matrix(nrow = dim(faceElem)[1], ncol = dim(faceElem)[2])
for (i in 1:dim(FcLm)[1]){
  FcLm[i,1] <- faceElem[i][1]
  FcLm[i,2] <- faceElem[i][2]
}

print("Calculate outer element face areas and indices...")
## Calculate the face areas and indices for the boundary faces of the Mesh
for (i in 1:dim(faceIndex)[1]) {
  # find how many faces do not have index
  zrID <- which(faceIndex[i,] == 0)
  if (MSH[i, 5] == 0){
    nonValid <- which(zrID == 4)
    zrID <- zrID[-nonValid]
  }
  if (!isempty(zrID)){
    # Find the faces that are associated with the element
    elemFCid <- which(FcLm == i,arr.ind = TRUE)
    
    # remove the faces that they are connected with another element
    dlt <- vector(mode = "integer")
    for (j in 1:dim(elemFCid)[1]){
      oth_el <- 2
      if (elemFCid[j,2] == 2)
        oth_el <- 1
      
      if (FcLm[elemFCid[j,1], oth_el] == 0)
        dlt <- append(dlt, j)
    }
    elemFCid <- elemFCid[dlt,]
    
    if (is.null(dim(elemFCid)[1])){
      if (!(length(elemFCid) == 2 & length(zrID) == 1)){
        print(paste0("Mismatch between faces for element: ",i) )
      }
      else{
        if (elemFCid[2] == 1){
          faceIndex[i,zrID[1]] = elemFCid[1]
        }else{
          faceIndex[i,zrID[1]] = -elemFCid[1]
        }
        ia = MSH[i,zrID[1]+1]
        ib = MSH[i,zrID[1]+2]
        if (ib == 0 | zrID[1]+1 == 5){
          ib = MSH[i,2]
        }
        L <- sqrt((XY$X[ia] - XY$X[ib])^2 + (XY$Y[ia] - XY$Y[ib])^2)
        xv <- c(0, L, L, 0)
        for (k in 1:4){
          yv <- c(strat[ia,k+2], strat[ib,k+2], strat[ib,k+1], strat[ia,k+1])
          faceArea[elemFCid[1], k] <- polyarea(xv,yv)
        }
      }
    } else{
      if (dim(elemFCid)[1] != length(zrID)){
        print(paste0("Mismatch between faces for element: ",i) )
      }else{
        for (j in 1:length(zrID)){
          if (elemFCid[j,2] == 1){
            faceIndex[i,zrID[j]] = elemFCid[j,1]
          }else{
            faceIndex[i,zrID[j]] = -elemFCid[j,1]
          }
          ia = MSH[i,zrID[j]+1]
          ib = MSH[i,zrID[j]+2]
          if (ib == 0 | zrID[j]+1 == 5){
            ib = MSH[i,2]
          }
          L <- sqrt((XY$X[ia] - XY$X[ib])^2 + (XY$Y[ia] - XY$Y[ib])^2)
          xv <- c(0, L, L, 0)
          for (k in 1:4){
            yv <- c(strat[ia,k+2], strat[ib,k+2], strat[ib,k+1], strat[ia,k+1])
            faceArea[elemFCid[j,1], k] <- polyarea(xv,yv)
          }
        }
      }
    }
  }
}

## Save data in R format
# save(XY, MSH, strat, faceIndex, faceArea, FcLm, elemArea, file = "FaceData.RData")




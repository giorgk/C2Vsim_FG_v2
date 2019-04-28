library(h5)
library(pracma)
library(xlsx)

# This script is used to extract and the flow field from any month
# Before useing this script you have to run the PreprocInputs and PreprocOutputs scripts to generate the
# PartTrackData.h5 file. THis file is also needed for the 

# Simulation Time series
simTime <- seq.Date(from = as.Date("1973/10/30"),to = as.Date("2015/09/30"),by = "month")

# Select the times to export the flow field
#------------------------------------------------------------------------------------------------------------------------
# Modify the following as needed.
# At the moment extracts all Aprils after the year 2000
selectTimes <- which(as.numeric(format(simTime,"%Y")) >= 2000 & as.numeric(format(simTime,"%m")) == 4, arr.ind = TRUE)
# -----------------------------------------------------------------------------------------------------------------------

# Load the flow field data from the PartTrackData.h5
hpart <- h5file(name =  "PartTrackData.h5", mode = "r")
hpartGroups <- list.groups(hpart)
hpartDataSets <- list.datasets(hpart)

XY <- hpart["/geodata/XY"]
MSH <- hpart["/geodata/MSH"]

# Calculate element barycenters
cc <- matrix(data = 0, nrow = dim(MSH)[1], ncol = 2)
quad_el <- which(MSH[,4] != 0, arr.ind = FALSE)
tri_el <- which(MSH[,4] == 0, arr.ind = FALSE)

x1 <- XY[MSH[quad_el,1],1]; y1 <- XY[MSH[quad_el,1],2]
x2 <- XY[MSH[quad_el,2],1]; y2 <- XY[MSH[quad_el,2],2]
x3 <- XY[MSH[quad_el,3],1]; y3 <- XY[MSH[quad_el,3],2]
x4 <- XY[MSH[quad_el,4],1]; y4 <- XY[MSH[quad_el,4],2]
cc[quad_el,1] <- (x1 + x2 + x3 +x4)/4
cc[quad_el,2] <- (y1 + y2 + y3 +y4)/4

x1 <- XY[MSH[tri_el,1],1]; y1 <- XY[MSH[tri_el,1],2]
x2 <- XY[MSH[tri_el,2],1]; y2 <- XY[MSH[tri_el,2],2]
x3 <- XY[MSH[tri_el,3],1]; y3 <- XY[MSH[tri_el,3],2]
cc[tri_el,1] <- (x1 + x2 + x3)/3
cc[tri_el,2] <- (y1 + y2 + y3)/3



# Calculate face centers and normals
FCEL <- hpart["/geodata/FCEL"]
nrmls <- matrix(data = 0, nrow = dim(FCEL)[1], ncol = 4)
inner_fc <- which(FCEL[,1] !=0 & FCEL[,2] !=0)
outer_fc <- which(FCEL[,1] ==0 | FCEL[,2] ==0)
nrmls[inner_fc,1] <- ( cc[FCEL[inner_fc,1],1] + cc[FCEL[inner_fc,2],1] )/2
nrmls[inner_fc,2] <- ( cc[FCEL[inner_fc,1],2] + cc[FCEL[inner_fc,2],2] )/2

dx <- cc[FCEL[inner_fc,2],1] - cc[FCEL[inner_fc,1],1]
dy <- cc[FCEL[inner_fc,2],2] - cc[FCEL[inner_fc,1],2]
ln <- sqrt(dx^2 + dy^2)
nrmls[inner_fc,3] <- dx/ln 
nrmls[inner_fc,4] <- dy/ln

hFI <- hpart["/geodata/FI"]
FI = matrix(nrow = dim(hFI)[1], ncol = dim(hFI)[2])
for (i in 1:dim(hFI)[2]) {
  FI[,i] = hFI[,i]
}


# For the outer faces use a loop
for (i in 1:length(outer_fc)) {
  # Find the element and the outer face index 
  iel <- which(abs(FI) == outer_fc[i], arr.ind = TRUE )
  if (dim(iel)[1] > 1){
    print("More than one elements are found")
  }
    
  a <- XY[MSH[iel[1],iel[2]][1],]
  if (MSH[iel[1],iel[2]+1][1] == 0 || iel[2] == 4){
    b <- XY[MSH[iel[1], 1][1],]
  }else{
    b <- XY[MSH[iel[1],iel[2]][1],]
  }
  c <- (a+b)/2
  nrmls[outer_fc[i],1:2] <- c
  
  ccel <- cc[iel[1],]
  
  # if the index is negative then the flow comes from outer to inside of element
  
  if (FI[iel] < 0){
    nn <- ccel - c
  }else{
    nn <- c - ccel
  }
  
  nn <- nn/sqrt(sum(nn^2))
  nrmls[outer_fc[i],3:4] <- nn
}

vflow <- hpart["flowdata/VFLOW"]
hflow <- hpart["flowdata/HFLOW"]

vf <- data.frame("CX" = cc[,1])
vf["CY"] <- cc[,2]
for (i in 1:length(selectTimes)) {
  vf[as.character(simTime[selectTimes[i]])] <- vflow[,selectTimes[i],1]
}


hf <- data.frame("CX" = nrmls[,1])
hf["CY"] <- nrmls[,2]
hf["NX"] <- nrmls[,3]
hf["NY"] <- nrmls[,4]
for (i in 1:length(selectTimes)) {
  hf[as.character(simTime[selectTimes[i]])] <- hflow[,selectTimes[i],1]
}


write.csv(vf, file = "VertflowData.csv")
write.csv(hf, file = "horflowData.csv")


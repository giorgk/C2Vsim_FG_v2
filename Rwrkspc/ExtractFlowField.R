library(h5)
library(pracma)

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

# For the outer faces use a loop







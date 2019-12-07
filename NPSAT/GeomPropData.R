library("akima")
source("../../C2VsimCG/Rwrkspc/c2vsim_io.R")
source("npsat_functions.R")
# Prepare Input file for bottom and Top -------------------------------------------
# However the top is not going to be used in the simulation.

# Read the node coordinates
XY <- c2vsim.readNodes(filename = "../c2vsimfg_beta2_publicrelease/C2VSimFG_BETA2_PublicRelease/Preprocessor/C2VSimFG_Nodes.dat", ND = 30179, Nskip = 90)

# Read the stratigraphy file
CVstrat <- c2vsim.readStrat("../c2vsimfg_beta2_publicrelease/C2VSimFG_BETA2_PublicRelease/Preprocessor/C2VSimFG_Stratigraphy.dat",
                 nSkip = 105, Nlay = 4, Nnodes = 30179)

TopElev <- CVstrat$ELV * c2vsim.units.ft2m()
BotElev <- TopElev - apply(CVstrat[,c(-1,-2)],1,sum) * c2vsim.units.ft2m()

# Read the expanded Mesh
ExpandedMesh <- npsat.readOBJmesh("ExpandedMesh.obj")

# Interpolate the top and bottom layer on the expanded mesh
TOP <- interpp(x = XY$X, y = XY$Y, z = TopElev, linear = FALSE, extrap = TRUE, duplicate = "mean", 
        xo = ExpandedMesh[[1]][,1], yo = ExpandedMesh[[1]][,2])

BOT <- interpp(x = XY$X, y = XY$Y, z = BotElev, linear = FALSE, extrap = TRUE, duplicate = "mean", 
               xo = ExpandedMesh[[1]][,1], yo = ExpandedMesh[[1]][,2])

npsat.input.WriteScattered(filename = "inputfiles/TopElev.npsat", PDIM = 2, TYPE = "HOR", MODE = "SIMPLE", 
                     DATA = cbind(TOP[[1]], TOP[[2]], TOP[[3]]))

npsat.input.WriteScattered(filename = "inputfiles/BotElev.npsat", PDIM = 2, TYPE = "HOR", MODE = "SIMPLE", 
                     DATA = cbind(BOT[[1]], BOT[[2]], BOT[[3]]))


# Prepare Input file for Hydraulic conductivity ---------------------------
# The hydraulic conductivity in C2Vsim Fine grid Beta version is defined from a parametric grid
GWFile <- "../c2vsimfg_beta2_publicrelease/C2VSimFG_BETA2_PublicRelease/Simulation/Groundwater/C2VSimFG_Groundwater1974.dat"
PMSH <-   read.table(file = GWFile,
                               header = FALSE, sep = "", skip = 354, nrows = 1419,
                               quote = "",fill = TRUE,
                               col.names = c("IE", "ND1", "ND2", "ND3", "ND4"))
tempLines <- readLines(GWFile)
tempLines <- tempLines[1794:7405]
iln <- seq(1, length(tempLines), 4)
PXY <- matrix(data = NA, nrow = 1403, ncol = 2)
PKH <- matrix(data = NA, nrow = 1403, ncol = 4)
PS <- matrix(data = NA, nrow = 1403, ncol = 4)
PN <- matrix(data = NA, nrow = 1403, ncol = 4)
PV <- matrix(data = NA, nrow = 1403, ncol = 4)
PL <- matrix(data = NA, nrow = 1403, ncol = 4)
ind <- 1
for (i in iln) {
  tmp <- strsplit(substr(tempLines[i],1,200), split = "\t")[[1]]
  tmp <- tmp[which(tmp != "")]
  tmp <- as.numeric(tmp)
  PXY[ind,] <- tmp[2:3]
  PKH[ind,1] <- tmp[4]
  PS[ind,1] <- tmp[5]
  PN[ind,1] <- tmp[6]
  PV[ind,1] <- tmp[7]
  PL[ind,1] <- tmp[8]
  for (j in 2:4) {
    tmp <- strsplit(substr(tempLines[i+j-1],1,200), split = "\t")[[1]]
    tmp <- tmp[which(tmp != "")]
    tmp <- as.numeric(tmp)
    PKH[ind, j] <- tmp[1]
    PS[ind, j] <- tmp[2]
    PN[ind, j] <- tmp[3]
    PV[ind, j] <- tmp[4]
    PL[ind, j] <- tmp[5]
  }
  ind <- ind + 1
}

# Write the parametric Mesh into obj file for visual inspection in houdini
npsat.writeMesh2obj(filename = "paraMesh.obj", XY = PXY, MSH = PMSH[,-1])
# Read the expanded mesh
pMSH_exp <- npsat.readOBJmesh("pMeshExpanded.obj")

# Interpolate the Parametric Values of Hydraulic conductivity to the actual elevation nodes
HK_Layer <- matrix(data = NA, nrow = dim(ExpandedMesh[[1]])[1], ncol = 4)
for (i in 1:4) {
  # first interpolate the parametric grid to the actual mesh
   temp <- interpp(x = PXY[,1], y = PXY[,2], z = PKH[,i], linear = FALSE, extrap = TRUE, duplicate = "mean", 
                    xo = XY$X, yo = XY$Y)[[3]]
   
  # Then extrapolate to the expanded grid
  HK_Layer[,i] <- interpp(x = XY$X, y = XY$Y, z = temp, linear = FALSE, extrap = TRUE, duplicate = "mean", 
                          xo = ExpandedMesh[[1]][,1], yo = ExpandedMesh[[1]][,2])[[3]]
}

ElevLay <- matrix(data = NA, nrow = dim(ExpandedMesh[[1]])[1], ncol = dim(HK_Layer)[2]-1)
# Calculate the elevations between the layers
for (i in 1:3) {
  temp <- TopElev - apply(CVstrat[,3:(i*2+2)], 1, sum) * c2vsim.units.ft2m()
  ElevLay[,i] <- interpp(x = XY$X, y = XY$Y, z = temp, linear = FALSE, extrap = TRUE, duplicate = "mean", 
          xo = ExpandedMesh[[1]][,1], yo = ExpandedMesh[[1]][,2])[[3]]
}

tmpHK <- cbind(ExpandedMesh[[1]][,1:2], HK_Layer[,1])
for (i in 1:3) {
  tmpHK <- cbind(tmpHK, ElevLay[,i], HK_Layer[,i+1])
}

npsat.input.WriteScattered(filename = "inputfiles/HK.npsat",PDIM = 2, TYPE = "FULL", MODE = "STRATIFIED", DATA = tmpHK)



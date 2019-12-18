source("../../C2VsimCG/Rwrkspc/c2vsim_io.R")
source("npsat_functions.R")

# Read the head values
HeadAll <- c2vsim.readHeadALL("../c2vsimfg_beta2_publicrelease/C2VSimFG_BETA2_PublicRelease/Results/C2VSimFG_GW_HeadAll.out")

# Last year average
View(HeadAll[[2]])
ids <- 494:505
Hav2015 <- c2vsim.avHead(HeadAll[[3]], ids)

HK <- npsat.Input.readScattered("inputfiles/HK.npsat")
VK <- npsat.Input.readScattered("inputfiles/VK.npsat")
POR <- npsat.Input.readScattered("inputfiles/POR.npsat")

# Read Node coordinates and Mesh
XY <- c2vsim.readNodes("../c2vsimfg_beta2_publicrelease/C2VSimFG_BETA2_PublicRelease/Preprocessor/C2VSimFG_Nodes.dat", ND = 30179, Nskip = 90)
MSH <- c2vsim.readMesh("../c2vsimfg_beta2_publicrelease/C2VSimFG_BETA2_PublicRelease/Preprocessor/C2VSimFG_Elements.dat", NE = 32537, Nskip = 142, Ncols = 6)



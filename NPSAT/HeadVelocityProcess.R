source("../../C2VsimCG/Rwrkspc/c2vsim_io.R")

# Read the head values
HeadAll <- c2vsim.readHeadALL("../c2vsimfg_beta2_publicrelease/C2VSimFG_BETA2_PublicRelease/Results/C2VSimFG_GW_HeadAll.out")

# Last year average
View(HeadAll[[2]])
ids <- 494:505
Hav2015 <- c2vsim.avHead(HeadAll[[3]], ids)



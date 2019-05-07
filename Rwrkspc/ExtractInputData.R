c2vsim_path <- "f:/UCDAVIS/C2VSIM_FG_OR/C2Vsim_FG_v2/C2VSimFG-BETA_PublicRelease/"

## Saturated Hydraulic properties
NlinesSkip <- 455
nNodes <- 30179*4
Ktemp <- read.table(file = paste0(c2vsim_path, "Simulation/Groundwater/C2VSimFG_Groundwater1974.dat"),
                 header = FALSE, sep = "", skip = NlinesSkip, nrows = nNodes,
                 quote = "",fill = TRUE,
                 col.names = c("ID", "PKH", "PS", "PN", "PV", "PL"))

# KXY and KZ are the matrices with the horizontal and vertical hydrulic conductivity values [Nnodes x Nlayers]
KXY <- matrix(data = NA, nrow = 30179, ncol = 4)
KZ <- matrix(data = NA, nrow = 30179, ncol = 4)
KXY[,1] <- Ktemp$PKH[seq(from = 1, to = 30179*4, by = 4)]
KZ[,1] <- Ktemp$PL[seq(from = 1, to = 30179*4, by = 4)]
for (i in 1:3) {
  KXY[,1+i] <- Ktemp$ID[seq(from = 1, to = 30179*4, by = 4)+i]
  KZ[,1+i] <- Ktemp$PV[seq(from = 1, to = 30179*4, by = 4)+i]
}


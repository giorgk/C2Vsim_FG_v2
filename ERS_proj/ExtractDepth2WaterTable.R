
c2vsim_path <-"f:/UCDAVIS/C2VSIM_FG_OR/C2Vsim_FG_v2/C2VSimFG-BETA_PublicRelease/"

# First we need the XY coordinates, the mesh elements and the Stratigraphy information

XY <- read.table(file = paste0(c2vsim_path, "Preprocessor/C2VSimFG_Nodes.dat"),
                 header = FALSE, sep = "", skip = 90, nrows = 30179,
                 quote = "",fill = TRUE,
                 col.names = c("ID", "X", "Y"))

MSH <- read.table(file = paste0(c2vsim_path, "Preprocessor/C2VSimFG_Elements.dat"),
                  header = FALSE, sep = "", skip = 142, nrows = 32537,
                  quote = "",fill = TRUE,
                  col.names = c("ID", "nd1", "nd2", "nd3", "nd4", "S"))

strat <- read.table(file = paste0(c2vsim_path, "Preprocessor/C2VSimFG_Stratigraphy.dat"),
                    header = FALSE, sep = "", skip = 105, nrows = 30179,
                    quote = "",fill = TRUE,
                    col.names = c("ID", "GSE", "A1", "L1", "A2", "L2", "A3", "L3", "A4", "L4"))

# Read the head file
HEAD <- read.table(file = paste0(c2vsim_path, "Results/C2VSimFG_GW_HeadAll.out"),
                   header = FALSE, sep = "", skip = 6, nrows = 504*4,fill = TRUE, comment.char="",
                   stringsAsFactors=FALSE)

Lay1_id <- seq(from = 1, to = 505*4, by = 4)
Lay1_HEAD <- HEAD[Lay1_id,]
LastYearHEAD <- Lay1_HEAD[494:505,]
LastYearHEAD <- LastYearHEAD[,-1]
LastYearHEADAv <- apply(LastYearHEAD,2,mean)
DGW2015 <- (strat$GSE - LastYearHEADAv)*0.3048

save(DGW2015, file = "f:/UCDAVIS/C2VSIM_FG_OR/C2Vsim_FG_v2/ERS_proj/DGW2015.Rdata")


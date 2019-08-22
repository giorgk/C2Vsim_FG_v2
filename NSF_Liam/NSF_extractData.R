library(pracma)
library(sf)

c2vsim_path <- "../C2VSimFG-BETA_PublicRelease/"

plss <- read_sf("Liam_legal_to_latlong.shp")
plss_geom <- st_geometry(plss)

NlinesSkip <- 90
nNodes <- 30179
XY <- read.table(file = paste0(c2vsim_path, "Preprocessor/C2VSimFG_Nodes.dat"),
                 header = FALSE, sep = "", skip = NlinesSkip, nrows = nNodes,
                 quote = "",fill = TRUE,
                 col.names = c("ID", "X", "Y"))

HEAD <- read.table(file = paste0(c2vsim_path, "Results/C2VSimFG_GW_HeadAll.out"),
                   header = FALSE, sep = "", skip = 6, nrows = 504*4,fill = TRUE, comment.char="",
                   stringsAsFactors=FALSE)

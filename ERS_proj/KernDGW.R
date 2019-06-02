library(pracma)
library(xlsx)
library(sf)

load(file = "f:/UCDAVIS/C2VSIM_FG_OR/C2Vsim_FG_v2/ERS_proj/DGW2015.Rdata")
script_path = "F:/UCDAVIS/C2VSIM_FG_OR/C2Vsim_FG_v2/ERS_proj/"
c2vsim_path <- "f:/UCDAVIS/C2VSIM_FG_OR/C2Vsim_FG_v2/C2VSimFG-BETA_PublicRelease/"

# Read the Kern subregions
kern <- read_sf(paste0(script_path, "gis_data/Kern_only.shp"))
kern_geom <- st_geometry(kern)

# Read the mesh elements and node coordinates
NlinesSkip <- 90
nNodes <- 30179
XY <- read.table(file = paste0(c2vsim_path, "Preprocessor/C2VSimFG_Nodes.dat"),
                 header = FALSE, sep = "", skip = NlinesSkip, nrows = nNodes,
                 quote = "",fill = TRUE,
                 col.names = c("ID", "X", "Y"))

NlinesSkip <- 142
nNodes <- 32537
MSH <- read.table(file = paste0(c2vsim_path, "Preprocessor/C2VSimFG_Elements.dat"),
                  header = FALSE, sep = "", skip = NlinesSkip, nrows = nNodes,
                  quote = "",fill = TRUE,
                  col.names = c("ID", "nd1", "nd2", "nd3", "nd4", "S"))

#calculate element barycenters
quad_el <- which(MSH$nd4 != 0, arr.ind = FALSE)
tri_el <- which(MSH$nd4 == 0, arr.ind = FALSE)
cc <- matrix(data = 0, nrow = dim(MSH)[1], ncol = 2)
cc[quad_el,1] <- (XY$X[MSH$nd1[quad_el]] + XY$X[MSH$nd2[quad_el]] + XY$X[MSH$nd3[quad_el]] + XY$X[MSH$nd4[quad_el]])/4
cc[quad_el,2] <- (XY$Y[MSH$nd1[quad_el]] + XY$Y[MSH$nd2[quad_el]] + XY$Y[MSH$nd3[quad_el]] + XY$Y[MSH$nd4[quad_el]])/4
cc[tri_el,1] <- (XY$X[MSH$nd1[tri_el]] + XY$X[MSH$nd2[tri_el]] + XY$X[MSH$nd3[tri_el]])/3
cc[tri_el,2] <- (XY$Y[MSH$nd1[tri_el]] + XY$Y[MSH$nd2[tri_el]] + XY$Y[MSH$nd3[tri_el]])/3

# for each subregion find the mesh element ids
kern_elems <- vector("list", dim(kern)[1]) 

for (i in 1:dim(kern)[1]) {
  kern_elems[[i]] <- which(inpolygon(cc[,1], cc[,2], kern_geom[[i]][[1]][[1]][,1], kern_geom[[i]][[1]][[1]][,2]) == TRUE, arr.ind = FALSE)
}


AvDGW <- vector(mode = "numeric",length = dim(kern)[1])

for (i in 1:dim(kern)[1]) {
  el_ids <- kern_elems[[i]]
  subRegNodes <- MSH[el_ids,2:5]
  tmpnds <- as.vector(t(subRegNodes))
  subRegNodes <- unique(tmpnds)
  i_zero <- which(subRegNodes == 0)
  if (length(i_zero) > 0)
    subRegNodes <- subRegNodes[-i_zero]
  AvDGW[i] <- mean(DGW2015[subRegNodes])
}


dgw_df <- data.frame(AvDGW, row.names = kern$KernDistri)
names(dgw_df)[names(dgw_df) == "AvDGW"] <- "DGW2015"
write.xlsx(x = dgw_df, file = paste0(script_path,"Depth2GW_2015.xlsx"), sheetName = "Kern", append = TRUE)

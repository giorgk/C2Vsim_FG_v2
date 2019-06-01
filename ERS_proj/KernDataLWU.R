library(h5)
library(sf)
library(pracma)
library(xlsx)

# Paths
script_path = "F:/UCDAVIS/C2VSIM_FG_OR/C2Vsim_FG_v2/ERS_proj/"
c2vsim_path <- "f:/UCDAVIS/C2VSIM_FG_OR/C2Vsim_FG_v2/C2VSimFG-BETA_PublicRelease/"

results_path <- paste0(c2vsim_path, "Results/")
preproc_path <- paste0(c2vsim_path, "Preprocessor/")
setwd(script_path)
NtimeSteps <- 504
Nyears <- NtimeSteps/12

# Read the Kern subregions
kern <- read_sf(paste0(script_path, "gis_data/Kern_only.shp"))
kern_geom <- st_geometry(kern)

# Read the mesh elements and node coordinates
NlinesSkip <- 90
nNodes <- 30179
XY <- read.table(file = paste0(preproc_path, "C2VSimFG_Nodes.dat"),
                 header = FALSE, sep = "", skip = NlinesSkip, nrows = nNodes,
                 quote = "",fill = TRUE,
                 col.names = c("ID", "X", "Y"))

NlinesSkip <- 142
nNodes <- 32537
MSH <- read.table(file = paste0(preproc_path, "C2VSimFG_Elements.dat"),
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

# Simulation Time series
#simTime <- seq.Date(from = as.Date("1973/10/30"),to = as.Date("2015/09/30"),by = "month")


# Read L&WU Budget file
LWU_BDGinfo <- h5file(name =  paste0(results_path, "C2VSimFG_L&WU_ZBudget.hdf"), mode = "r")
hfileGroups <- list.groups(LWU_BDGinfo)
hfileDataSets <- list.datasets(LWU_BDGinfo)

# Read data and aggregate per Kern region
data_ids <- 15:53
#agg_data <- matrix(data = NA, nrow = dim(kern)[1], ncol = length(data_ids))
agg_vec <- array(data = NA, dim = c(Nyears, length(data_ids), dim(kern)[1] ))


for (i in 1:length(data_ids)) {
  temp <- LWU_BDGinfo[hfileDataSets[data_ids[i]]]
  if (dim(temp)[2] == 0)
    next
  
  splited_string <- strsplit(hfileDataSets[data_ids[i]], split = ' ')[[1]]
  
  
  for (j in 1:dim(kern)[1]){
    el_ids <- kern_elems[[j]] # In this script el_ids are the actual number of element ids
    #Extract the data for the subregion only [Ntimes X Nelem]
    subregData <- temp[,el_ids]
    # For each element average the values per water year
    YearlySubregData <- matrix(data = NA, nrow = Nyears, ncol = dim(subregData)[2])
    for (k in 1:Nyears) {
      ks <- 1+(k-1)*12
      ke <- k*12
      if (length(el_ids) == 1)
        YearlySubregData[k,] <- mean(subregData[ks:ke])
      else 
        YearlySubregData[k,] <- apply(subregData[ks:ke,],2,mean)
    }
    
    if (splited_string[length(splited_string)] == "Area"){
      agg_vec[,i,j] <- apply(YearlySubregData, 1, sum)
    }
    else{
      agg_vec[,i,j] <- apply(YearlySubregData, 1, mean)
    }
  }
  
  #temp_mat <- matrix(data = NA, nrow = dim(temp)[2], ncol = dim(temp)[1])
  #for (j in 1:dim(temp)[1]) {
  #  temp_mat[,j] <- temp[j,]
  #}
  
  
  #for (j in 1:dim(kern)[1]) {
  #  agg_data[j,i] <- sum(temp_mat[kern_elems[[j]],])
  #}
}

agg_data <- matrix(data = NA, nrow = Nyears, ncol = length(data_ids))
agg_df <- data.frame(agg_data, row.names = 1974:2015)
for (i in 1:length(data_ids)) {
  names(agg_df)[names(agg_df) == paste0("X", as.character(i))] <- strsplit(hfileDataSets[data_ids[i]], split = '/')[[1]][3]
}

for (i in 1:dim(kern)[1]){
    agg_df[,] <- agg_vec[,,i]
  if (i == 1)
    write.xlsx(x = agg_df, file = "KernAnnualAvLWU.xlsx", sheetName = kern$KernDistri[i], append = FALSE)
  else
    write.xlsx(x = agg_df, file = "KernAnnualAvLWU.xlsx", sheetName = kern$KernDistri[i], append = TRUE)
}

## convert data to data frame and write them to excel
#agg_df <- data.frame(agg_data, row.names = kern$KernDistri)
#for (i in 1:length(data_ids)) {
#  names(agg_df)[names(agg_df) == paste0("X", as.character(i))] <- strsplit(hfileDataSets[data_ids[i]], split = '/')[[1]][3]
#}
#
#write.csv(agg_df, file = paste0(script_path, "KernLWU.csv"))


# Paths
script_path = "F:/UCDAVIS/C2VSIM_FG_OR/C2Vsim_FG_v2/ERS_proj/"
c2vsim_path <- "f:/UCDAVIS/C2VSIM_FG_OR/C2Vsim_FG_v2/C2VSimFG-BETA_PublicRelease/"
results_path <- paste0(c2vsim_path, "Results/")

# Kings elements subregions are defined in a spreadsheet
Kings <- read.csv(file = paste0(script_path, "KingsSubregions.csv"))
subreg_list <- factor(Kings$Subregion)


# Read L&WU Budget file
LWU_BDGinfo <- h5file(name =  paste0(results_path, "C2VSimFG_L&WU_ZBudget.hdf"), mode = "r")
hfileGroups <- list.groups(LWU_BDGinfo)
hfileDataSets <- list.datasets(LWU_BDGinfo)


data_ids <- 15:53
agg_data <- matrix(data = NA, nrow = length(levels(subreg_list)), ncol = length(data_ids))
for (i in 1:length(data_ids)) {
  temp <- LWU_BDGinfo[hfileDataSets[data_ids[i]]]
  if (dim(temp)[2] == 0)
    next
  
  temp_mat <- matrix(data = NA, nrow = dim(temp)[2], ncol = dim(temp)[1])
  for (j in 1:dim(temp)[1]) {
    temp_mat[,j] <- temp[j,]
  }
  
  for (j in 1:length(levels(subreg_list))) {
    el_ids <- which(Kings$Subregion == levels(subreg_list)[j], arr.ind = FALSE)
    agg_data[j,i] <- sum(temp_mat[el_ids,])
  }
}

# convert data to data frame and write them to excel
agg_df <- data.frame(agg_data, row.names = levels(subreg_list))
for (i in 1:length(data_ids)) {
  names(agg_df)[names(agg_df) == paste0("X", as.character(i))] <- strsplit(hfileDataSets[data_ids[i]], split = '/')[[1]][3]
}

write.csv(agg_df, file = paste0(script_path, "KingsLWU.csv"))

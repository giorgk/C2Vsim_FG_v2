library(pracma)
library(xlsx)

load(file = "f:/UCDAVIS/C2VSIM_FG_OR/C2Vsim_FG_v2/ERS_proj/DGW2015.Rdata")
script_path = "F:/UCDAVIS/C2VSIM_FG_OR/C2Vsim_FG_v2/ERS_proj/"
c2vsim_path <- "f:/UCDAVIS/C2VSIM_FG_OR/C2Vsim_FG_v2/C2VSimFG-BETA_PublicRelease/"

# Kings elements subregions are defined in a spreadsheet
Kings <- read.csv(file = paste0(script_path, "KingsSubregions.csv"))
subreg_list <- factor(Kings$Subregion)

MSH <- read.table(file = paste0(c2vsim_path, "Preprocessor/C2VSimFG_Elements.dat"),
                  header = FALSE, sep = "", skip = 142, nrows = 32537,
                  quote = "",fill = TRUE,
                  col.names = c("ID", "nd1", "nd2", "nd3", "nd4", "S"))

AvDGW <- vector(mode = "numeric",length = length(levels(subreg_list)))

for (i in 1:length(levels(subreg_list))){
  el_ids <- which(Kings$Subregion == levels(subreg_list)[i], arr.ind = FALSE)
  subRegNodes <- MSH[Kings$ElementID[el_ids],2:5]
  tmpnds <- as.vector(t(subRegNodes))
  subRegNodes <- unique(tmpnds)
  i_zero <- which(subRegNodes == 0)
  if (length(i_zero) > 0)
    subRegNodes <- subRegNodes[-i_zero]
  AvDGW[i] <- mean(DGW2015[subRegNodes])
}


# Convert to dataframe and write to excel
dgw_df <- data.frame(AvDGW, row.names = levels(subreg_list))
names(dgw_df)[names(dgw_df) == "AvDGW"] <- "DGW2015"
write.xlsx(x = dgw_df, file = paste0(script_path,"Depth2GW_2015.xlsx"), sheetName = "Kings", append = FALSE)

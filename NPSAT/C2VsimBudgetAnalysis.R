library(hdf5r)
library(pracma)
library(Hmisc)

c2vsim_path <- "../../C2Vsim_FG_v2/C2VSimFG-BETA_PublicRelease/"
results_path <- paste0(c2vsim_path, "Results/")
c2vsim_tm <- seq.Date(from = as.Date("1973/10/1"),to = as.Date("2015/09/1"),by = "month")
ndays <- monthDays(c2vsim_tm)

GW_BDGinfo <- H5File$new(filename = paste0(results_path, "C2VSimFG_GW_ZBudget.hdf"),mode = "r+")
info <- GW_BDGinfo$ls(recursive = T)
Attribute_names <- GW_BDGinfo[[info$name[6]]][]
# The units in the hdf files are in ft, ft^2, ft^3

# Layer 1 --------------
# indices matrix
lay1_ids <- GW_BDGinfo[[info$name[7]]][,]

### PUMPING 
PMP <- array(data = 0, dim = c(32537, 504))

#### Element Pumping ++++++++++++
# this is the row number for the 	Pumping by Element_Outflow (-) in attribute names
elem_pmp_att_id <- 20
# The element i will get its pumping amount from the row ele_pmp_ids[i]
ele_pmp_ids <- lay1_ids[,elem_pmp_att_id]

tmp_id <- which(ele_pmp_ids != 0)
# the row number in info for Layer_1/Pumping by Element_Outflow (-)
ele_pmp_data_id <- 34
PMP[tmp_id,] <- PMP[tmp_id,] + GW_BDGinfo[[info$name[ele_pmp_data_id]]][,]

#### Well Pumping ++++++++++++
# this is the row number for the Pumping by Well_Outflow (-) in attribute names
well_pmp_att_id <- 22
# The element i will get its pumping amount from the row ele_pmp_ids[i]
well_pmp_ids <- lay1_ids[,well_pmp_att_id]

tmp_id <- which(well_pmp_ids != 0)
# the row number in info for 	Layer_1/Pumping by Well_Outflow (-)
well_pmp_data_id <- 36
PMP[tmp_id,] <- PMP[tmp_id,] + GW_BDGinfo[[info$name[well_pmp_data_id]]][,]

### STREAM INFLOWS 
STRM <- array(data = 0, dim = c(32537, 504))
#### Streams_Inflow (+)  ++++++++++++
strm_P_att_id <- 3
strm_P_ids <- lay1_ids[,strm_P_att_id]

tmp_id <- which(strm_P_ids != 0)
strm_P_ids_data_id <- 44
STRM[tmp_id,] <- STRM[tmp_id,] + GW_BDGinfo[[info$name[strm_P_ids_data_id]]][,]

#### Streams_Inflow (-)  ++++++++++++
strm_M_att_id <- 4
strm_M_ids <- lay1_ids[,strm_M_att_id]

tmp_id <- which(strm_M_ids != 0)
strm_M_ids_data_id <- 45
STRM[tmp_id,] <- STRM[tmp_id,] - GW_BDGinfo[[info$name[strm_M_ids_data_id]]][,]

### RECHARGE 
RCH <- array(data = 0, dim = c(32537, 504))
# Deep percolation
rch_att_id <- 7 # Deep Percolation_Inflow (+) 
rch_ids <- lay1_ids[,rch_att_id]

tmp_id <- which(rch_ids != 0)
rch_ids_data_id <- 26
RCH[tmp_id,] <- RCH[tmp_id,] + GW_BDGinfo[[info$name[rch_ids_data_id]]][,]

#Small Watershed Percolation_Inflow (+)
#swp_att_id <- 13 # Small Watershed Percolation_Inflow (+)
#swp_ids <- lay1_ids[,swp_att_id]

#tmp_id <- which(swp_ids != 0)
#swp_ids_data_id <- 41
#RCH[tmp_id,] <- RCH[tmp_id,] + GW_BDGinfo[[info$name[swp_ids_data_id]]][,]





# Layer 2 --------------
# indices matrix
lay2_ids <- GW_BDGinfo[[info$name[8]]][,]

#### Element Pumping ++++++++++++
# this is the row number for the 	Pumping by Element_Outflow (-) in attribute names
elem_pmp_att_id <- 20
# The element i will get its pumping amount from the row ele_pmp_ids[i]
ele_pmp_ids <- lay2_ids[,elem_pmp_att_id]

tmp_id <- which(ele_pmp_ids != 0)
# the row number in info for Layer_1/Pumping by Element_Outflow (-)
ele_pmp_data_id <- 62
PMP[tmp_id,] <- PMP[tmp_id,] + GW_BDGinfo[[info$name[ele_pmp_data_id]]][,]

#### Well Pumping ++++++++++++
# this is the row number for the Pumping by Well_Outflow (-) in attribute names
well_pmp_att_id <- 22
# The element i will get its pumping amount from the row ele_pmp_ids[i]
well_pmp_ids <- lay2_ids[,well_pmp_att_id]

tmp_id <- which(well_pmp_ids != 0)
# the row number in info for 	Layer_1/Pumping by Well_Outflow (-)
well_pmp_data_id <- 64
PMP[tmp_id,] <- PMP[tmp_id,] + GW_BDGinfo[[info$name[well_pmp_data_id]]][,]

# Layer 3 --------------
# indices matrix
lay3_ids <- GW_BDGinfo[[info$name[9]]][,]

#### Element Pumping ++++++++++++
# this is the row number for the 	Pumping by Element_Outflow (-) in attribute names
elem_pmp_att_id <- 20
# The element i will get its pumping amount from the row ele_pmp_ids[i]
ele_pmp_ids <- lay3_ids[,elem_pmp_att_id]

tmp_id <- which(ele_pmp_ids != 0)
# the row number in info for Layer_1/Pumping by Element_Outflow (-)
ele_pmp_data_id <- 90
PMP[tmp_id,] <- PMP[tmp_id,] + GW_BDGinfo[[info$name[ele_pmp_data_id]]][,]

#### Well Pumping ++++++++++++
# this is the row number for the Pumping by Well_Outflow (-) in attribute names
well_pmp_att_id <- 22
# The element i will get its pumping amount from the row ele_pmp_ids[i]
well_pmp_ids <- lay3_ids[,well_pmp_att_id]

tmp_id <- which(well_pmp_ids != 0)
# the row number in info for 	Layer_1/Pumping by Well_Outflow (-)
well_pmp_data_id <- 92
PMP[tmp_id,] <- PMP[tmp_id,] + GW_BDGinfo[[info$name[well_pmp_data_id]]][,]

# Layer 4 --------------
# indices matrix
lay4_ids <- GW_BDGinfo[[info$name[10]]][,]

#### Element Pumping ++++++++++++
# this is the row number for the 	Pumping by Element_Outflow (-) in attribute names
elem_pmp_att_id <- 20
# The element i will get its pumping amount from the row ele_pmp_ids[i]
ele_pmp_ids <- lay4_ids[,elem_pmp_att_id]

tmp_id <- which(ele_pmp_ids != 0)
# the row number in info for Layer_1/Pumping by Element_Outflow (-)
ele_pmp_data_id <- 118
PMP[tmp_id,] <- PMP[tmp_id,] + GW_BDGinfo[[info$name[ele_pmp_data_id]]][,]

#### Well Pumping ++++++++++++
# this is the row number for the Pumping by Well_Outflow (-) in attribute names
well_pmp_att_id <- 22
# The element i will get its pumping amount from the row ele_pmp_ids[i]
well_pmp_ids <- lay4_ids[,well_pmp_att_id]

tmp_id <- which(well_pmp_ids != 0)
# the row number in info for 	Layer_1/Pumping by Well_Outflow (-)
well_pmp_data_id <- 120
PMP[tmp_id,] <- PMP[tmp_id,] + GW_BDGinfo[[info$name[well_pmp_data_id]]][,]



# CALCULATE ANNUAL BUDGETS ----------

rch_monthly <- colSums(RCH)
pmp_monthly <- colSums(PMP) 
strm_monthly <- colSums(STRM)


rch_yearly <- colSums(array(data = rch_monthly, dim = c(12,42)))
pmp_yearly <- colSums(array(data = pmp_monthly, dim = c(12,42)))
strm_yearly <- colSums(array(data = strm_monthly, dim = c(12,42)))



CVBD <- c2vsim.readGWBUD(filename = paste0(results_path, "C2VSimFG_GW_Budget.bud"),
                         Nsub = 21, Nskip = c(8,rep(9,20)), NtimeSteps = 504,CG = F)

cumCVBD <- c2vsim.cumGWBUD(CVBD)

sum(colSums(array(data = cumCVBD$PMP, dim = c(12,42) ))/1000000)

CVBDCG <- c2vsim.readGWBUD(filename = "f:/UCDAVIS/C2VsimCG/c2vsim_cg_1921ic_r374_rev/C2VSim_CG_1921IC_R374_rev/Results/CVground.BUD",
                         Nsub = 21, Nskip = 8, NtimeSteps = 1056,CG = T)

cumCVBDCG <- c2vsim.cumGWBUD(CVBDCG)

plot(colSums(array(cumCVBDCG$NDP,c(12,88)))/1000000)


Time <- 1922:2015
pmp_FG_BD <- matrix(data = NA, nrow = length(Time), ncol = 1)
pmp_FG_BD[53:94] <- colSums(array(data = cumCVBD$PMP, dim = c(12,42) ))/1000000
rch_FG_BD <- matrix(data = NA, nrow = length(Time), ncol = 1)
rch_FG_BD[53:94] <- colSums(array(data = cumCVBD$DP, dim = c(12,42) ))/1000000
strm_FG_BD <- matrix(data = NA, nrow = length(Time), ncol = 1)
strm_FG_BD[53:94] <- colSums(array(data = cumCVBD$GFS, dim = c(12,42) ))/1000000


pmp_FG_my <- matrix(data = NA, nrow = length(Time), ncol = 1)
pmp_FG_my[53:94] <- pmp_yearly/43559.9/1000000
rch_FG_my <- matrix(data = NA, nrow = length(Time), ncol = 1)
rch_FG_my[53:94] <- rch_yearly/43559.9/1000000
strm_FG_my <- matrix(data = NA, nrow = length(Time), ncol = 1)
strm_FG_my[53:94] <- strm_yearly/43559.9/1000000

pmp_CG_BD <- matrix(data = NA, nrow = length(Time), ncol = 1)
pmp_CG_BD[1:88] <- colSums(array(data = cumCVBDCG$P, dim = c(12,88) ))/1000000
rch_CG_BD <- matrix(data = NA, nrow = length(Time), ncol = 1)
rch_CG_BD[1:88] <- colSums(array(data = cumCVBDCG$NDP, dim = c(12,88) ))/1000000
strm_CG_BD <- matrix(data = NA, nrow = length(Time), ncol = 1)
strm_CG_BD[1:88] <- colSums(array(data = cumCVBDCG$GFS, dim = c(12,88) ))/1000000

plot_df <- data.frame(Time,pmp_FG_BD, rch_FG_BD, strm_FG_BD,
                           pmp_FG_my, rch_FG_my, strm_FG_my,
                           pmp_CG_BD, rch_CG_BD, strm_CG_BD)

ggplot(plot_df, aes(x = Time),linetype = variable )+
  geom_step(aes(y = pmp_FG_BD, color = "Pumping (Fine)", linetype="Pumping (Fine)"), size = 0.8)+
  geom_step(aes(y = rch_FG_BD, color = "Recharge (fine)", linetype="Recharge (fine)"), size = 0.8)+
  geom_step(aes(y = strm_FG_BD, color = "Stream leackage (Fine)", linetype="Stream leackage (Fine)"), size = 0.8)+
  #geom_step(aes(y = pmp_FG_my, color = "C2VSIM Fine Pumping V1"), size = 1.2)+
  #geom_step(aes(y = rch_FG_my, color = "C2VSIM Fine Recharge v1"), size = 1.2)+
  #geom_step(aes(y = strm_FG_my, color = "C2VSIM Fine Stream leackage v1"), size = 1.2)+
  geom_step(aes(y = pmp_CG_BD, color = "Pumping (Coarse)", linetype="Pumping (Coarse)"), size = 1.2) + #, linetype="dotted"
  geom_step(aes(y = rch_CG_BD, color = "Recharge (Coarse)", linetype="Recharge (Coarse)"), size = 1.2)+
  geom_step(aes(y = strm_CG_BD, color = "Stream leackage (Coarse)", linetype="Stream leackage (Coarse)"), size = 1.2)+
  scale_color_manual("", breaks = c("Pumping (Coarse)", "Pumping (Fine)", 
                                    "Recharge (Coarse)", "Recharge (fine)",
                                    "Stream leackage (Coarse)", "Stream leackage (Fine)"),  
                     values = c('#619cff','#436db2','#f8766d','#bd5a53','#00bfc4','#00999d'))+
  scale_linetype_manual("",breaks = c("Pumping (Coarse)", "Pumping (Fine)", 
                                   "Recharge (Coarse)", "Recharge (fine)",
                                   "Stream leackage (Coarse)", "Stream leackage (Fine)"),
                        values=rep(c("dotted", "solid"),3))+
  #theme_minimal()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(colour = "#E3E3E3"),
        panel.grid.minor = element_line(colour = "grey", linetype = 3),
        plot.title = element_text(size = 24, hjust = 0.5, vjust=2.12),
        axis.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, hjust = 0.5, vjust = 1.12),
        axis.title.y = element_text(size = 16, hjust = 0.5, vjust = 2.12 , color = "black"),
        legend.position = c(0.20, 0.77),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1.5,"line"),
        legend.key=element_blank(),
        legend.justification = "center")+
  labs(x="Year", y="MAF")

ggsave(filename = "C2VSIM_Budget_annual.png", plot = g, device="png")  

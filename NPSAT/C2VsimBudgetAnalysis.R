library(hdf5r)
library(pracma)
library(Hmisc)
library(rgdal)
library(plotly)

c2vsim_path <- "../c2vsimfg_beta2_publicrelease/C2VSimFG_BETA2_PublicRelease/"
results_path <- paste0(c2vsim_path, "Results/")
c2vsim_tm <- seq.Date(from = as.Date("1973/10/1"),to = as.Date("2015/09/1"),by = "month")
ndays <- monthDays(c2vsim_tm)

GW_BUD <- H5File$new(filename = paste0(results_path, "C2VSimFG_GW_ZBudget.hdf"),mode = "r+")
info <- GW_BUD$ls(recursive = T)
Attribute_names <- GW_BUD[[info$name[6]]][]
Attribute_names <- trimws(Attribute_names,which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
info_names <- info$name[]
# The units in the hdf files are in ft, ft^2, ft^3

getDatafromBudget <- function(lay, name){
  out = vector(mode = "list", length = 2)
  lay_ids <- GW_BUD[[paste0("Attributes/Layer",lay,"_ElemDataColumns")]][,]
  elem_att_id <- which(Attribute_names == name)
  ids <- lay_ids[,elem_att_id]
  tmp_id <- which(ids != 0)
  out[[1]] <- tmp_id
  data_id <- which(info_names == paste0("Layer_", lay, "/", name))
  out[[2]] <- GW_BUD[[ paste0("Layer_", lay, "/", name) ]][,]
  return(out)
}

### RECHARGE 
RCH <- array(data = 0, dim = c(32537, 504))
tmp <- getDatafromBudget(1, "Deep Percolation_Inflow (+)")
RCH[tmp[[1]],] <- RCH[tmp[[1]],] + tmp[[2]] 
C2VSIM_RCH <- RCH
save(list = c("C2VSIM_RCH"), file = "C2VSIM_RCH.RData")

### STREAMS 
STRM <- array(data = 0, dim = c(32537, 504))
tmp <- getDatafromBudget(1, "Streams_Inflow (+)")
STRM[tmp[[1]],] <- STRM[tmp[[1]],] + tmp[[2]] 
tmp <- getDatafromBudget(1, "Streams_Outflow (-)")
STRM[tmp[[1]],] <- STRM[tmp[[1]],] - tmp[[2]] 


### PUMPING 
PMP <- array(data = 0, dim = c(32537, 504))
for (i in 1:4) {
  tmp <- getDatafromBudget(i, "Pumping by Element_Outflow (-)")
  PMP[tmp[[1]],] <- PMP[tmp[[1]],] + tmp[[2]] 
  tmp <- getDatafromBudget(i, "Pumping by Well_Outflow (-)")
  PMP[tmp[[1]],] <- PMP[tmp[[1]],] + tmp[[2]] 
}



# CALCULATE ANNUAL BUDGETS ----------

rch_monthly <- colSums(RCH)
pmp_monthly <- colSums(PMP) 
strm_monthly <- colSums(STRM)


rch_yearly <- colSums(array(data = rch_monthly, dim = c(12,42)))
pmp_yearly <- colSums(array(data = pmp_monthly, dim = c(12,42)))
strm_yearly <- colSums(array(data = strm_monthly, dim = c(12,42)))


source("../../../C2VsimCG/Rwrkspc/c2vsim_io.R")
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
  #geom_step(aes(y = strm_FG_my, color = "C2VSIM Fine Stream leackage v1"), size = 1.2)
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


ggplot(plot_df, aes(x = Time))+
  geom_step(aes(y = strm_FG_my), size = 1.2,color='blue')+
  geom_step(aes(y = strm_FG_BD), size = 1, color="red")

C2VSIM_CG_BUD <- data.frame(1922:2009,pmp_CG_BD[1:88],rch_CG_BD[1:88],strm_CG_BD[1:88])
colnames(C2VSIM_CG_BUD) <- c("Time", "Pmp", "Rch", "Strm")
C2VSIM_FG_BUD <- data.frame(1974:2015, pmp_FG_BD[53:94], rch_FG_BD[53:94], strm_FG_BD[53:94])
colnames(C2VSIM_FG_BUD) <- c("Time", "Pmp", "Rch", "Strm")
save(list = c("C2VSIM_CG_BUD","C2VSIM_FG_BUD"), file = "C2VSIMBUD_MAF.RData")

# ------- Caclulate recharge rate per year  -------
# load the mesh file with the area
c2vsim_mesh <- readOGR(dsn = "../gis_data/C2Vsim_Elements.shp")
c2vsim_mesh_polys <- c2vsim_mesh@polygons
# Calculate the area of Central Valley using the shapefile to make sure that this is in m
elem_area <- vector(mode = "numeric", length = length(c2vsim_mesh_polys))
for (i in 1:length(c2vsim_mesh_polys)) {
  elem_area[i] <- c2vsim_mesh_polys[[i]]@area
}

RCH <- array(data = 0, dim = c(32537, 504))
tmp <- getDatafromBudget(1, "Deep Percolation_Inflow (+)")
RCH[tmp[[1]],] <- RCH[tmp[[1]],] + tmp[[2]]

# Convert Cubic feet to cubic meter
RCH <- RCH/35.3147
RCH_YR <- array(data = 0, dim = c(dim(RCH)[1], length(1974:2015)))
sind <- 1
for (i in 1:dim(RCH_YR)[2]) {
  eind <- sind + 11
  RCH_YR[,i] <- rowSums(RCH[,sind:eind])
  sind <- sind + 12
}



RCH_YR <- 1000*(RCH_YR/elem_area)

non_zero_rch <- which(RCH_YR[,1] != 0)
tmp <- hist(RCH_YR[non_zero_rch,1],nclass = 200)
df <- data.frame(tmp$mids,tmp$counts)
colnames(df) <- c("mids","counts")

p <- plot_ly(data = df)
for (i in 1:dim(RCH_YR)[2]) {
  non_zero_rch <- which(RCH_YR[,i] > 3 & RCH_YR[,i] < 1000)
  tmp <- hist(RCH_YR[non_zero_rch,i],nclass = 75)
  df <- data.frame(tmp$mids,tmp$counts)
  colnames(df) <- c("mids","counts")
  df$counts <- 100*(tmp$counts /sum(tmp$counts))
  p <- add_trace(p, data = df, x=~mids,y=~counts,
                 type = 'scatter', mode = 'none', fill = 'tonexty')
  #p <- add_trace(p, data = df, x=rep(1961+i,length(tmp$mids)),y=~mids,z=~counts,
  #               type = 'scatter3d', mode = 'lines', fill = 'tonexty')
}

p1 <- p %>%
  layout(title = 'C2VSIM FG recharge distributions',
         xaxis = list(title = "Recharge [mm/year]", titlefont = list(size = 18), zeroline = F),
         yaxis = list(title = "Percentage of grid cells [%]", titlefont = list(size = 18)), 
         showlegend = FALSE)
p1

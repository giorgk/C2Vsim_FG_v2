library(akima)
library(plotly)
library(fitdistrplus)

source("../../../C2VsimCG/Rwrkspc/c2vsim_io.R")
source("npsat_functions.R")

# Read the head values
HeadAll <- c2vsim.readHeadALL("../c2vsimfg_beta2_publicrelease/C2VSimFG_BETA2_PublicRelease/Results/C2VSimFG_GW_HeadAll.out")

# Last year average
View(HeadAll[[2]])
ids <- 494:505
Hav2015 <- c2vsim.avHead(HeadAll[[3]], ids)

HK <- npsat.Input.readScattered("inputfiles/HK.npsat")
VK <- npsat.Input.readScattered("inputfiles/VK.npsat")
POR <- npsat.Input.readScattered("inputfiles/POR.npsat")

# Read Node coordinates and Mesh
XY <- c2vsim.readNodes("../c2vsimfg_beta2_publicrelease/C2VSimFG_BETA2_PublicRelease/Preprocessor/C2VSimFG_Nodes.dat", ND = 30179, Nskip = 90)
MSH <- c2vsim.readMesh("../c2vsimfg_beta2_publicrelease/C2VSimFG_BETA2_PublicRelease/Preprocessor/C2VSimFG_Elements.dat", NE = 32537, Nskip = 142, Ncols = 6)
STRAT <- c2vsim.readStrat("../c2vsimfg_beta2_publicrelease/C2VSimFG_BETA2_PublicRelease/Preprocessor/C2VSimFG_Stratigraphy.dat", Nnodes = 30179,nSkip = 105,Nlay = 4 )
# convert Stratigraphy variable to absolute elevations in m
ft2m <- 0.3048
Z <- matrix(data = NA, nrow = dim(STRAT)[1], ncol = 5)
Z[,1] <- STRAT$ELV*ft2m
Z[,2] <- Z[,1] - (STRAT$L11 + STRAT$L12)*ft2m
Z[,3] <- Z[,2] - (STRAT$L21 + STRAT$L22)*ft2m
Z[,4] <- Z[,3] - (STRAT$L31 + STRAT$L32)*ft2m
Z[,5] <- Z[,4] - (STRAT$L41 + STRAT$L42)*ft2m

# In the GeomPropData.R script run the section that creates the PXY, PKH and PN or PS
# PN is the Specific yield and PS is the specific storage. We use Specific yield as proxy for porosity



# Shape functions for hexahedral elements
# N1=0.125*(1-ksi)*(1-eta)*(1-zta); N2=0.125*(1+ksi)*(1-eta)*(1-zta);
# N3=0.125*(1+ksi)*(1+eta)*(1-zta); N4=0.125*(1-ksi)*(1+eta)*(1-zta);
# N5=0.125*(1-ksi)*(1-eta)*(1+zta); N6=0.125*(1+ksi)*(1-eta)*(1+zta);
# N7=0.125*(1+ksi)*(1+eta)*(1+zta); N8=0.125*(1-ksi)*(1+eta)*(1+zta);

# Shape functions for prism elements
# N1=ksi*(1-zta)/2; N2=eta*(1-zta)/2; N3=(1-ksi-eta)*(1-zta)/2;
# N4=ksi*(1+zta)/2; N5=eta*(1+zta)/2; N6=(1-ksi-eta)*(1+zta)/2;

# Shape function derivatives for prism


prism_quad_pnts <- matrix( data =
  c(2/3, 1/6, 1/sqrt(3), 1/6,
    1/6, 2/3, 1/sqrt(3), 1/6,
    1/6, 1/6, 1/sqrt(3), 1/6,
    2/3, 1/6, -1/sqrt(3), 1/6,
    1/6, 2/3, -1/sqrt(3), 1/6,
    1/6, 1/6, -1/sqrt(3), 1/6),
  nrow = 6, ncol = 4,byrow = T)
hex_quad_pnts <- matrix(data=
  c(-1/sqrt(3), -1/sqrt(3), -1/sqrt(3), 1/8,
    -1/sqrt(3),  1/sqrt(3), -1/sqrt(3), 1/8,
     1/sqrt(3),  1/sqrt(3), -1/sqrt(3), 1/8,
     1/sqrt(3), -1/sqrt(3), -1/sqrt(3), 1/8,
    -1/sqrt(3), -1/sqrt(3),  1/sqrt(3), 1/8,
    -1/sqrt(3),  1/sqrt(3),  1/sqrt(3), 1/8,
     1/sqrt(3),  1/sqrt(3),  1/sqrt(3), 1/8,
     1/sqrt(3), -1/sqrt(3),  1/sqrt(3), 1/8),
  nrow = 8, ncol = 4, byrow = T)



velPnts <- vector(mode = "list", length = dim(MSH)[1]*4*8)
cnt_pnt <- 1
for (i in 1:dim(MSH)[1]) { # Loop through the Mesh elements
  print(paste("Element: ",i))
  flush.console()
  if (MSH$ND4[i] == 0){
    ii <- 2:4
    quad_pnts <- prism_quad_pnts
  }
  else{
    ii <- 2:5
    quad_pnts <- hex_quad_pnts
  }
  nodeid <- as.numeric(MSH[i,ii])
  
  xp <- XY[nodeid, 2]
  yp <- XY[nodeid, 3]
  
  
  
  for (k in 1:4) { # Loop through the layers
    zt <- Z[nodeid, k]
    zb <- Z[nodeid, k+1]
    xyz <- matrix(c(xp, xp, yp, yp, zt, zb), nrow = 2*length(xp), ncol=3, byrow = F)
    
    
    
    for (j in 1:dim(quad_pnts)[1]) { # Loop through the quadrature points
      if(MSH$ND4[i] == 0){
        gn <- shapeFunDeriv_Prism(quad_pnts[j,1:3])
        N <- shapeFun_prism(quad_pnts[j,1:3])
      }
      else{
        gn <- shapeFunDeriv_Hex(quad_pnts[j,1:3])
        N <- shapeFun_hex(quad_pnts[j,1:3])
      }
      
      p_q <- t(xyz)%*%N
      kx <- interp(x = PXY[,1], y = PXY[,2], z = PKH[,k], xo = p_q[1], yo = p_q[2])
      kz <- interp(x = PXY[,1], y = PXY[,2], z = PV[,k], xo = p_q[1], yo = p_q[2])
      por <- interp(x = PXY[,1], y = PXY[,2], z = PN[,k], xo = p_q[1], yo = p_q[2])
      
      
      
      J <- gn %*% xyz
      invJ <- inv(J)
      
      vall <- matrix(data = NA, nrow = length(HeadAll[[3]]), ncol = 3)
      for (ih in 1:length(HeadAll[[3]])) { # Loop through the head time steps
        if(k == 4){
          # THis is temporary until I find how to compute the head as the last layer
          ht <- HeadAll[[3]][[ih]][nodeid,3]*ft2m
          hb <- HeadAll[[3]][[ih]][nodeid,4]*ft2m
          hb <- hb - (ht - hb)
          ht <- HeadAll[[3]][[ih]][nodeid,4]*ft2m
        }
        else{
          ht <- HeadAll[[3]][[ih]][nodeid,k]*ft2m
          hb <- HeadAll[[3]][[ih]][nodeid,k+1]*ft2m
        }
        
        
        dH <- invJ %*% (gn%*%c(ht,hb))
        vall[ih,] <- -c(por$z)*c(c(kx$z), c(kx$z), c(kz$z))*dH
      }
      velPnts[[cnt_pnt]] <- list(p = p_q, v=vall)
      cnt_pnt <- cnt_pnt + 1
    }
  } 
}
velPntsLast10yr <- vector(mode = "list", length = cnt_pnt-1)
for (i in 1:length(velPntsLast10yr)) {
  velPntsLast10yr[[i]] <- list(p=velPnts[[i]]$p, v=velPnts[[i]]$v[386:505,])
}
save(velPntsLast10yr,file = "velPntsLast10yr.RData")



iii <- round(993392*rand(1))
ix <- velPnts[[iii]]$v[266:505,1]
iy <- velPnts[[iii]]$v[266:505,2]
iz <- velPnts[[iii]]$v[266:505,3]

ix01 <- fit01(ix)
iy01 <- fit01(iy)
iz01 <- fit01(iz)

p <- plot_ly()
p <- add_trace(p, x = ix, y = iy, z=iz, type='scatter3d',mode = 'markers' )
p

gg <- linspace(0,1,6)
print(gg)
count <- array(data = NA, dim = c(length(gg)^3,4))
ii <- 1
for (i in gg) {
  for (j in gg) {
    for (k in gg) {
      dst <- sqrt((i-ix01)^2 + (j-iy01)^2 + (k-iz01)^2)
      count[ii,] <- c(i,j,k,length(which(dst < 0.1)))
      ii <- ii+1
    }
  }
}

ww <- which(count[,4]>1)
print(length(ww))

p <- plot_ly()
p <- add_trace(p, x = count[ww,1], y = count[ww,2], z=count[ww,3], color=count[ww,4]/max(count), type='scatter3d',mode = 'markers' )
p


ixx <- apply(matrix(data = ix, nrow = 12, ncol = 20), 2, mean)
iyy <- apply(matrix(data = iy, nrow = 12, ncol = 20), 2, mean)
izz <- apply(matrix(data = iz, nrow = 12, ncol = 20), 2, mean)

p <- plot_ly()
p <- add_trace(p, x = ixx, y = iyy, z=izz, type='scatter3d',mode = 'markers' )
p

p <- plot_ly()
p <- add_trace(p, x=1:length(ix), y = ix,type = 'scatter', mode = 'lines')
p <- add_trace(p, x=seq(1,length(ix),12), y = ixx, type = 'scatter', mode = 'lines')
p


p <- add_trace(p, x=1:length(ix), y = fit01(iy),type = 'scatter', mode = 'lines')
p <- add_trace(p, x=1:length(ix), y = fit01(iz),type = 'scatter', mode = 'lines')
p

plot(iz)



p <- plot_ly()
p <- add_trace(p, x = ixx, y = iyy, type = 'scatter', mode = 'markers' )
p
p <- plot_ly()
p <- add_trace(p, x = ix, y = iz, type = 'scatter', mode = 'markers' )
p

p <- plot_ly()
p <- add_trace(p, x = iy, y = iz, type = 'scatter', mode = 'markers' )
p


p <- plot_ly()
p <- add_trace(p, x = ix, y = iy, z=iz, type='scatter3d',mode = 'markers' )

# http://www.learnbymarketing.com/tutorials/linear-regression-in-r/
fit_lm <- lm(izz ~ ix + iy)
summary(fit_lm)
cc <- fit_lm$coef

bb <-  data.frame(x = c(min(ix), min(ix), max(ix), max(ix), min(ix)),
              y = c(min(iy), max(iy), max(iy), min(iy), min(iy)))
bb$z <- cc[1] + bb$x*cc[2] + bb$y*cc[3]
p <- add_trace(p, x = bb$x, y = bb$y, z=bb$z, type='scatter3d', mode = 'line' )
p



g <- vector2Angle(velPnts[[iii]]$v[266:505,])
ix <- fit01(g[,1])
iy <- fit01(g[,2])
iz <- fit01(g[,3])
#p <- add_trace(p, x = log10(g[,1]), y = log10(g[,2]), color=log10(g[,3]), type='scatter', mode = 'markers')
#p <- add_trace(p, x = log10(g[,1]), y = log10(g[,2]), z=log10(g[,3]))

ix <- fit01(velPnts[[iii]]$v[266:505,1])
iy <- fit01(velPnts[[iii]]$v[266:505,2])
iz <- fit01(velPnts[[iii]]$v[266:505,3])
# 3d data
# 2D data
p <- plot_ly()
p <- add_trace(p, x = ix, y = iy, color = iz, type='scatter', mode = 'markers' )
p

p <- plot_ly()
p <- add_trace(p, x = g[,1], y = g[,2], color = g[,3], type='scatter', mode = 'markers' )
p



linspace(min(ix), max(ix),20)


hist(fit_lm$resid)
qqnorm(fit_lm$resid)
qqline(fit_lm$resid)


# See the link for the analysis
# http://www.di.fc.ul.pt/~jpn/r/distributions/fitting.html#choosing-which-distribution-to-fit
my_data <- g[386:505,2] # last 20 years -> 266:505, last 10 years -> 386:505
fit <- fitdistr(my_data, densfun="normal")
hist(my_data, pch=20, breaks=25, prob=TRUE, main="")
curve(dnorm(x, fit$estimate[1], fit$estimate[2]), col="red", lwd=2, add=T)

plotdist(my_data, histo = TRUE, demp = TRUE)
descdist(my_data, discrete=FALSE, boot=500)

fit_n <- fitdist(my_data, "norm")
fit_w  <- fitdist(my_data, "Weibull")
fit_g  <- fitdist(my_data, "gamma")
fit_ln <- fitdist(my_data, "lnorm")
summary(fit_ln)

par(mfrow=c(2,2))
plot.legend <- c("normal", "Weibull", "lognormal", "gamma")
denscomp(list(fit_n, fit_w, fit_g, fit_ln), legendtext = plot.legend)
cdfcomp (list(fit_n, fit_w, fit_g, fit_ln), legendtext = plot.legend)
qqcomp  (list(fit_n, fit_w, fit_g, fit_ln), legendtext = plot.legend)
ppcomp  (list(fit_n, fit_w, fit_g, fit_ln), legendtext = plot.legend)


shapeFunDeriv_Prism <- function(n){
  xi <- n[1]
  eta <- n[2]
  zta <- n[3]
  gn <- matrix(data = NA, nrow = 3, ncol = 6) 
  gn[1,1] <- 1/2 - zta/2;
  gn[1,2] <- 0;
  gn[1,3] <- zta/2 - 1/2;
  gn[1,4] <- zta/2 + 1/2;
  gn[1,5] <- 0;
  gn[1,6] <- zta/2 - 1/2;
  
  gn[2,1] <- 0;
  gn[2,2] <- 1/2 - zta/2;
  gn[2,3] <- zta/2 - 1/2;
  gn[2,4] <- 0;
  gn[2,5] <- zta/2 + 1/2;
  gn[2,6] <- -zta/2 - 1/2;
  
  gn[3,1] <- -xi/2;
  gn[3,2] <- -eta/2;
  gn[3,3] <- eta/2 + xi/2 - 1/2;
  gn[3,4] <- xi/2;
  gn[3,5] <- eta/2;
  gn[3,6] <- 1/2 - xi/2 - eta/2;
  
  return(gn)
}

shapeFunDeriv_Hex <- function(n){
  ksi <- n[1]
  eta <- n[2]
  zta <- n[3]
  gn <- matrix(data = NA, nrow = 3, ncol = 8)
  gn[1,1] <- -((eta - 1)*(zta - 1))/8  
  gn[1,2] <-  ((eta - 1)*(zta - 1))/8
  gn[1,3] <- -((eta + 1)*(zta - 1))/8
  gn[1,4] <-  ((eta + 1)*(zta - 1))/8 
  gn[1,5] <-  ((eta - 1)*(zta + 1))/8
  gn[1,6] <- -((eta - 1)*(zta + 1))/8
  gn[1,7] <-  ((eta + 1)*(zta + 1))/8
  gn[1,8] <- -((eta + 1)*(zta + 1))/8
  
  gn[2,1] <- -((ksi - 1)*(zta - 1))/8
  gn[2,2] <-  ((ksi + 1)*(zta - 1))/8
  gn[2,3] <- -((ksi + 1)*(zta - 1))/8
  gn[2,4] <-  ((ksi - 1)*(zta - 1))/8
  gn[2,5] <-  ((ksi - 1)*(zta + 1))/8
  gn[2,6] <- -((ksi + 1)*(zta + 1))/8
  gn[2,7] <-  ((ksi + 1)*(zta + 1))/8
  gn[2,8] <- -((ksi - 1)*(zta + 1))/8
  
  gn[3,1] <- -((eta - 1)*(ksi - 1))/8
  gn[3,2] <-  ((eta - 1)*(ksi + 1))/8
  gn[3,3] <- -((eta + 1)*(ksi + 1))/8
  gn[3,4] <-  ((eta + 1)*(ksi - 1))/8
  gn[3,5] <-  ((eta - 1)*(ksi - 1))/8
  gn[3,6] <- -((eta - 1)*(ksi + 1))/8
  gn[3,7] <-  ((eta + 1)*(ksi + 1))/8
  gn[3,8] <- -((eta + 1)*(ksi - 1))/8
  
  return(gn)
}

shapeFun_prism <- function(n){
  # N1=ksi*(1-zta)/2; N2=eta*(1-zta)/2; N3=(1-ksi-eta)*(1-zta)/2;
  # N4=ksi*(1+zta)/2; N5=eta*(1+zta)/2; N6=(1-ksi-eta)*(1+zta)/2;
  N <- vector(mode = "numeric", length = 6)
  xi <- n[1]
  eta <- n[2]
  zta <- n[3]
  
  N[1] <- xi*(1-zta)/2
  N[2] <- eta*(1-zta)/2
  N[3] <- (1-xi-eta)*(1-zta)/2
  N[4] <- xi*(1+zta)/2
  N[5] <- eta*(1+zta)/2
  N[6] <- (1-xi-eta)*(1+zta)/2
  return(N)
}

shapeFun_hex <- function(n){
  # N1=0.125*(1-ksi)*(1-eta)*(1-zta); N2=0.125*(1+ksi)*(1-eta)*(1-zta);
  # N3=0.125*(1+ksi)*(1+eta)*(1-zta); N4=0.125*(1-ksi)*(1+eta)*(1-zta);
  # N5=0.125*(1-ksi)*(1-eta)*(1+zta); N6=0.125*(1+ksi)*(1-eta)*(1+zta);
  # N7=0.125*(1+ksi)*(1+eta)*(1+zta); N8=0.125*(1-ksi)*(1+eta)*(1+zta);
  N <- vector(mode = "numeric", length = 8)
  ksi <- n[1]
  eta <- n[2]
  zta <- n[3]
  N[1] <- 0.125*(1-ksi)*(1-eta)*(1-zta)
  N[2] <- 0.125*(1+ksi)*(1-eta)*(1-zta)
  N[3] <- 0.125*(1+ksi)*(1+eta)*(1-zta)
  N[4] <- 0.125*(1-ksi)*(1+eta)*(1-zta)
  N[5] <- 0.125*(1-ksi)*(1-eta)*(1+zta)
  N[6] <- 0.125*(1+ksi)*(1-eta)*(1+zta)
  N[7] <- 0.125*(1+ksi)*(1+eta)*(1+zta)
  N[8] <- 0.125*(1-ksi)*(1+eta)*(1+zta)
  return(N)
}

vector2Angle <- function(V){
  ph <- rad2deg(atan2(V[,1], V[,2]))
  neg <- which(ph<0)
  ph[neg] <- 360 + ph[neg]
  th <- rad2deg(atan2(sqrt(V[,1]^2 + V[,2]^2), V[,3]))
  m <- sqrt(V[,1]^2 + V[,2]^2 + V[,2]^2)
  return(cbind(ph, th, m))
}

fit01 <- function(x){(x-min(x))/(max(x)-min(x))}


p <- plot_ly()
p <- add_trace(p,x = xp, y = yp, z = zt)
p <- add_trace(p,x = xp, y = yp, z = zb)
p <- add_trace(p,x = p_q[1], y = p_q[2], z = p_q[3])
p


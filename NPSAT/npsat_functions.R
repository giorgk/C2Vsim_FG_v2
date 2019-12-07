npsat.writeMesh2obj <- function(filename, XY, MSH){
  con <- file(filename, open = "w")
  if(dim(XY)[2] == 2)
    XY <- cbind(XY, 0)
  write.table(cbind("v", XY), file = filename, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  id_qd = which(MSH[,4] != 0 )
  id_tr = which(MSH[,4] == 0 )
  write.table(cbind("f", MSH[id_qd,]), file = filename, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
  write.table(cbind("f", MSH[id_tr,1:3]), file = filename, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
  close(con)
}

npsat.readOBJmesh <- function(filename, maxNgon = 4){
  allLines <- readLines(filename)
  XYZ = matrix(data = NA, nrow = length(allLines), ncol = 3)
  MSH = matrix(data = NA, nrow = length(allLines), ncol = maxNgon)
  ivrt <- 1
  iel <- 1
  
  for (i in 1:length(allLines)) {
    if (strcmp(substr(allLines[i],1,1),"v")){
      # obj files write 3d coordinates
      XYZ[ivrt,] <- scan(text = substr(allLines[i],2,100), n = 3, quiet = TRUE) 
      ivrt <- ivrt + 1
    }
    if (strcmp(substr(allLines[i],1,1),"f")){
      # we dont know 
      tmp <- scan(text = substr(allLines[i],2,100), quiet = TRUE)
      MSH[iel,1:length(tmp)] <- tmp 
      iel <- iel + 1
    }
  }
  out <- vector(mode = "list", length = 2)
  ivs <- -(ivrt:length(allLines))
  ils <- -(iel:length(allLines))
  out[[1]] <- XYZ[ivs,]
  out[[2]] <- MSH[ils,]
  return(out)
}

npsat.Input.WriteMesh <- function(filename, XY, MSH){
  con <- file(filename, open = "w")
  write(paste(dim(XY)[1], dim(MSH)[1]), file = con)
  if (dim(XY)[2] == 2){
    XY <- cbind(XY, 0)
  }
  write.table(XY, file = filename, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
  write.table(MSH-1, file = filename, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
  close(con)
}

#' npsat.input.WriteScattered
#'
#' @param filename The name of the file
#' @param PDIM is the number of columns that correspond to coordinates.  This is 1 for 1D points of 2 for 2D points.
#' @param TYPE Valid options for type are FULL, HOR or VERT
#' @param MODE Valid options for mode are SIMPLE or STRATIFIED
#' @param DATA the data to be printed. Data should have as many columns as needed.
#' For example it can be :
#' [x v]
#' [x v1 z1 v2 z2 ... vn-1 zn-1 vn]
#' [x y v]
#' [x y v1 z1 v2 z2 ... vn-1 zn-1 vn]
#'
#' @return
#' @export
#'
#' @examples
#' For 2D interpolation such as top, bottom elevation or recharge
#' npsat.input.WriteScattered(filename, 2, "HOR", "SIMPLE", data)
npsat.input.WriteScattered <- function(filename, PDIM, TYPE, MODE, DATA){
  write("SCATTERED", file = filename, append = FALSE)
  write(TYPE, file = filename, append = TRUE)
  write(MODE, file = filename, append = TRUE)
  Ndata <- dim(DATA)[2] - PDIM
  write(paste(dim(DATA)[1], Ndata), file = filename, append = TRUE)
  write.table(DATA, file = filename, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
}
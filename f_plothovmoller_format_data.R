# Format data for image plot

f_plothovmoller_format_data <- function(L=list, xn=xname, yn=yname, zn=varname, xr=xrange, yr=yrange, zr=zrange, zlog=logtrans) {

  #print(zn)

  # Populate list and transpose z variable matrix if needed
  p <- list()
  p$x <- L[[xn]]
  p$y <- L[[yn]]
  p$z <- L[[zn]]
  p$xlim <- xr
  p$ylim <- yr
  p$zlim <- zr
  
  # Select center of pisces1d 3x3 horizontal grid if needed
  if (length(dim(p$z))==2) {
    #print("z variable is 2D")
  } else if (length(dim(p$z))==4 & dim(p$z)[1]==3 & dim(p$z)[2]==3) {
    #print("z variable is 4D with pisces1D 3x3 horizontal grid, selecting central cell")
    p$z <- p$z[2,2,,]
  } else {
    stop("Unexpected z variable dimensions")
  }
  
  # Transpose z if needed
  if (dim(p$z)[1]==length(p$x) & dim(p$z)[2]==length(p$y)) {
    #print("Dimensions match")
  } else if (dim(p$z)[2]==length(p$x) & dim(p$z)[1]==length(p$y)) {
    p$z <- t(p$z)
    #print("Dimensions transposed")
  } else {
    stop("Dimensions don't match")
  }
  
  # Crop if needed
  jx <- p$x >= xr[1] & p$x <= xr[2]
  jy <- p$y >= yr[1] & p$y <= yr[2]
  p$x <- p$x[jx]
  p$y <- p$y[jy]
  p$z <- p$z[jx,jy]
  
  # Log-transform z if needed
  if (zlog) {p$z <- log10(p$z); p$zlim <- log10(p$zlim)}
  
  # Cap z variable according to plot z ranges
  p$z[!is.na(p$z) & !is.nan(p$z) & p$z<=min(p$zlim)] <- min(p$zlim)
  p$z[!is.na(p$z) & !is.nan(p$z) & p$z>=max(p$zlim)] <- max(p$zlim)
  p$z[is.infinite(p$z)] <- NaN
  
  return(p)
}

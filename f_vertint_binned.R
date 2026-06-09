# Function to compute vertical integrals of binned data

# Input is
# 3D array with dimensions x = lon, y = lat, z = depth
#   depth bounds of the vertical integrationzmi:zMi
#   depth (vertical) grid with bin center and bounds (these expected in 2 x nlevels)
#   indices of selected x-y on the horizontal grid
# Operates only on vectors (depth dimension)
# The time dimension is currently not allowed, and managed through apply(), or similar, in script that calls this function

#  Martí Gali Tapias, November 2019

f_vertint_binned <- function(v = binnedvar, zm = zmin, zM = zmax, zc = zcenter, zb = zbounds, xy_sel = selected_xy_indices) {
  
  v <- as.vector(v)
  
  # Add check to ensure depth axis is positive?
  
  # Select x-y indices (usually, center of 3x3 horizontal grid)
  if (!is.null(xy_sel)) {
    v <- v[xy_sel$xi,xy_sel$yi,]
  }
  
  # Correct for weird intervals shifted by zcenter in PISCES deptht vector, if needed
  if (min(zb) < 0) {
    zb <- zb - min(zb, na.rm = T) 
  }
  
  # Indices of first and last bin. WARNING: Using bin bounds!
  zmi <- max(which(zb[1,] <= zm), na.rm = F)
  zMi <- min(which(zb[2,] >= zM), na.rm = F)
  # print(c(zmi,zMi))
  
  if (is.na(zmi) | is.na(zMi)) {
    
    vertint_out <- NA
    warning("f_vertint_binned did not find indices within desired depth bounds")
    
  } else if (length(v)==0) {
    
    vertint_out <- NA
    
  } else {
    
    # Check for NA
    zmiss <- as.numeric(is.na(v))
    countzmiss <- sum(zmiss[zmi:zMi])   # count missing values in depth interval
    lenzmiss <- length(zmiss[zmi:zMi])  # length of depth interval
    fraczmiss <- countzmiss / lenzmiss  # fraction of missing values in depth interval
    fracmax <- 0.2                      # maximum fraction of missing values allowed
 
    # If fraction of missing values is larger than threshold, do not integrate. Else, disregard NA and 
    # compensate for them when computing integrals (equivalent to linearly interpolating, but more efficient)
    
    if(fraczmiss > fracmax) {
      
      vertint_out <- NA
      warning(paste0("f_vertint_binned found >",round(as.numeric(fraczmiss*100),0),"% depth bins with NA, ",countzmiss," out of ",lenzmiss,": no output produced"))
      
    } else {
      
      # Apply depth weighting according to layer thickness
      zthick <- diff(zb)
      vthick <- v*zthick
      
      # Correct for fraction of bin comprised in selected zmin-zmax interval
      zweights <- rep(x = 0, length(v))                     # initialize vector
      zweights[(zmi+1):(zMi-1)] <- 1
      zweights[zmi] <- (zb[2,zmi] - zm) / diff(zb[,zmi])    # fraction of top depth bin included in zmin-zmax  bounds
      zweights[zMi] <- (zM - zb[1,zMi]) / diff(zb[,zMi])    # fraction of bottom depth bin included in zmin-zmax  bounds
      
      # Compute scaling factor to account for bin with NA if pertinent
      zmzw <- zmiss * zweights
      sf <- sum(zweights[zmi:zMi]) / sum(zmzw[zmi:zMi], na.rm = T)
      
      # Compute integral
      vweighted <- vthick*zweights
      vertint_out <- sum(vweighted, na.rm = T)
    }
    
  }
  
  return(vertint_out)
}

# Function to compute stats for mod (xvec) vs. obs (yvec) vectors
# NOTE: ST statistic (Jollif et al., 2009 JMS) can be computed from the outputs

f_skillstats_xvec_yvec <- function(xvar, yvar) {
  
  # Remove NA
  imatch <- !is.na(xvar) & !is.na(yvar)
  n <- sum(imatch)

  if (n) {
    x <- xvar[imatch]
    y <- yvar[imatch]
    
    # Correlation coefficient
    rp <- cor(x, y, use = "pairwise", method = "pearson")
    
    # Root mean squared error (if computed in log space, gives reliability index = 10^rmse)
    rmse <- sqrt( mean( (x - y)^2, na.rm = TRUE ) )
    
    # Model/Obs quotient (related to reliability index). Reliability index is 10^(log10 space RMSE), or 10^(sqrt(mean((model/data)^2))
    ri <- mean( (x / y), na.rm = TRUE )
    
    # Median absolute percentage error or deviation (mapd)
    mapd <- 100*median( abs((x-y)/y), na.rm = T )
    
    # Normalized standard deviation
    sdstar <- sd(x, na.rm = T) / sd(y, na.rm = T)
    
    # Bias
    bias <- mean(x - y, na.rm = T)
   
    # Relative bias (NOTE: not mean relative bias!)
    rbias <- bias / mean(y , na.rm = T)
    
    # Modelling efficiency (doi:10.1016/0022-1694(70)90255-6)
    mef <- 1 - rmse / sd(y, na.rm = T)
    
    # Kling-Gupta efficiency (doi:10.1029/2011WR010962)
    muy <- mean(y, na.rm = T)
    mux <- mean(x, na.rm = T)
    kge <- 1 - sqrt( (rp - 1)^2 + (sdstar - 1)^2 + (mux/muy - 1)^2 )
    
    # ST from Jollif et al. (2009) JMS: calculate a posteriori because it depends on normalization to maximum Bias
    
    # CV around the mean for observations only
    cvy <- sd(y, na.rm = T)/muy

    return(list(n=n, bias=bias, rbias=rbias, rp=rp, sdstar=sdstar, rmse=rmse, ri=ri, mapd=mapd, mef=mef, kge=kge, mux=mux, muy=muy, cvy=cvy))

  } else {

    return(list(n=n, bias=NA, rbias=NA, rp=NA, sdstar=NA, rmse=NA, ri=NA, mapd=NA, mef=NA, kge=NA, mux=mux, muy=muy, cvy=cvy))
  } 
}

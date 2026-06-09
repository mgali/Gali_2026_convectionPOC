# Despike profile
# Based on f_despike_profile
# x is spike time series (or vertical profile)
# k is running window width: 11 in Briggs2020 (smoother results and higher bPOC!), 7 in Briggs2011 and Lacour2019

f_despike_generic <- function(x, k, st, varname) {
  
  if ( sum(is.na(x)) == length(x) | k > length(x) ) {
    
    out <- data.frame(despiked = NA, Spike = NA, SpikeFreq = NA)
    
  } else {
    
    xrm <- caTools::runmin(x, k = k, endrule = "min") # SEE NOTE
    despiked <- caTools::runmax(xrm, k = k, endrule = "keep") # needed because otherwise layers with sharp but consistent peaks would be treated as spikes
    despiked[is.na(x)] <- NA
    # ensure spike signal is positive and within lower-upper spike thresholds
    inospike <- ( x < (despiked + st$lower) ) | ( x > (despiked + st$upper) )
    inospike[is.na(inospike)] <- FALSE # important to avoid NA indices
    x[inospike] <- despiked[inospike]
    Spike <- x - despiked
    SpikeFreq <- as.numeric(as.logical(Spike)) # spike frequency
    out <- data.frame(despiked=despiked, Spike=Spike, SpikeFreq=SpikeFreq)
    
  }
  
  # return(out)
  return(
    setNames(out,
             c(paste0(varname,"_bdespiked"),c(paste0(varname,"_bspike")),c(paste0(varname,"_bSpikeFreq")))
    )
  )
  
}

# ----
# NOTE
# endrule could be median, too. Don't know how to do it with endrule = "func".
# Doing it manually:
# xrm[(length(x) - k %/% 2):length(x)] <- median(tail(x, n = k %/% 2), na.rm = T)
# However, median will be sligthly distored if one spike is present
# Using "min" is a minor issue because measurements get then bin-averaged and boxcar-smoothed

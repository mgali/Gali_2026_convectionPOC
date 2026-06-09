# Function to compute deep convection (>1000m) onset according to BGC-Argo trajectory data
# WARNING: Not adapted for Southern hemisphere
# Threshold criterion tcrit based on bio-optical variable jump. Usually works for bbp, but not robust enough
# Avoid "false positives" by applying noise threshold criterion "bocut"

f_conv1000_onset <- function(bo, bodate, bocut) {
  
  cyear <- median(year(bodate))
  startdate <- decimal_date(bodate[1]) - cyear
  enddate <- decimal_date(bodate[length(bodate)]) - cyear
  if ( startdate > (335/365) | enddate < (90/365) ) {
    date_conv <- NA # don't compute onset of deep convection with inappropriate measurement period
  } else {
    tout <- seq(bodate[1], bodate[length(bodate)], "1 day")
    cc <- approx(bodate, bo, tout, method = "l", rule = 2, ties = mean) # Interpolate to daily, no NA in output
    csmooth <- runmed(cc$y, k = 3, ) # boxcar median smoothing
    
    if ( sum(csmooth>bocut) ) {
      date_conv <- min(cc$x[csmooth>bocut] - as.duration("1 day"), na.rm=T) # additional criterion of high values persistence is not necessary with median smoothing aplied above
      if (date_conv > force_tz(as.POSIXct(paste0(cyear,"-05-01")), "GMT") & date_conv < force_tz(as.POSIXct(paste0(cyear,"-12-01")), "GMT")) {date_conv <- NA}
    } else {
      date_conv <- NA
    }
  }
  return( date_conv )
}

# Function to compute stratification onset
f_strat_onset_v2 <- function(zmix, zmixdate, zcut) {
  
  # # Uncomment for testing
  # yi <- which(mdf$year==2016)
  # zmix = mdf[yi,mldvar]
  # zmixdate = mdf$date[yi]
  # zcut = 100
  
  zmax <- max(zmix, na.rm = T)
  syear <- median(year(zmixdate))
  startdate <- decimal_date(zmixdate[1]) - syear
  enddate <- decimal_date(zmixdate[length(zmixdate)]) - syear
  if ( zmax < 40 | startdate > (45/365) | enddate < (121/365) ) {
    date_strat <- NA # don't compute onset of stratification in case of permanent shallow stratification
  } else {
    tout <- seq(zmixdate[1], zmixdate[length(zmixdate)], "1 day")
    decdate_noyear <- decimal_date(tout)-year(tout)
    zz <- approx(zmixdate, zmix, tout, method = "l", rule = 2, ties = mean) # Interpolate to daily, no NA in output
    zmix <- zz$y
    zsmooth <- runmed(zmix, k = 31) # boxcar median smoothing (tested 7, 15 and 31)
    ihalfyear <- which.min(abs(decdate_noyear - 0.5))
    zsmooth_filled <- zsmooth
    zsmooth_filled[is.na(zsmooth) & (1:(length(zsmooth)) < ihalfyear)] <- 1000 # fill with 1000 when NA occur in first half of year
    # zcut <- median(zsmooth_filled, na.rm = T) # median of regularly interpolated MLD is adaptive threshold to compute stratification (better than fixed criterion of eg 40 or 50 m depth)
    # print(zcut)
    
    date_deep <- min(tout[which.max(zsmooth[1:ihalfyear])]) # ensure date_deep occurs during first half of year (missing MLD data in winter, due to bgc-argo z range limited to 0-1000, could result in spurious annual max MLD in fall)
    date_strat <- tout[min(which(zsmooth<=zcut & tout>date_deep))] # ensure it occurs later than deepest annual mixing
  }
  return(date_strat)
}

# # DEBUGGING
# plot(zmixdate, zmix, pch=20, cex=0.5)
# lines(tout, zsmooth, col="gray")
# abline(v = date_strat)

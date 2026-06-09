# Function to estimate date of convection termination based on latest profile with stratification weaker than a given MLD criterion (MLD input: zmix)
f_conv_termin <- function(zmix, zmixdate) {
  
  # Subset months < 6
  jm <- which(month(zmixdate) < 6)
  if (length(jm)) {
    zmix <- zmix[jm]
    zmixdate <- zmixdate[jm]
    zmax <- max(zmix, na.rm = T)
    syear <- median(year(zmixdate))
    startdate <- decimal_date(zmixdate[1]) - syear
    enddate <- decimal_date(zmixdate[length(zmixdate)]) - syear
    if ( zmax < 40 | startdate > (45/365) | enddate < (121/365) ) {
      date_strat <- NA # don't compute onset of stratification in case of permanent shallow stratification
    } else {
      zna <- which(is.na(zmix))
      date_deep <- zmixdate[!is.na(zmix) & zmix==zmax]
      if (length(date_deep)>1) {date_deep <- date_deep[1]}
      if (length(zna)) {
        date_strat <- zmixdate[max(zna)] # version 2025-05-06
        # date_strat <- zmixdate[max(zna)+1] # version 2025-05-05 and before
        if (date_strat<date_deep) {date_strat <- date_deep}
      } else {
        date_strat <- NA
      }
    }
  } else {
    date_strat <- NA
  }
  return(date_strat)
}

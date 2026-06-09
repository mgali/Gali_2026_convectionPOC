# Calculate POC/bbp700 conversion factor along the vertical profile as continuous function of bbp700 and chl
# 
# Martí Gali Tapias, May 2024
# 
# Rationale: Koestner et al. 2024 FMARS (K24)
# Inputs: vectors of bbp700 and chla
# Optional feature: if units variable supplied (currently mmolC/m2), output is converted to these units

f_K24 <- function(bbpin, chlin, units) {

  # Fixed parameters for surface-layer POC from bbp700 and chl
  # According to Table 3 of Koestner, Stramski and Reynolds (2024). Improved multivariable algorithms for estimating 
  # oceanic particulate organic carbon concentration from optical backscattering and chlorophyll-a measurements
  # Front. Mar. Sci., URL PLACEHOLDER
  k1 <- 52.82
  k2 <- 0.1353
  k3 <- 0.8849
  k4 <- 0.2268
  # POWER BIAS CORRECTION
  e1 <- 1.469
  e2 <- -0.734
  pocmin <- 36.8                                                              # mgC/m3
  
  # Play with CHLA (adjust for "Roesler factor" relating fluorescence to Chl in each region or globally)
  # chlin <- chlin*2
  
  # Calculate conversion factor profile: cfvec
  chl2bbp <- chlin/bbpin
  chl2bbp[!is.na(chl2bbp) & chl2bbp<1] <- 1
  chl2bbp[!is.na(chl2bbp) & chl2bbp>2000] <- 2000
  pocstar <- k1 * (bbpin^k2) * (chl2bbp^k3) * ( chl2bbp^(k4*log10(bbpin)) )   # mgC/m3
  ipocmin <- which(pocstar < pocmin)
  pocok <- pocstar
  pocok[ipocmin] <- 10^e2 * ((pocstar[ipocmin])^e1)                           # mgC/m3
  # Optional feature: units
  # pocok <- pocstar                                                            # do not apply correction: designed for epipelagic 
  if (units=="mmolC.m-3") {pocok <- pocok/12.011}
  return(pocok)
}


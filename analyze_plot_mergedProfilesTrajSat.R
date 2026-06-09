# Merge Argo profile and trajectory data with satellite matchup data (GlobCoLour), all of which have been previously preprocessed
# Compare equivalent or related variables across datasets
# Section each year according to convection start and end dates and produce primary convperiod* files for posterior manual editing
# Save merged datasets "data*Rda"
# If enabled, export data for years with convection events, with date centered on convection start date
#
# Input arguments: float number (FWMO), region, exporteventyear
# Input files: Preprocessed files generated directly from Coriolis-downloaded (BGC)Argo dataset using mergeByFloat_preprocess_binZ_binT.R
# Outputs: plot*summary, merged prof/traj/sat data with or without regular time interval ("regularized)
# 
# Marti Gali Tapias, April 2024, marti.gali.tapias@gmail.com, mgali@icm.csic.es

# -----------------------------------------------------------------------------------------------------
# Input arguments
# Typical combinations (usually with tgrid = "5 days")
#   To make "plotsummary": regularize <- T or F, exporteventyear <- 9999
#   To make "plotintsummary": regularize <- T, exporteventyear <- 9999
#   To export data without event centering, al years: savedata <- T, regularize <- F, plots = F
#   To export convection event centered data, year-wise: savedata <- T, regularize <- T, plots = F unless they are used for checking that time shif is correct
savedata <- F            # Enable saving data
regularize <- F          # Default: FALSE
biospike <- ""           # Spike treatment in traj files. Default = "" (Tukey's criterion), Briggs-style = "b"
docorrchl <- T           # Default: TRUE (use CORRCHL); Set to FALSE to use uncorrected Chla flourescence baseline
plotsummary <- T         # Figures S3 to S9 of convection paper SM
plotintsummary <- F      # Plot vertically-integrated data (rough plots, not included in convection paper)

# Run code in testing or batch processing mode (terminal command: Rscript analyze_plot_mergedProfilesTrajSat.R arg1 arg2)
testing <- T
if (testing) {
  fwmo <- 6901527                           # 6901480, 6901486, 6901516, 6901524, 6901527, 6901472 (STG!)
  region <- "NASPG"                         # "NASPG" (all convection study floats except 6901472), "STG" (6901472); other options are "MED" and "SO"
  exporteventyear <- 9999                   # Either 9999 or the year to be exported, with time centered on convection start date (-3 and + 9 months)
} else {
  args <- commandArgs(trailingOnly=TRUE)
  fwmo <- as.numeric(args[1])
  region <- args[2]
  exporteventyear <- 9999    # Either 9999 or the year to be exported, with time centered on convection start date (-3 and + 9 months)
}
zgrid_out <- "A18"                          # previously tested values, no longer available: L31, L75
zgrid_bin <- "A18"                          # previously tested values, no longer available: L31, L75
Torig <- "1day"
pocmethod <- "G22"                          # either G22 or K24 (K24zeu not enabled: Zeu input needed). If NULL, older files (with G22) are loaded
hgrid <- "orca1"

# -----------------------------------------------------------------------------------------------------
# Arguments used for producing output at desired time intervals
tgrid <- "5 days"                           # possible values: "Torig" and periods allowed by lubridate::floor_date: "day", "5 days", "month"... (longer than month is useless!)
tgridname <- paste(unlist(strsplit(tgrid, split = " ")), collapse = "")

# -----------------------------------------------------------------------------------------------------
# Define paths and names of input variable files
spath <- "~/Desktop/Gali_2026_convectionPOC/input_data/globcolour_matchups/"
apath <- "~/Desktop/Gali_2026_convectionPOC/input_data/"      # BGC-Argo profile files
tpath <- "~/Desktop/Gali_2026_convectionPOC/input_data/"      # trajectory files
# opath <- "~/Desktop/Gali_2026_convectionPOC/"                 # output
ppath <- "~/Desktop/Gali_2026_convectionPOC/output/"          # plots

# -----------------------------------------------------------------------------------------------------
# Libraries and functions
library(lubridate)
library(caTools)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(oce)
library(ocedata)
library(fields)

source("~/Desktop/Gali_2026_convectionPOC/f_vertint_binned.R")
source('~/Desktop/Gali_2026_convectionPOC/f_plothovmoller_format_data.R')
source('~/Desktop/Gali_2026_convectionPOC/f_strat_onset_v2.R')
source('~/Desktop/Gali_2026_convectionPOC/f_conv1000_onset.R')
source('~/Desktop/Gali_2026_convectionPOC/f_conv_termin.R')
source('~/Desktop/Gali_2026_convectionPOC/f_K24.R')
source('~/Desktop/Gali_2026_convectionPOC/f_skillstats_xvec_yvec.R')
source('~/Desktop/Gali_2026_convectionPOC/f_myNAstats.R')

# -----------------------------------------------------------------------------------------------------
# Load data (all pre-processed)
# PROFILES AND MATCHUPS
load(paste0(spath,"matchups_CATS-",region,"_L3bin.Rda"))
dfsat <- mm[mm$fwmo==fwmo,]
ifelse(is.null(pocmethod),
       apattern <- paste("Mprof",fwmo,zgrid_out,"from",zgrid_bin,Torig, sep = "_"),
       apattern <- paste("Mprof",pocmethod,fwmo,zgrid_out,"from",zgrid_bin,Torig, sep = "_")) 
afiles <- grep(pattern = apattern,
               x = list.files(path = apath, pattern = ".Rda"),
               value = T)
if (length(afiles)>1) {
  stop("Only one matching file should be found. FWMO=",fwmo)
} else if (length(afiles)==0) {
  stop("No matching file found. FWMO=",fwmo)
}
load(paste0(apath,afiles[1]))

# TRAJECTORY (1000m continuous record: daily and potentially raw too; not monthly because of unwanted behaviour of Tukey-approach spike processing)
tlist <- list.files(path = tpath, pattern = paste0("Mtraj_",fwmo,"_1000m_day_stats_"), full.names = T)
ifelse(docorrchl,
       selectcorrchl <- grepl("CORRCHL", tlist), # STANDARD PROCESSING
       selectcorrchl <- !grepl("CORRCHL", tlist) # PROCESSING WITHOUT CORRECTION FOR EXAMPLE SM FIGURE
)
load(grep("day", tlist[selectcorrchl], value = T)); dbday <- df.bin
idate <- which(grepl("DATE", names(dbday$mean)))
if (length(idate) > 1) {
  idate <- idate[2:length(idate)]
  dbday$mean[,idate] <- NULL
} 
rm(df.bin)

#  BATHY AND MLD DATASETS FOR MAPPING (same code as in fig_conv_events.R)
load("~/Desktop/Gali_2026_convectionPOC/input_data/gebco_08_05degr.Rda") # Load and preprocess bathymetry (NOTE: See preprocess_gebco.R)
data("coastlineWorldMedium"); coastline <- coastlineWorldMedium
load("~/Desktop/Gali_2026_convectionPOC/input_data/subsampled_opera_a67o_v20230630_omlda.omldamax.mlotst.mlotstmax.tos.sos.taum.Rda")

dfa <- as.data.frame(edata)
dimvars <- c("month","year","level","ncell","decdate")
vname <- names(dfa)[!(names(dfa)%in%c("month","year","level","ncell","decdate"))]
vcell <- unique(dfa$ncell)
dfa$decdate <- dfa$year + (dfa$month - 0.5)/12
vtime <- unique(dfa$decdate)
Adata <- unlist(sapply(vname, function(x) return(dfa[[x]]))) # Put data in 3D array format: define array dimensions and populate
A <- array(Adata, dim = c(length(vtime), length(vcell), length(vname)))
Ar <- aperm(A, perm = c(2,3,1))
odim <- dim(Ar) # orig dims
A4 <- array(Ar, dim = c(odim[1],odim[2],12,odim[3]/12))
Amax <- apply(A4, c(1,2,4), max, na.rm=T)

# %%%%%%%%%%%%%%%%%%%%%%%%% COPIED FROM map_profile_counts.NASPG.R %%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load and preprocess bathymetry (NOTE: See preprocess_gebco.R)
if (hgrid=="orca1") {
  load("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/ORCA1_hgrid.Rda")
} else if (hgrid=="lon1lat1") {
  nav_lat <- t(array(data = rep(seq(-89.5,89.5,1), 360), dim = c(180, 360)))
  nav_lon <- array(data = rep(seq(-179.5,179.5,1), 180), dim = c(360, 180))
}

# Function for cropping ORCA1 data for image.plot (can take irregular/non-rectangular grid)----
f_crop4map <- function(xmat, ymat, zmat, zmin, zmax, xindvec, yindvec) {
  xmat <- xmat[ xindvec , yindvec ]
  ymat <- ymat[ xindvec , yindvec ]
  zmat <- zmat[ xindvec , yindvec ]
  zmat[zmat<=zmin] <- NA
  zmat[zmat>zmax] <- zmax # very important when adjusting z scale manually! otherwise values >zmax displayed as NA
  return(list(xmat=xmat, ymat=ymat, zmat=zmat))
}
# Function for interpolating MLD data for contour (cannot take non-rectangular grid: both x and y must be monotonically increasing)----
f_prep4contour1deg <- function(xin, yin, zin, dlondeg, dlatdeg) {
  iinclude <- !(zin==0 | is.na(zin) | abs(zin)==Inf)
  fld <- akima::interp(x = xin[iinclude], y = yin[iinclude], z = zin[iinclude],
                       xo = seq( min(floor(xin)), max(ceiling(xin)), dlondeg),
                       yo = seq( min(floor(yin)), max(ceiling(yin)), dlatdeg),
                       linear = T, extrap = F)
}

# -----------------------------------------------------------------------------------------------------
# Merge matchups into Mprof$prof based on cycle number (more robust than using year-doy because of date format conversions)
if (dim(dfsat)[1]) {
  pdummy <- t(as.data.frame(strsplit(as.character(dfsat$file), split = "[_,.nc]"), col.names = seq(1,dim(dfsat)[1])))
  dfsat$cycle <- as.numeric(pdummy[,3]) # cycle numbers suffixed with D (descending profile) will be turned NA here
  dfsat <- dfsat[!is.na(dfsat$cycle),]  # removing cycle numbers == NA
  mdf <- merge.data.frame(x = Mprof$prof, y = dfsat, by = "cycle", all.x = T, all.y = F)
  rm(pdummy)
} else {
  mdf <- Mprof$prof
  mdf$year <- year(mdf$date)
}
yearS <- unique(floor(mdf$year))
yearS <- yearS[!is.na(yearS)] # NOTE: unnecessary safeguard in the convection paper dataset

# NOTE1: mdf could directly overwrite Mprof$prof, but I keep it as a separate df because it is lighter to handle (Mprof weighs several MB)
# NOTE2: avoid repeated column names other than those used in merging (ie, avoid shared column name "date")

# -----------------------------------------------------------------------------------------------------
# Correct units (sat POC, Argo PAR) and apply spectral conversion factors for bbp
# Calculate 1000m POC with different approaches

if (dim(dfsat)[1]) {
  
  # Satellite POC: from mg to mmol
  jpocsat <- grepl("POC", names(mdf), ignore.case = F)
  mdf[,jpocsat] <- mdf[,jpocsat] / 12
  
  # Satellite BBP: from 443 to 700
  jbbpsat <- grepl("BBP", names(mdf), ignore.case = F)
  etha <- -0.41 # Cetinic 2012. NOTE negative sign
  mdf[,jbbpsat] <- mdf[,jbbpsat] * (700/443)^etha
}

# BGC-Argo PAR: from µmol photons m-2 s-1 to mol photons m-2 d-1. Use index of adjusted PAR if available (jpar)
cfpar <- 86400/1e6
ifelse(sum(!is.na(Mprof$data[,,Mprof$varnameS=="downwelling_par_adjusted"])) > 0,
       jpar <- which(Mprof$varnameS=="downwelling_par_adjusted"),
       jpar <- which(Mprof$varnameS=="downwelling_par"))
Mprof$data[,,jpar] <- Mprof$data[,,jpar] * cfpar
jmdfpar <- grep("par0", names(mdf))
mdf[,jmdfpar] <- mdf[,jmdfpar] * cfpar

# Calculate POC [mmol/m3] with G22 and K24 conversion factors (cf) in different fractions
bbpvarS <- c("BBP700","BBP700_despiked","BBP700_spike","BBP700_bdespiked","BBP700_bspike")
chlvarS <- c("CHLA_ADJUSTED","CHLA_ADJUSTED_despiked","CHLA_ADJUSTED_spike","CHLA_ADJUSTED_bdespiked","CHLA_ADJUSTED_bspike")
cfg22 <- 1000
for (jv in 1:length(bbpvarS)) {
  bv <- bbpvarS[jv]
  cv <- chlvarS[jv]
  dbday$mean[[paste0("POC.G22.",bv)]] <- cfg22 * dbday$mean[[bv]]
  dbday$mean[[paste0("POC.K24.",bv)]] <- f_K24(bbpin = dbday$mean[[bv]], chlin = dbday$mean[[cv]], units="mmolC.m-3")
  dbday$mean[[paste0(cv,"_over_",bv)]] <- dbday$mean[[cv]] / dbday$mean[[bv]]
  dbday$mean[[paste0("CF.K24.",bv)]] <- dbday$mean[[paste0("POC.K24.",bv)]] / dbday$mean[[bv]]
}

# -----------------------------------------------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROCESS PAR IRRADIANCE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -----------------------------------------------------------------------------------------------------
# Check and correct PAR profiles using deepest bin
poffset <- min(Mprof$data[,,jpar], na.rm = T)
if (poffset > -1 & poffset < 1) {
  Mprof$data[,,jpar] <- Mprof$data[,,jpar] - poffset
} else {
  stop("PAR offset exceeds 1 mol photons/m2/d, inspect profiles")
}
# Define default PAR0- variable. Best results so far with par0_exp
jpardef <- grep("par0_exp", names(mdf)) # "par0_dir", "par0_exp", "par0_lin"
par0m <- mdf[,jpardef]

# -----------------------------------------------------------------------------------------------------
# Fit sinusoidal curve to seasonal cycle of the ratio: satellite PAR_0+ / Argo PAR_0-
if (dim(dfsat)[1]) {
  f_doy2rad <- function(x) (x+10)*(pi/365)
  f_sinfit <- function(x,y) {
    jfit <- !is.na(x) & !is.na(y) & y>0 & y<4
    x <- x[jfit]; y <- y[jfit]
    fsin <- nls(formula = y ~ a + b * (sin(x)), model = T, start = list(a = 0, b = 1))
    # print( paste0("R2 of linear fit is ", (cor( y, predict(object = fsin, newdata = x[jfit]) ))^2 ) )
    return(fsin)
  }
  fsin1 <- f_sinfit(f_doy2rad(mdf$doy), mdf$day_1.PAR / par0m)
  fsin2 <- f_sinfit(f_doy2rad(mdf$doy), mdf$day8_5.PAR / par0m)
}
# -----------------------------------------------------------------------------------------------------
# Ad-hoc "correction" of profiles with low PAR (e.g. floats 6901524, 6901527 during 2013)
# These profiles usually acquired at nighttime during productive season (modified schedule to study diel variability?)
# Detect low surface PAR and interpolate using mean of current and [previous,posterior] good profiles
# NOTE: this will smooth out "bad" profiles and potentially spurious features caused by low PAR before further scaling (next step)
# The correction results in much higher correlation with satellite PAR for affected floats
mdf$localtime <- lubridate::local_time(mdf$date, units = "hours", tz = "UTC") + mdf$longitude/15
lowpar <- (mdf$localtime<10 | mdf$localtime>14) & !is.na(mdf$localtime)
jlowpar <- which(lowpar)
jokpar <- which(!lowpar)
for (j in jlowpar) {if (j>1 & j<length(mdf$date)) {
  jprev <- max(c(1,jokpar[jokpar < j]))
  jpost <- min(jokpar[jokpar > j])
  Mprof$data[1,j,jpar] <- (Mprof$data[1,jprev,jpar] + Mprof$data[1,j,jpar] + Mprof$data[1,jpost,jpar])/3
  mdf[j,jmdfpar] <- (mdf[jprev,jmdfpar] + mdf[j,jmdfpar] + mdf[jpost,jmdfpar])/3
}}
par0m <- mdf[,jpardef]

# -----------------------------------------------------------------------------------------------------
# Scale downwelling PAR profile to satellite PAR0. Calculate scaling factor using satellite PAR daily, then 8D, then fill gaps with fsin1 interpolation
# NOTE: not appropriate with low-resolution z grid
if (dim(dfsat)[1]) {
  scalepar <- predict(object = fsin1, newdata = list(x=f_doy2rad( mdf$doy )) )
  mdf[,jmdfpar] <- mdf[,jmdfpar] * scalepar
  # Finally fill gaps, if any, with satellite data
  mdf[is.na(par0m),jmdfpar] <- mdf$day_1.PAR[is.na(par0m)]
  mdf[is.na(par0m),jmdfpar] <- mdf$day8_5.PAR[is.na(par0m)]
  par0m <- mdf[,jpardef]
  SPAR <- matrix(data = scalepar, nrow = length(Mprof$zcenter), ncol = length(scalepar), byrow = T)
  Mprof$data[,,jpar] <- Mprof$data[,,jpar] * SPAR
}

# -----------------------------------------------------------------------------------------------------
# Calculate irradiance exposure metrics
# MLD-mean irradiance, mdf$parmld. Clamp mld to 1000 when sigma-t criterion not met
mld4par <- mdf$tmld_SIGMAT_0.005
mld4par[is.na(mld4par) & mdf$doy<150] <- 1000
mdf$parmld <- par0m * (1 / (mdf$Kdpar_exp * mld4par)) * (1 - exp(-mdf$Kdpar_exp * mld4par))
# Euphotic layer depth based on irradiance threshold, mdf$zeu_corr
# Define bottom of euphotic layer as first level where EdPAR<0.1 mol/m2/d, as per Lacour et al. 2017
mdf$zeu_corr <- sapply(seq(1,dim(mdf)[1]), function(jprof) {
  if (!is.na(par0m[jprof]) & !is.na(mdf$Kdpar_exp[jprof])) {
    z <- 10^(seq(log10(3),log10(300),0.01))
    parz <- par0m[jprof] * exp(-mdf$Kdpar_exp[jprof] * z)
    return( max(z[parz>0.1]) )
  } else {
    return(NA)
  }
})


# -----------------------------------------------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTE PARTICLE SIZE ESTIMATES AND TPOC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -----------------------------------------------------------------------------------------------------
# ESTIMATES BASED ON THE METHOD OF BRIGGS ET AL. 2013 AO AS APPLIED BY REMBAUVILLE ET AL. 2017 JGR
# We take advantage that the method relies on bins means vs STD of detrended data (similar to spike signal, but not equivalent!)
bbpSD <- Mprof$dataSD[,,Mprof$varnameS=="bbp700_detrended"]
bbpMN <- Mprof$data[,,Mprof$varnameS=="bbp700"]
chlSD <- Mprof$dataSD[,,Mprof$varnameS=="chla_adjusted_detrended"]
chlMN <- Mprof$data[,,Mprof$varnameS=="chla_adjusted"]
V <- 10 # mL
Qbb <- 0.024
tstamp <- 1 # s, integration time
tres <- 0.1 # s, residence time
theta <- tres/tstamp
alpha <- ifelse( theta > 1 , 1 - (3*theta)^-1, theta - (theta^2)/3)
Abbp <- (bbpSD/bbpMN)*(10/Qbb)*(1/alpha)
Achl <- (chlSD/chlMN)*(1/alpha)
ESDbbp <- 2*sqrt(Abbp/pi)
SIchl <- 2*sqrt(Achl/pi)

TPOC <- Mprof$data[,,Mprof$varnameS=="spoc"] + Mprof$data[,,Mprof$varnameS=="bpoc"]

Mprof$varnameS <- c( Mprof$varnameS, "ESDbbp" , "SIchl" , "tpoc" )
TMP <- array(NA, dim = c(dim(Mprof$data)[1],dim(Mprof$data)[2],dim(Mprof$data)[3]+3))
TMP[,,1:dim(Mprof$data)[3]] <- Mprof$data
TMP[,,dim(Mprof$data)[3]+1] <- ESDbbp
TMP[,,dim(Mprof$data)[3]+2] <- SIchl
TMP[,,dim(Mprof$data)[3]+3] <- TPOC
Mprof$data <- TMP

# -----------------------------------------------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%% COMPUTE VERTICAL INTEGRALS/MEANS OVER RELEVANT LAYERS %%%%%%%%%%%%%%%%%%%%%%%%
# -----------------------------------------------------------------------------------------------------
svarS <- c("bbp700","bbp700_despiked","spoc","bpoc","tpoc","chla","chla_adjusted","chla_adjusted_spike",
           "chla_spikefreq","bbp700_spike","bbp700_spikefreq","ESDbbp","SIchl")
# NOTE: see layer means (instead of integrals) in "analyze_float_matchups.R"
# TO DO: put loop in external function?
zbounds <- t(data.frame(Mprof$zbreaks[1:(length(Mprof$zbreaks)-1)], Mprof$zbreaks[2:(length(Mprof$zbreaks))] ))
for (vv in svarS) {
  
  # Top 100m integral
  mdf[,paste0(vv,"_","iz0.100")] <- sapply(seq(1,dim(mdf)[1]), function(j) {
    vi <- f_vertint_binned(v = Mprof$data[,j,Mprof$varnameS==vv], # vertical profile, binned onto vertical grid (now zgrid_out)
                           zm = 0,
                           zM = 100,
                           zc = Mprof$zcenter, zb = zbounds, xy_sel = NULL)
    return(vi)
  })
  # Top 200m integral
  mdf[,paste0(vv,"_","iz0.200")] <- sapply(seq(1,dim(mdf)[1]), function(j) {
    vi <- f_vertint_binned(v = Mprof$data[,j,Mprof$varnameS==vv], # vertical profile, binned onto vertical grid (now zgrid_out)
                           zm = 0,
                           zM = 200,
                           zc = Mprof$zcenter, zb = zbounds, xy_sel = NULL)
    return(vi)
  })
  # Mixed layer depth integral (zmax = MLD as defined below)
  mdf[,paste0(vv,"_","iz0.mld")] <- sapply(seq(1,dim(mdf)[1]), function(j) {
    zmax <- min(c(mdf$tmld_SIGMAT_0.03[j], 1000), na.rm = T)
    if (is.na(zmax) | is.null(zmax)) {
      return(NA)
    } else {
      vi <- f_vertint_binned(v = Mprof$data[,j,Mprof$varnameS==vv], # vertical profile, binned onto vertical grid (now zgrid_out)
                             zm = 0,
                             zM = zmax,
                             zc = Mprof$zcenter, zb = zbounds, xy_sel = NULL)
      return(vi)
    }
  })
  # 100-500 m integral
  mdf[,paste0(vv,"_","iz100.500")] <- sapply(seq(1,dim(mdf)[1]), function(j) {
    vi <- f_vertint_binned(v = Mprof$data[,j,Mprof$varnameS==vv], # vertical profile, binned onto vertical grid (now zgrid_out)
                           zm = 100,
                           zM = 500,
                           zc = Mprof$zcenter, zb = zbounds, xy_sel = NULL)
    return(vi)
  })
  # 500-900 m integral
  mdf[,paste0(vv,"_","iz500.900")] <- sapply(seq(1,dim(mdf)[1]), function(j) {
    vi <- f_vertint_binned(v = Mprof$data[,j,Mprof$varnameS==vv], # vertical profile, binned onto vertical grid (now zgrid_out)
                           zm = 500,
                           zM = 900,
                           zc = Mprof$zcenter, zb = zbounds, xy_sel = NULL)
    return(vi)
  })
  # 0-900 m integral
  mdf[,paste0(vv,"_","iz0.900")] <- sapply(seq(1,dim(mdf)[1]), function(j) {
    vi <- f_vertint_binned(v = Mprof$data[,j,Mprof$varnameS==vv], # vertical profile, binned onto vertical grid (now zgrid_out)
                           zm = 0,
                           zM = 900,
                           zc = Mprof$zcenter, zb = zbounds, xy_sel = NULL)
    return(vi)
  })
  # Integral over euphotic layer depth (several criteria available based on either PAR threshold or Kd, similar results)
  mdf[,paste0(vv,"_","iz0.zeu")] <- sapply(seq(1,dim(mdf)[1]), function(j) {
    zmax <- ifelse(!is.na(mdf$Kdpar_exp[j]), 4.6/mdf$Kdpar_exp[j], 4.6/mdf$day_1.KDPAR[j]) # use satellite KdPAR if in situ not available
    if (is.na(zmax) | is.null(zmax)) {
      return(NA)
    } else {
      vi <- f_vertint_binned(v = Mprof$data[,j,Mprof$varnameS==vv], # vertical profile, binned onto vertical grid (now zgrid_out)
                             zm = 0,                                # surface
                             zM = zmax,                             # FOD based on PAR
                             zc = Mprof$zcenter, zb = zbounds, xy_sel = NULL)
      return(vi)
    }
  })
  # Integral of the productive layer "zprod" defined as the deepest of euphotic or mixing layer
  mdf[,paste0(vv,"_","izprod")] <- sapply(seq(1,dim(mdf)[1]), function(j) {
    zmax <- max(4.6/mdf$Kdpar_exp[j], mdf$gmld_CHLA_ADJUSTED[j], na.rm=T)
    if (is.na(zmax) | is.null(zmax)) {
      return(NA)
    } else {
      vi <- f_vertint_binned(v = Mprof$data[,j,Mprof$varnameS==vv], # vertical profile, binned onto vertical grid (now zgrid_out)
                             zm = 0,                                # surface
                             zM = zmax,                             # FOD based on PAR
                             zc = Mprof$zcenter, zb = zbounds, xy_sel = NULL)
      return(vi)
    }
  })
}

# --------------------------------------------------------------------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SEGMENTING YEARS BY MIXING REGIMES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------------------------------------------------------------------------------------------------------------------------
# Segmentation based on combination of criteria for detecting the onset and cessation of deep convection and stable stratification
# NOTE: criteria and metrics tuned for North Atlantic SPG floats. Adaptation needed before applying to e.g. Mediterranean convection
# Correct Chl offset using trajectory file data from stable no-convection periods and estimate LOD
# Compute running mean and SD of daily-binned despiked means
# Use minimum monthly mean for chl baseline
# Use first quartile of monthly running std to calculate noise cutoff (~LOD) as 10*SD
# This cutoff represents the noise of stable ("blank") periods (typically at least 1/4 of the annual cycle)
# Allows robust detection of signal emergence from noise. Can help discard spurious detection of convection onset from Chl
# Histograms of crunsd and crunsm are very illustrative of the robustness of this approach
# Use baseline and noise levels directly in plots rather than "correcting" stored trajectory data
tint <- seq(dbday$mean$DATE[1], dbday$mean$DATE[length(dbday$mean$DATE)], "1 day")
dbintchl <- approx(dbday$mean$DATE, dbday$mean$CHLA_ADJUSTED_despiked, tint, method = "l", rule = 2, ties = mean) # Interpolate to daily, no NA in output
crunsd <- runsd(dbintchl$y, k = 31, endrule=c("NA"), align = c("center"))
crunsm <- runmean(dbintchl$y, k = 31, endrule=c("NA"), align = c("center"))
dbintbbp <- approx(dbday$mean$DATE, dbday$mean$BBP700, tint, method = "l", rule = 2, ties = mean) # Interpolate to daily, no NA in output

# Calculate chl noise threshold and baseline
chlnoise <- quantile(crunsd, .25, na.rm=T) * 10
print(paste0("Chla noise = ",chlnoise,"mg m-3"))
chlbase <- quantile(crunsm, .25, na.rm=T) # do not use min: does not work if negative chl excursions present. First quartile is OK

# Same for bbp (noise only: no such thing as "blanks" in bbp record). Do not attempt drift correction
# TO DO: Use 31-day avergaed bbp, whole signal (not despiked!) to find peak of gravitationally-supplied stock at 1000 m
dbintbbp <- approx(dbday$mean$DATE, dbday$mean$BBP700, tint, method = "l", rule = 2, ties = mean) # Interpolate to daily, no NA in output
brunsd <- runsd(dbintbbp$y, k = 31, endrule=c("NA"), align = c("center"))
brunsm <- runmean(dbintbbp$y, k = 31, endrule=c("NA"), align = c("center"))
bbpnoise <- quantile(brunsd, .25, na.rm=T) * 10

# -----------------------------------------------------------------------------------------------------
# Define deep convection and stratification periods according to diverse criteria
# Deep convection onset based on trajectory file Chl (could also be based on first annual exceedance of a given MLD threshold)
cdates <- as_datetime( sapply(yearS, function(x) {
  f_conv1000_onset(bo = dbday$mean$CHLA_ADJUSTED_despiked[which(year(dbday$mean$DATE)==x)],
                   bodate = dbday$mean$DATE[which(year(dbday$mean$DATE)==x)],
                   bocut = chlnoise)
}) )
# Deep convection onset based on trajectory file bbp
bdates <- as_datetime( sapply(yearS, function(x) {
  f_conv1000_onset(bo = dbday$mean$BBP700_despiked[which(year(dbday$mean$DATE)==x)],
                   bodate = dbday$mean$DATE[which(year(dbday$mean$DATE)==x)],
                   bocut = bbpnoise*2/3)
}) )
# Convection termination based on latest detection of MLD0.03>1000 (NA)
# Initially used 0.005 and 0.01 sigma-t thresholds
edates <- as_datetime( sapply(yearS, function(x) {
  yi <- which(mdf$year==x)
  f_conv_termin(zmix = mdf[yi,"tmld_SIGMAT_0.03"], zmixdate = mdf$date[yi])
}) )
# Stratification onset based on MLD time series and euphotic depth criterion
# Initially used 0.005 and 0.01 sigma-t thresholds
zeucut <- max(4.6/mdf$Kdpar_exp, na.rm=T) # permanent stratification occurs when 31-day smoothed MLD crosses the deepest annual euphotic layer depth
sdates <- as_datetime( sapply(yearS, function(x) {
  yi <- which(mdf$year==x)
  f_strat_onset_v2(zmix = mdf[yi,"tmld_SIGMAT_0.005"], zmixdate = mdf$date[yi], zcut = zeucut)
}) )

# -----------------------------------------------------------------------------------------------------
# Define convective and non-convective periods based on consensus between different initiation-termination criteria
convperiod <- data.frame(cbind(cdates, bdates, edates, sdates))
convperiod$startdates <- sapply(1:length(cdates), function(x) {
  dvec <- convperiod[x,c("cdates","bdates")]
  ifelse(sum(is.na(dvec))<2, return(max(dvec, na.rm=T)), return(NA))
})
convperiod$enddates <- edates
convperiod.a <- lapply(convperiod, date) # Need to re-format dates prior to export. #View(as.data.frame(convperiod.a))

# Export csv, edit manually the file for each float USING TEXT EDITOR, then load *_manu*csv file and overwrite start and end dates
# Currently writing these files to subdirectory "test_preprocessing" to avoid overwriting correct files used later in the workflow
# Editing includes deciding whether there is event or not, in which case both start and end dates should be NA
write.csv(convperiod.a, file = paste0(apath,test_preprocessing,"/convperiod_",fwmo,"_auto.csv"), row.names = F)
if (!file.exists(paste0(apath,"convperiod_",fwmo,"_manu.csv"))) {
  convperiod.m <- convperiod.a
  write.csv(convperiod.m, file = paste0(apath,test_preprocessing,"/convperiod_",fwmo,"_manu.csv"), row.names = F)
} else {
  convperiod.m <- read.csv(file = paste0(apath,"v0_convperiod_",fwmo,"_manu.csv"))
}
convperiod.m <- as.data.frame(lapply(convperiod.m, function(x) {
  x <- force_tz(as.POSIXct(as.Date(x), tzone = "GMT")) # View(as.data.frame(convperiod.m))
}))

# -----------------------------------------------------------------------------------------------------
# Calculate the statistics for convective vs. non-c periods => ystats data frame (see NOTE on alternative coding solution without for loop)
# In parallel, simply append cperiods to dbyy to be able to make boxplots by variables, years and periods
row.names(convperiod.m) <- as.character(yearS)
d1 <- as.duration("1 day")
ystats <- list()
dbdp <- list()

# Select only 12 "official" events
convperiods <- read.csv("~/Desktop/Gali_2026_convectionPOC/input_data/convperiods_all.csv")
convperiod.m <- convperiod.m[year(convperiod.m$startdates) %in% convperiods$year[convperiods$fwmo==fwmo],]

for (yy in yearS) {
  
  cp <- convperiod.m[ as.character(yy),] # subset by year
  
  if (!is.na(cp$startdates) & !is.na(cp$enddates)) {
    
    # Define vectors of pre-, conv and bloom (post-conv) periods
    dbyy <- dbday$mean[year(dbday$mean$DATE)==yy,]
    dbyymed <- dbday$med[year(dbday$med$DATE)==yy,]
    cperiods <- cut(x = dbyy$DATE,
                    breaks = c(as.POSIXct(paste0(yy,"-01-01")), cp$startdates-d1, cp$enddates, dbyy$DATE[length(dbyy$DATE)]),
                    labels = c("pre","conv","post"))
    
    # Ad-hoc correction to very extreme daily mean values of bbp700 (over 995% quantile of the bbp700 mean/median ratio): replace by median
    rbbp <- dbyy$BBP700/dbyymed$BBP700
    rout <- quantile(rbbp, 0.995, na.rm=T)
    iout <- which(rbbp > rout)
    dbyy$BBP700[iout] <- dbyymed$BBP700[iout]
    
    # Compute period statistics: mean, sd, day counts and summary statistics (so that mean and maximum event amplitude can be computed)
    ystats[[as.character(yy)]] <-
      lapply(dplyr::select(dbyy, -c(DATE,JULD,CYCLE_NUMBER,doy)), function(yyvv) {
        iok <- !is.na(yyvv) & !is.na(cperiods)
        pstats <- merge(
          array2DF(tapply(yyvv, cperiods, function(x) {
            data.frame(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T))    # mean and sd functions
          })),
          array2DF(tapply(yyvv[iok], cperiods[iok], summary)),      # summary function
          by = "Var1")
        porder <- factor(pstats$Var1, ordered = T, levels = unique(cperiods[iok]))
        pstats <- pstats[order(porder),]                            # key to get matching period stats and observed day counts (below)
        pstats$counts <- table(cperiods)                            # observed days per period
        return(pstats)
      })
    dbdp[[as.character(yy)]] <- cbind(dbyy, cperiods)
  }
}
if (length(ystats)) {
  
  # Convert list of lists (yearS, varnames) to data frame for export
  Ystats <- lapply(ystats, function(l) {
    # use "fill=T" in case data frames in the list have different columns (which shouldn't occur after removing NA's prior to calling the summary function)
    data.table::rbindlist(l, use.names = T, fill = T, idcol = "varname")
  })
  ystats <- data.table::rbindlist(Ystats, use.names = T, fill = F, idcol = "year")
  ystats <- dplyr::rename(ystats, cperiod=Var1)
  dbdp <- data.table::rbindlist(dbdp, use.names = T, fill = F, idcol = "year")
  dbdp <- dbdp[!is.na(dbdp$cperiods),]
}

# -----------------------------------------------------------------------------------------------------
# GENERAL PLOT SETTINGS
# Color settings
ncol <- 100
ncol.chl <- ncol # colors, palettes
col.chl <- oce.colorsChlorophyll(ncol.chl)[1:(ncol.chl*0.8)]
col.chltraj <- col.chl[20] # similar to "yellowgreen"

pal.chl1000 <- colorRampPalette(c("black","green")) # custom1
col.chl1000 <- pal.chl1000(ncol.chl)

pch.years <- c(0,15,1,16,2,17)[1:length(yearS)]
names(pch.years) <- yearS

ncol.poc <- ncol
pal.poc <- colorRampPalette(c("beige","wheat3","cadetblue4","mediumpurple4")) # custom1
col.poc <- pal.poc(n = ncol.poc)
col.poctraj <- col.poc[30] # similar to "tan"

col.zeu <- "magenta2"
col.tmld03 <- "black"
col.tmld01 <- "black"
col.tmld005 <- "black"
col.gmldchl <- "magenta2"
col.MLDtemp <- "gray"

col.month <- brewer.pal(12, "Paired")
col.bathy <- gray.colors(100)[1:75]

# Additional axis settings
dateticks1 <- seq(floor_date(min(mdf$date), unit = "3 months"), floor_date(max(mdf$date), unit = "3 months"), by = "3 months")
dateticks2 <- seq(floor_date(min(mdf$date), unit = "months"), floor_date(max(mdf$date), unit = "months"), by = "months")

# -----------------------------------------------------------------------------------------------------
# Center data on convection start date for each year and event (if "exporteventyear" enabled")
if (savedata) {
  if (exporteventyear %in% (seq(2013,2024))) {
    
    # Get start date of convection event
    startdate <- (convperiod.m[ as.character(exporteventyear),"startdates"]) # subset by year
    # Calculate time shift to center data on convection start date 1 February
    shiftdays <- as.duration( startdate - force_tz(as.POSIXct(paste0(exporteventyear,"-02-01"), tzone = "GMT")) )
    # Apply time shift to all time variables
    yearorigin <- as.POSIXct(paste0(exporteventyear-1,"-12-31"))
    mdf$date <- mdf$date - shiftdays
    mdf$doy <- floor(julian(x = mdf$date, origin = yearorigin))
    dbdp$DATE <- dbdp$DATE - shiftdays
    dbdp$doy <- floor(julian(x = dbdp$DATE, origin = yearorigin))
    lapply(dbday, function(dlist) {
      if ("DATE" %in% names(dlist)) {dlist$DATE <- dlist$DATE - shiftdays}
      # NOTE: do not shift "doy" nor "JULD" because these variables are not used to plot elements of dbday list
    })
  }
}

# -----------------------------------------------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%% TIME BINNING OF MERGED DATA FRAME (mdf) AND TRAJ (dbday) %%%%%%%%%%%%%%%%%%%%%%
# -----------------------------------------------------------------------------------------------------
# Code section inherited from mergeByFloat_*R

# Store non-binned data in ".orig" variables
Mprof.orig <- Mprof
mdf.orig <- mdf

if (tgrid != "Torig") {
  
  # GENERAL VARIABLES
  # -----------------
  ndays <- as.numeric(as.duration(tgrid))/86400
  
  # PROFILE DATA WITH SATELLITE MATCHUPS
  # ------------------------------------
  # Create vector b4vec for binning data frame columns, and df b4mat for 2D binning
  # Use DOY for binning (see note 4 at bottom of mergeByFloat_binZ_selected_floats.R)
  ddate <- decimal_date(mdf$date) - year(mdf$date[1])
  b4vec <- floor(ddate*365/ndays)*ndays
  b4mat <- data.frame(DATEBIN = sort(rep(b4vec, length(Mprof$zcenter))), DEPTH = rep(Mprof$zcenter, length(mdf$date)))
  
  # Temporal binning of profile time series data frame
  # NOTE: To get CYCLE as integer, take only the min of the different cycles included in a given time bin
  mdf.minsVars <- aggregate(subset.data.frame(mdf, select = c(cycle, date, doy)), by = list(DATEBIN = b4vec), function(x) min(x, na.rm = T), simplify = TRUE)
  mdf.meansVars <- aggregate(subset.data.frame(mdf, select = -c(cycle, date, doy)), by = list(DATEBIN = b4vec), function(x) mean(x, na.rm = T), simplify = TRUE)
  
  # Temporal binning of merged observations data frame in long format
  # Format time-binned data as 3D arrays with implicit dimensions: depth x time x variable, plus summarized profile data
  for (dd in c("data","dataSD","dataN")) {
    arrd <- dim(Mprof[[dd]])
    TOBIN <- array(Mprof[[dd]], dim = c(arrd[1]*arrd[2], arrd[3])) # re-arrange in long format prior to binning
    TMP <- aggregate(TOBIN, by = list(DEPTH = b4mat$DEPTH, DATEBIN = b4mat$DATEBIN), function(x) nanmean(x), simplify = TRUE)
    TMP <- subset.data.frame(TMP, select = -c(DEPTH, DATEBIN))
    Mprof[[dd]] <- array(unlist(TMP), dim = c(length(Mprof$zcenter), dim(mdf.minsVars)[1], length(Mprof$varnameS)) )
  }
  mdf.minsVars$DECDATE <- year(mdf$date[1]) + mdf.minsVars$DATEBIN/365
  mdf <- cbind.data.frame(mdf.minsVars, mdf.meansVars[,2:dim(mdf.meansVars)[2]]) # keep DATEBIN for posterior gap-filling or interpolation
  
  # Force time stamps (doy and date) to "tgrid" intervals (not if exporteventyear enabled: causes error)
  if (exporteventyear %in% (seq(2013,2024))) {mdf$doy <- mdf$doy - (mdf$doy%%ndays) + ceiling(ndays/2)}
  mdf$date <- date_decimal(mdf$DECDATE) + as.duration(tgrid)/2
  
  # Remove Mprof$prof: no longer useful. Mprof.orig$prof holds original profile list, and mdf files hold the same merged with extra-data
  Mprof$prof <- NULL
  
  # TRAJECTORY DATA WITH DEFAULT DAILY RESOLUTION: SELECT LIST ELEMENTS TO BE BINNED
  # --------------------------------------------------------------------------------
  ddate <- decimal_date(dbday$mean$DATE) - year(dbday$mean$DATE[1])
  b4traj <- floor(ddate*365/ndays)*ndays
  dbbin <- list()
  tvarS <- c("mean","med")
  for (tt in tvarS) {
    db.minsVars <- aggregate(subset.data.frame(dbday[[tt]], select = c(CYCLE_NUMBER, DATE, JULD)), by = list(DATEBIN = b4traj), function(x) min(x, na.rm = T), simplify = TRUE)
    db.meansVars <- aggregate(subset.data.frame(dbday[[tt]], select = -c(CYCLE_NUMBER, DATE, JULD)), by = list(DATEBIN = b4traj), function(x) mean(x, na.rm = T), simplify = TRUE)
    db.minsVars$DECDATE <- year(dbday[[tt]][["DATE"]][1]) + db.minsVars$DATEBIN/365
    dbbin[[tt]] <- cbind.data.frame(db.minsVars, db.meansVars[,2:dim(db.meansVars)[2]])
  }
  # Force time stamp to "tgrid" intervals
  dbbin[[tt]][["DATE"]] <- date_decimal(dbbin[[tt]][["DECDATE"]]) + as.duration(tgrid)/2
}

# -----------------------------------------------------------------------------------------------------
# Regularize time-binned data, filling unobserved time bins with NA in both mdf and Mprof$data* arrays
if (regularize) {
  obsbins <- mdf$DATEBIN
  totbins <- seq(min(mdf$DATEBIN), max(mdf$DATEBIN), ndays)
  y1 <- min(mdf$year, na.rm=T)
  if (length(obsbins)!=length(totbins)) {
    tmatch <- which(totbins%in%obsbins)
    tmpdf <- as.data.frame(matrix( data = NA, nrow = length(totbins), ncol = dim(mdf)[2]) )
    names(tmpdf) <- names(mdf)
    tmpdf[ tmatch, ] <- mdf
    tmpdf$DATEBIN <- totbins
    tmpdf$DECDATE <- min(mdf$year, na.rm=T) + tmpdf$DATEBIN/365
    tmpdf$date <- date_decimal(tmpdf$DECDATE) + as.duration(tgrid)/2
    tmpdf$doy <- tmpdf$DATEBIN - 365*(floor(tmpdf$DECDATE) - min(mdf$year, na.rm=T)) + ceiling(ndays/2)
    tmpdf$year <- year(tmpdf$date)
    mdf <- tmpdf; rm(tmpdf)
    for (dd in c("data","dataSD","dataN")) {
      TMPARR <- array( data = NA, dim = c( dim(Mprof[[dd]])[1], length(totbins), dim(Mprof[[dd]])[3]) )
      TMPARR[,tmatch,] <- Mprof[[dd]]
      Mprof[[dd]] <- TMPARR; rm(TMPARR)
    }
  }
}

# -----------------------------------------------------------------------------------------------------
# Write out all the data by float and event
if (savedata) {
  ifelse(regularize, reg <- "reg", reg <- "noreg")
  ifelse(exporteventyear %in% (seq(2013,2024)), Tshift <- "Tshift", Tshift <- "noTshift")
  ifelse(exporteventyear %in% (seq(2013,2024)), Tshift <- "Tshift", Tshift <- "noTshift")
  ifelse(
    docorrchl,
    fname <- paste("~/Desktop/Gali_2026_convectionPOC/input_data/test_preprocessing/data",fwmo,Tshift,exporteventyear,paste0(reg,".Rda"), sep = "_"),
    fname <- paste("~/Desktop/Gali_2026_convectionPOC/input_data/test_preprocessing/data",fwmo,Tshift,exporteventyear,paste0(reg,".NOCORRCHL.Rda"), sep = "_")
  )
  # if (!file.exists(fname)) {
    save(dbdp, ystats, mdf, Mprof, mdf.orig, Mprof.orig, file = fname)
  # }
}

# stop("AAA")
# -----------------------------------------------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%% MAKE PLOTS INTEGRATING PROFILES, TRAJECTORY AND SATELLITE %%%%%%%%%%%%%%%%%%%%%
# -----------------------------------------------------------------------------------------------------

# Plot parameters and map region
selperiod <- force_tz(as.POSIXct(c(paste0(min(yearS),"-01-01"),paste0(max(yearS),"-12-31"))), tzone = "GMT")
pyearS <- yearS
if (fwmo==6901524) {
  selperiod <- force_tz(as.POSIXct(c(paste0(2014,"-01-01"),paste0(max(yearS),"-12-31"))), tzone = "GMT")  # remove year 2013 (instrumental issues?)
  pyearS <- yearS[2:(length(yearS))]
}
if (fwmo==6901472) {
  selperiod <- force_tz(as.POSIXct(c(paste0(2013,"-01-01"),paste0(max(yearS),"-12-31"))), tzone = "GMT")  # remove year 2013 (instrumental issues?)
  pyearS <- yearS[2:(length(yearS))]
}
pw <- 40*( length(pyearS)/5 )
ph <- 23
if (fwmo %in% c(6901480,6901481,6901482,6901484,6901485,6901486,6901489,6901521,6901523,6901524,6901527)) {
  rlon <- c(-62,-23)
  rlat <- c(50,66)
} else if (fwmo==6901516) { # if automatic, -48 -10 and 46 68
  rlon <- c(-50,-10)
  rlat <- c(50,65)
} else if (fwmo==6901472) { # if automatic, -54 -23 and 13 28
  rlon <- c(-75,-5)
  rlat <- c(2,38)
} else {
  rlon <- c(floor(min(mdf$longitude, na.rm=T)-10),ceiling(max(mdf$longitude, na.rm=T)+10))
  rlat <- c(floor(min(mdf$latitude, na.rm=T)-5),ceiling(max(mdf$latitude, na.rm=T)+5))
  if (fwmo %in% c(6901516,6901472)) print(rlon); print(rlat)
}

# Z variable ranges for Hovmoller plots
zr.chl <- 10^c(-2,1)		# Epi. Max is 1 µg L-1 for NASTG and 10 for NASPG
zr.chls <- 10^c(-3,0)		# chl_spike
zr.chl <- 10^c(-2.5,1); zr.chls <- zr.chl
zr.bbp <- 10^c(-4,-2)     # bbp700_despiked			
zr.bbps <- 10^c(-5,-3)    # bbp700_spike
zr.bbp <- 10^c(-4.5,-2); zr.bbps <- zr.bbp
zr.spoc <- 10^c(-1,1)			# Epi Max is 10 µM for NASTG and 100 for NASPG
zr.bpoc <- 10^c(-2,1)

# Custom function to get current plot coordinates and store in a variable for later use (eg, overlay a custom-made colorbar)
f_get_pcoords <- function() {
  u <- par("usr") 
  v <- c(grconvertX(u[1:2], "user", "ndc"), grconvertY(u[3:4], "user", "ndc"))
  return(v)
}

# --------------------------------------------------------------------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMBINED HOVMOELLER AND TRAJ PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------------------------------------------------------------------------------------------------------------------------
if (plotsummary) {
  p <- paste0( "~/Desktop/Gali_2026_convectionPOC/output/Fig_S3-S9_", fwmo , "_", paste(selperiod, collapse = "_"),"_",biospike,".png" )
  png(filename = p, width = pw, height = ph, units = 'cm', pointsize = 8, bg = "white", res = 600, type = "quartz")
  
  
  # Multipanel setup
  m1 <- list()
  yyp1 <- length(pyearS)+0 # decide here whether to add (+1) or not (+0) an extra panel for color bar (months) in the top row
  pw1res <- pw%%yyp1
  pw1 <- (pw-pw1res)/yyp1
  for (jy in 1:(yyp1+1)) {
    ifelse( jy==(yyp1+1), m1[[jy]] <- matrix(data = jy, nrow = 4, ncol = pw1res), m1[[jy]] <- matrix(data = jy, nrow = 4, ncol = pw1))
  }
  m1 <- Reduce(cbind, m1)
  m2 <- matrix(data = max(m1)+1, nrow = 3, ncol = pw)
  m3 <- m2+1
  m4 <- m2+2
  m5 <- m2+3
  m6 <- matrix(data = m5+1, nrow = 3, ncol = pw)
  m7 <- m6+1
  m8 <- matrix(data = m7+1, nrow = 3, ncol = pw)
  layout(rbind(m1, m2, m3, m4, m5, m6, m7, m8)) # layout.show(seq(1,10))
  par(oma = c(1,2,0,1.5))
  
  # -----------------------------------------------------------------------------------------------------
  # 1ST ROW: Maps (as many as pyearS) of float trajectory overlaid on bathymetry. Dots colored by month. Contour of maximum March MLD (reanalysis?)
  gebco$zmat[gebco$zmat>0] <- NA
  gebco$zmat[gebco$zmat<(-8000)] <- NA
  Amlotstmax <- Amax[,which(vname=="mlotstmax.m"),unique(dfa$year)%in%pyearS]
  
  for (jy in 1:length(pyearS)) {
    
    yearmonth <- pyearS[jy]+2.5/12
    ymdf <- mdf.orig[ mdf.orig$year==pyearS[jy] , ]
    
    ifelse( yyp1 > length(pyearS), par(mar = c(1,4,1,0)), par(mar = c(1,4,1,6)) ) 
    image(gebco$lonvec, gebco$latvec, gebco$zmat, xlab = "", ylab = "", xlim = rlon, ylim = rlat, col = col.bathy,
          xaxt = "n", yaxt = "n", main = "")
    plot(coastline, clon = 0, clat = 0, span = c(length(gebco$lon), length(gebco$lon)),
         col = "black", bg = rgb(1,1,1,0), lwd = 0.5, axes = F, add = T)
    contour(gebco$lonvec,gebco$latvec,gebco$zmat, levels = c(-1000,-2000), lwd = 0.7, lty = 1, labcex = 1, col = "cornsilk4", drawlabels = F, add = T)
    zoomm <- f_prep4contour1deg(xin = nav_lon[vcell], yin = nav_lat[vcell],
                                zin = Amlotstmax[,jy],                                                 # Annual maximum MLD contours
                                # zin = A[which(vtime==yearmonth),,which(vname=="mlotstmax.m")],       # March MLD contours
                                dlondeg = 0.5, dlatdeg = 0.5)
    contour(zoomm$x, zoomm$y, zoomm$z, levels = c(500,1000), lwd = c(1,2), lty = 1, cex = 1, col = "black", drawlabels = F, add = T)
    points(ymdf$lon, ymdf$lat, col = col.month[ month(ymdf$date) ], pch = 20, cex = 2.5)
    title(main = pyearS[jy], line = -2.5, cex.main = 3.5, col.main = "white", adj = 0.1) # cex.main = 3
    
  }
  if(yyp1 > length(pyearS)) {plot.new()} # additional year (added through yyp1 to make room for color bar)
  if (pw1res!=0) plot.new()
  
  # -----------------------------------------------------------------------------------------------------
  # Hovmoeller diagrams
  
  # 2ND ROW: CHLA
  p2 <- f_plothovmoller_format_data(L = list(date=mdf$date, depth = Mprof$zcenter, chla_adjusted_despiked = Mprof$data[,,Mprof$varnameS=="chla_adjusted_despiked"]),
                                    xn = "date", yn = "depth", zn = "chla_adjusted_despiked",
                                    xr = selperiod, yr = c(5,1000), zr = zr.chl, zlog = T)
  par(mar = c(0,5,0,5 ))
  image(x = p2$x, y = p2$y, z = p2$z, log = "y",
        xlim = p2$xlim, ylim = rev(p2$ylim), zlim = p2$zlim, col = col.chl,
        xaxt = "n", xlab = "", ylab = "", cex.axis = 2) # cex.axis = 1.5
  # Alternatives to zeu: zeu_dir, zeu_corr (BEST), , 4.6/Kdpar_exp (second best), satellite (eg day8_5.ZEU), 4.6/Kdpar_lin
  lines(mdf.orig$date, mdf.orig[["zeu_corr"]], lwd = 2, col = col.zeu)
  # lines(mdf.orig$date, mdf.orig[["gmld_CHLA_ADJUSTED"]], lwd = 2, col = col.gmldchl, lty = 3)
  # lines(mdf.orig$date, mdf.orig[["tmld_SIGMAT_0.03"]], lwd = 1, col = col.tmld03)
  lines(mdf.orig$date, mdf.orig[["tmld_CTEMP_0.03"]], lwd = 2, col = col.MLDtemp) # dsigmat=0.005 at 1000 m in Lab Sea corresponds to dT≈0.024
  lines(mdf.orig$date, mdf.orig[["tmld_SIGMAT_0.005"]], lwd = 1, col = col.tmld005, lty = 3)
  abline(v = dateticks1, lty = 3, lwd = 1, col = "darkgray")
  abline(h = c(20,100), lty = 3, lwd = 1, col = "darkgray")
  mtext(text = "Depth (m)", side = 2, line = 3.2, cex = 2) # cex = 1.3
  vp2 <- f_get_pcoords() # get main plot coordinates to overlay colorbar
  
  # 3RD ROW: CHLA spikes
  p3 <- f_plothovmoller_format_data(L = list(date=mdf$date, depth = Mprof$zcenter,
                                             chla_adjusted_spike = Mprof$data[,,Mprof$varnameS=="chla_adjusted_spike"]),
                                    xn = "date", yn = "depth", zn = "chla_adjusted_spike",
                                    xr = selperiod, yr = c(5,850), zr = zr.chls, zlog = T)
  par(mar = c(0,5,0,5))
  image(x = p3$x, y = p3$y, z = p3$z, log = "y",
        xlim = p3$xlim, ylim = c(1000,5), zlim = p3$zlim, col = col.chl,
        xaxt = "n", xlab = "", yaxt = "n", ylab = "", cex.axis = 10)
  # axis(side = 4, cex.axis = 2)
  lines(mdf.orig$date, mdf.orig[["zeu_corr"]], lwd = 2, col = col.zeu)
  # lines(mdf.orig$date, mdf.orig[["gmld_CHLA_ADJUSTED"]], lwd = 2, col = col.gmldchl, lty = 3)
  # lines(mdf.orig$date, mdf.orig[["tmld_SIGMAT_0.03"]], lwd = 1, col = col.tmld03)
  lines(mdf.orig$date, mdf.orig[["tmld_CTEMP_0.03"]], lwd = 2, col = col.MLDtemp) # dsigmat=0.005 at 1000 m in Lab Sea corresponds to dT≈0.024
  lines(mdf.orig$date, mdf.orig[["tmld_SIGMAT_0.005"]], lwd = 1, col = col.tmld005, lty = 3)
  abline(v = dateticks1, lty = 3, lwd = 1, col = "darkgray")
  abline(h = c(20,100), lty = 3, lwd = 1, col = "darkgray")
  # mtext(text = "Depth (m)", side = 4, line = 4, cex = 2)
  legend("bottomright", cex = 2.0,
         legend = c(expression(paste(Z[eu])),
                    # expression(MLD[fluo]),
                    expression(MLD[0.005]),
                    expression(MLD[Temp0.03])),
         col = c(col.zeu, col.tmld005, col.MLDtemp), lwd = c(2,1,2), lty = c(1,3,1), bg = "#FFFFFF95", box.col = "#FFFFFF00")
  
  # 4TH ROW: bbp despiked
  p4 <- f_plothovmoller_format_data(L = list(date=mdf$date, depth = Mprof$zcenter, bbp700_despiked = Mprof$data[,,Mprof$varnameS=="bbp700_despiked"]),
                                    xn = "date", yn = "depth", zn = "bbp700_despiked",
                                    xr = selperiod, yr = c(5,1000), zr = zr.bbp, zlog = T)
  par(mar = c(0,5,0,5))
  image(x = p4$x, y = p4$y, z = p4$z, log = "y",
        xlim = p4$xlim, ylim = rev(p4$ylim), zlim = p4$zlim, col = col.poc,
        xaxt = "n", xlab = "", ylab = "", cex.axis = 2)
  lines(mdf.orig$date, mdf.orig[["tmld_SIGMAT_0.03"]], lwd = 1, col = col.tmld03)
  lines(mdf.orig$date, mdf.orig[["tmld_SIGMAT_0.01"]], lwd = 1, col = col.tmld01, lty = 2)
  lines(mdf.orig$date, mdf.orig[["tmld_SIGMAT_0.005"]], lwd = 1, col = col.tmld005, lty = 3)
  abline(v = dateticks1, lty = 3, lwd = 1, col = "darkgray")
  abline(h = c(20,100), lty = 3, lwd = 1, col = "darkgray")
  mtext(text = "Depth (m)", side = 2, line = 3.2, cex = 2)
  vp4 <- f_get_pcoords() # get main plot coordinates to overlay colorbar
  
  # 5TH ROW: bbp spike
  p5 <- f_plothovmoller_format_data(L = list(date=mdf$date, depth = Mprof$zcenter, bbp700_spike = Mprof$data[,,Mprof$varnameS=="bbp700_spike"]),
                                    xn = "date", yn = "depth", zn = "bbp700_spike",
                                    xr = selperiod, yr = c(5,850), zr = zr.bbps, zlog = T)
  par(mar = c(0,5,0,5))
  image(x = p5$x, y = p5$y, z = p5$z, log = "y",
        xlim = p5$xlim, ylim = c(1000,5), zlim = p5$zlim, col = col.poc,
        xaxt = "n", xlab = "", yaxt = "n", ylab = "", cex.axis = 2)
  # axis(side = 4, cex.axis = 2, at = c(5,10,20,50,100,200,500,1000), labels = rep("",8))
  lines(mdf.orig$date, mdf.orig[["tmld_SIGMAT_0.03"]], lwd = 1, col = col.tmld03)
  lines(mdf.orig$date, mdf.orig[["tmld_SIGMAT_0.01"]], lwd = 1, col = col.tmld01, lty = 2)
  lines(mdf.orig$date, mdf.orig[["tmld_SIGMAT_0.005"]], lwd = 1, col = col.tmld005, lty = 3)
  abline(v = dateticks1, lty = 3, lwd = 1, col = "darkgray")
  abline(h = c(20,100), lty = 3, lwd = 1, col = "darkgray")
  # mtext(text = "Depth (m)", side = 4, line = 4, cex = 2)
  legend("bottomright",
         legend = c(expression(MLD[0.005]),
                    expression(MLD[0.01]),
                    expression(MLD[0.03])),
         col = c(col.tmld005, col.tmld03), lwd = c(1,1,1), lty = c(3,2,1), bg = "#FFFFFF95", box.col = "#FFFFFF00", cex = 2)
  
  # -----------------------------------------------------------------------------------------------------
  # Trajectory 1000 m plots
  
  # 6TH ROW: Trajectory data for CHLA
  par(mar = c(0,5,0,5))
  ytot <- dbday$mean$CHLA_ADJUSTED - chlbase # chlbase is very small and used only for visual adjustment of mean detrended Fchl on 0
  plot(dbday$mean$DATE, ytot, type="l", col="gray", lwd=0.1, ylog = "T", cex.axis = 2, bty = "n", bg = rgb(1,1,1,1),
       ylim = c(0, 0.12), #ylim = c(0.020, 0.12),
       # ylim=c(0.99, 1.1)*quantile(ytot, c(0.001,0.999), na.rm=T),
       xlim=p2$xlim, xaxt = "n", xlab = "", ylab = "", xaxs="i")
  lines(dbday$mean$DATE, dbday$mean[[paste0("CHLA_ADJUSTED_",biospike,"despiked")]] - chlbase, lwd=2, col=col.chltraj)
  lines(dbday$mean$DATE, ytot, lwd=0.5, col="black")
  abline(v = dateticks1, lty = 3, lwd = 1, col = "darkgray")
  points(x = convperiod.m$startdates, y = rep(0.95*par("usr")[4], length(convperiod.m$startdates)), pch=6, col="blue", lwd=4, ljoin=1)
  points(x = convperiod.m$enddates, y = rep(0.95*par("usr")[4], length(convperiod.m$enddates)), pch=2, col="blue", lwd=4, ljoin=1)
  abline(v = c(convperiod.m$startdates,convperiod.m$enddates), lty = 2, lwd = 1, col = "blue")
  abline(h = chlnoise, lty = 2, lwd = 1, col = rgb(0,1,0,0.8))
  mtext(text = expression(paste(Delta,"FChl",italic("a")," (mg ",m^-3,")")), side = 2, line = 3.2, cex = 2) # cex = 1.3
  legend(ifelse(fwmo==6901524,"topright","topleft"), cex = 2.0,
         legend = c("daily","daily (despiked)","noise threshold"),
         col = c("black",col.chltraj,rgb(0,1,0,0.8)),
         lwd = c(1,2,1),
         lty = c(1,1,2),
         seg.len = c(1,1,1),
         bg = "#FFFFFF95", box.col = "#FFFFFF00")
  
  
  # 7TH ROW: Trajectory data for BBP700
  par(mar = c(0,5,0,5))
  ytot <- dbday$mean$BBP700 * 1e3
  plot(dbday$mean$DATE, ytot, type="l", col="gray", lwd=0.1, ylog = "F", cex.axis = 1.5, bty = "n", bg = rgb(1,1,1,1),
       ylim = c(0.00,0.50), #ylim = c(0.00015,0.00035),
       # ylim=c(0.99, 1.05)*range(ytot, na.rm=T),
       xlim=p2$xlim, xaxt = "n", yaxt = "n", xlab = "", ylab = "", xaxs="i")
  axis(side = 4, cex.axis = 2, padj = 0.5) # , at = seq(0,0.0005,0.0001), labels = c(0,"",0.2,"",0.4,"")
  lines(dbday$mean$DATE, dbday$mean[[paste0("BBP700_",biospike,"despiked")]] * 1e3, lwd=2, col=col.poctraj)
  lines(dbday$mean$DATE, ytot, lwd=0.5, col="black")
  abline(v = dateticks1, lty = 3, lwd = 1, col = "darkgray")
  abline(v = c(convperiod.m$startdates,convperiod.m$enddates), lty = 2, lwd = 1, col = "blue")
  mtext(text = expression(paste("b"[bp700]," (",km^-1,")")), side = 4, line = 5, cex = 2)
  legend("topright", cex = 2.0, 
         legend = c("daily","daily (despiked)"),
         col = c("black",col.poctraj),
         lwd = c(1,2),
         lty = c(1,1),
         seg.len = c(1,1),
         bg = "#FFFFFF95", box.col = "#FFFFFF00")
  
  # 8TH ROW: Trajectory data for TEMP
  par(mar = c(4,5,0,5))
  ytot <- dbday$mean$CTEMP
  plot(dbday$mean$DATE, ytot, type="l", col="gray", lwd=0.1, ylog = "T", cex.axis = 2, bty = "n", bg = rgb(1,1,1,1),
       ylim=c(0.99, 1.01)*range(ytot, na.rm=T),
       xlim=p2$xlim, xaxt = "n", xlab = "", ylab = "", xaxs="i")
  # lines(dbday$mean$DATE, dbday$mean$TEMP_ADJUSTED_despiked, lwd=10, col=rgb(0.7,0,0,alpha = 0.7))
  lines(dbday$mean$DATE, ytot, lwd=0.5, col="black")
  abline(v = dateticks1, lty = 3, lwd = 1, col = "darkgray")
  abline(v = c(convperiod.m$startdates,convperiod.m$enddates), lty = 2, lwd = 1, col = "blue")
  mtext(text = expression(paste("Temp. (",Theta, ", ºC)")), side = 2, line = 3.2, cex = 2)
  legend("topleft", cex = 2.0,
         legend = c("daily"),
         col = c("black"),
         lwd = c(1),
         lty = c(1),
         seg.len = c(1),
         bg = "#FFFFFF95", box.col = "#FFFFFF00")
  
  # Common x axis for plot rows 2-8
  axis.POSIXct(side = 1, at = dateticks1, format = "%Y-%m", cex.axis = 2, tck=-0.06, padj = 1) # cex.axis = 1.5
  axis.POSIXct(side = 1, at = dateticks2, labels = "", tck=-0.03)
  
  # -----------------------------------------------------------------------------------------------------
  # Colorbars: with this method they have to be called after all plot windows have been filled
  
  # Colorbar for chlorophyll (rows 2-3)
  v23 <- c( 0.94*vp2[2], 0.95*vp2[2], 0.855*vp2[3], 0.935*vp2[4] )
  par( fig=v23, new=TRUE, mar=c(0,0,0,0) )
  cbar.chl <- seq(p2$zlim[1],p2$zlim[2], length.out = ncol.chl)
  cbar.tick <- c(-3,-2,-1,0,1)
  image(t(cbar.chl), cbar.chl, col = col.chl, xaxt = "n", yaxt = "n", bg = "white") 
  axis(4, cex.axis=1.5, mgp = c(0, 1, 0), at = cbar.tick, labels = 10^cbar.tick, las = 1)
  mtext(side = 4, text = expression(paste("FChl",italic("a")," (mg ",m^-3,")")), las = 0, line = 5, cex = 1.4) #cex = 1.3
  box()
  
  # Colorbar for bbp700 or POC (rows 4-5)
  v45 <- v23 - c(0,0,0.24,0.24) # For some reason, better scaling when recycling vp2 than when using vp4
  par( fig=v45, new=TRUE, mar=c(0,0,0,0) )
  cbar.poc <- seq(p4$zlim[1],p4$zlim[2], length.out = ncol.poc)
  cbar.tick <- c(-4,-3,-2,-1)
  image(t(cbar.poc), cbar.poc, col = col.poc, xaxt = "n", yaxt = "n", bg = "white") 
  axis(4, cex.axis=1.5, mgp = c(0, 1, 0), at = cbar.tick, labels = 10^cbar.tick, las = 1)
  mtext(side = 4, text = expression(paste('b'[bp700],' (',m^-1*')')), line = 6, las = 0, cex = 1.5) #cex = 1.3
  box()
  
  # Colorbar for months (row 1)
  v1 <- v23 + c(0.05,0.05,0.15,0.15)
  par( fig=v1, new=TRUE, mar=c(0,0,0,0) )
  cbar.month <- seq(1, 12)
  cbar.tick <- seq(2, 12, 2)
  image(t(cbar.month), cbar.month, col = rev(col.month), xaxt = "n", yaxt = "n", bg = "white") 
  axis(4, cex.axis=1.5, mgp = c(0, 1, 0), at = cbar.tick, labels = rev(c("Jan","Mar","May","Jul","Sep","Nov")), las = 1, cex = 1.5) #cex = 1.3
  box()
  
  dev.off()
  # stop("END OF PLOTTING")
}

# # --------------------------------------------------------------------------------------------------------------------------
# # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMBINED VERTICAL INTEGRALS AND TRAJ PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # --------------------------------------------------------------------------------------------------------------------------
# # NOTE: Not shown in convection paper. UNCOMMENT SECTION BELOW FOR PLOTTING.

if (plotintsummary) {
  p <- paste0( "~/Desktop/Gali_2026_convectionPOC/output/intsummary_", fwmo , "_", paste(selperiod, collapse = "_"),"_",biospike,".png" )
  png(filename = p, width = pw, height = ph, units = 'cm', pointsize = 8, bg = "white", res = 600, type = "quartz")


  # Multipanel setup
  m1 <- list()
  yyp1 <- length(pyearS)+0 # decide here whether to add (+1) or not (+0) an extra panel for color bar (months) in the top row
  pw1res <- pw%%yyp1
  pw1 <- (pw-pw1res)/yyp1
  for (jy in 1:(yyp1+1)) {
    ifelse( jy==(yyp1+1), m1[[jy]] <- matrix(data = jy, nrow = 4, ncol = pw1res), m1[[jy]] <- matrix(data = jy, nrow = 4, ncol = pw1))
  }
  m1 <- Reduce(cbind, m1)
  m2 <- matrix(data = max(m1)+1, nrow = 3, ncol = pw)
  m3 <- m2+1
  m4 <- m2+2
  m5 <- m2+3
  m6 <- matrix(data = m5+1, nrow = 3, ncol = pw)
  m7 <- m6+1
  m8 <- matrix(data = m7+1, nrow = 3, ncol = pw)
  layout(rbind(m1, m2, m3, m4, m5, m6, m7, m8)) # layout.show(seq(1,10))
  par(oma = c(1,1,0.5,0.5))


  # -----------------------------------------------------------------------------------------------------
  # 1ST ROW: Maps (as many as pyearS) of float trajectory overlaid on bathymetry. Dots colored by month. Contour of maximum March MLD (reanalysis?)
  for (jy in 1:length(pyearS)) {

    yearmonth <- pyearS[jy]+2.5/12
    ymdf <- mdf.orig[ mdf.orig$year==pyearS[jy] , ]

    ifelse( yyp1 > length(pyearS), par(mar = c(1,4,1,0)), par(mar = c(1,4,1,6)) )
    image(gebco$lonvec, gebco$latvec, gebco$zmat, xlab = "", ylab = "", xlim = c(-61,-24), ylim = c(52,65), col = col.bathy,
          xaxt = "n", yaxt = "n", main = "")
    plot(coastline, clon = 0, clat = 0, span = c(length(gebco$lon), length(gebco$lon)),
         col = "black", bg = rgb(1,1,1,0), lwd = 0.5, axes = F, add = T)
    ifelse(jy==1, dl <- T, dl <- F)
    ifelse(jy==1, cl <- 1, cl <- 0.01)
    contour(gebco$lonvec,gebco$latvec,gebco$zmat, levels = c(-1000,-3000), lwd = 0.5, lty = 1, labcex = 1, col = col.bathy[1], drawlabels = dl, add = T)
    zoomm <- f_prep4contour1deg(xin = nav_lon[vcell], yin = nav_lat[vcell],
                                zin = A[which(vtime==yearmonth),,which(vname=="mlotstmax.m")],        # March MLD contours: 2015
                                dlondeg = 0.5, dlatdeg = 0.5)
    contour(zoomm$x, zoomm$y, zoomm$z, levels = c(1000,2500), lwd = c(1,2), lty = 1, cex = 1, col = "black", drawlabels = T, add = T)
    points(ymdf$lon, ymdf$lat, col = col.month[ month(ymdf$date) ], pch = 20, cex = 2.5)
    title(main = pyearS[jy], line = -2, cex.main = 3, col.main = "white", adj = 0.1)

  }
  if(yyp1 > length(pyearS)) {plot.new()} # additional year (added through yyp1 to make room for color bar)
  if (pw1res!=0) plot.new()

  # -----------------------------------------------------------------------------------------------------
  # Vertical integrals

  mdf2plot <- mdf # either mdf.orig or mdf (which may include time binning and regularization)
  mchlfact <- 1
  SIfact <-  mean(mdf2plot$ESDbbp_iz0.100, na.rm=T) / mean(mdf2plot$SIchl_iz0.100, na.rm=T)

  # 2ND ROW: CHLA
  par(mar = c(0,5,0,5))
  ylim <- c(0, 1.1*quantile(mdf[,c("chla_adjusted_iz0.100","chla_adjusted_iz100.500","chla_adjusted_iz500.900")], 0.999, na.rm=T))
  plot(mdf2plot$date, mdf2plot$chla_adjusted_iz0.100, type="l", col="darkgreen", lwd=2, cex.axis = 1.5, bty = "n", bg = rgb(1,1,1,1), # ylog = "T",
       ylim=ylim,
       xlim=selperiod, xaxt = "n", xlab = "", ylab = "", xaxs="i")
  # lines(mdf2plot$date, mdf2plot$chla_adjusted_izprod, col="black", lwd=2)
  lines(mdf2plot$date, mdf2plot$chla_adjusted_iz100.500 * mchlfact, col="limegreen", lwd=2)
  lines(mdf2plot$date, mdf2plot$chla_adjusted_iz500.900 * mchlfact, col=col.chltraj, lwd=3)
  abline(v = dateticks1, lty = 3, lwd = 1, col = "darkgray")
  text(x = dateticks1[1]-as.duration("75 days"), y = 0.8*max(ylim), labels = expression(paste(sum(Chl[despiked]),' (mg ',m^-2,')')), cex = 2, adj = 0)
  vp2 <- f_get_pcoords() # get main plot coordinates to overlay colorbar
  legend("right",
         legend = c(expression(paste(sum(Chl[0-100]))),
                    expression(paste(sum(Chl[100-500]))),
                    expression(paste(sum(Chl[500-900])))),
         col = c("darkgreen","limegreen",col.chltraj), lwd = c(2,2,3), lty = c(1,1,1), bty = "n", bg = "white", cex = 2)

  # 3RD ROW: CHLA spikes
  par(mar = c(0,5,0,5))
  ylim <- c(0, quantile(mdf[,c("chla_adjusted_spike_iz0.100","chla_adjusted_spike_iz100.500","chla_adjusted_spike_iz500.900")], 0.999, na.rm=T))
  plot(mdf2plot$date, mdf2plot$chla_adjusted_spike_iz0.100, type="l", col="darkgreen", lwd=2, cex.axis = 1.5, bty = "n", bg = rgb(1,1,1,1), # ylog = "T",
       ylim=ylim,
       xlim=selperiod, xaxt = "n", xlab = "", yaxt = "n", ylab = "", xaxs="i")
  # lines(mdf2plot$date, mdf2plot$chla_adjusted_spike_izprod, col="black", lwd=2)
  lines(mdf2plot$date, mdf2plot$chla_adjusted_spike_iz100.500 * mchlfact, col="limegreen", lwd=2)
  lines(mdf2plot$date, mdf2plot$chla_adjusted_spike_iz500.900 * mchlfact, col=col.chltraj, lwd=3)
  points(mdf2plot$date, (mdf2plot$SIchl_iz0.100/100) * SIfact, col = "darkgreen", pch = 20, cex = 2)
  points(mdf2plot$date, (mdf2plot$ESDbbp_iz0.100/100), col = "tan4", pch = 1, cex = 1)
  axis(side = 4, cex.axis = 1.5)
  abline(v = dateticks1, lty = 3, lwd = 1, col = "darkgray")
  text(x = dateticks1[1]-as.duration("75 days"), y = 0.8*max(ylim), labels = expression(paste(sum(Chl[spikes]),' (mg ',m^-2,')')), cex = 2, adj = 0)
  legend("right",
         legend = c(expression(paste(ESD[bbp,0-100],' (µm)')),
                    expression(paste(SI[Chl,0-100],' (-)'))),
         col = c("tan4","darkgreen"), pch = c(1,20), bty = "n", bg = "white", cex = 2)

  # 4TH ROW: bbp despiked
  par(mar = c(0,5,0,5))
  ylim <- c(0, 1.1*quantile(mdf[,c("bbp700_iz0.100","bbp700_iz100.500","bbp700_iz500.900")], 0.999, na.rm=T))
  plot(mdf2plot$date, mdf2plot$bbp700_iz0.100, type="l", col="tan4", lwd=2, cex.axis = 1.5, bty = "n", bg = rgb(1,1,1,1), # ylog = "T",
       ylim=ylim,
       xlim=selperiod, xaxt = "n", xlab = "", ylab = "", xaxs="i")
  # lines(mdf2plot$date, mdf2plot$bbp700_izprod, col="black", lwd=2)
  lines(mdf2plot$date, mdf2plot$bbp700_iz100.500, col="tan3", lwd=2)
  lines(mdf2plot$date, mdf2plot$bbp700_iz500.900, col=col.poctraj, lwd=3)
  abline(v = dateticks1, lty = 3, lwd = 1, col = "darkgray")
  text(x = dateticks1[1]-as.duration("75 days"), y = 0.8*max(ylim), labels = expression(paste(sum(bbp[despiked]),' (',m^-3,')')), cex = 2, adj = 0)
  legend("right",
         legend = c(expression(paste(sum(bbp[0-100]))),
                    expression(paste(sum(bbp[100-500]))),
                    expression(paste(sum(bbp[500-900])))),
         col = c("tan4","tan3",col.poctraj), lwd = c(2,2,3), lty = c(1,1,1), bty = "n", bg = "white", cex = 2)

  # 5TH ROW: bbp spikes
  par(mar = c(0,5,0,5))
  ylim <- c(0, quantile(mdf[,c("bbp700_spike_iz0.100","bbp700_spike_iz100.500","bbp700_spike_iz500.900")], 0.999, na.rm=T))
  plot(mdf2plot$date, mdf2plot$bbp700_spike_iz0.100, type="l", col="tan4", lwd=2, cex.axis = 1.5, bty = "n", bg = rgb(1,1,1,1), # ylog = "T",
       ylim=ylim,
       xlim=selperiod, xaxt = "n", xlab = "", yaxt = "n", ylab = "", xaxs="i")
  # lines(mdf2plot$date, mdf2plot$bbp700_spike_izprod, col="black", lwd=2)
  lines(mdf2plot$date, mdf2plot$bbp700_spike_iz100.500, col="tan3", lwd=2)
  lines(mdf2plot$date, mdf2plot$bbp700_spike_iz500.900, col=col.poctraj, lwd=3)
  axis(side = 4, cex.axis = 1.5)
  abline(v = dateticks1, lty = 3, lwd = 1, col = "darkgray")
  text(x = dateticks1[1]-as.duration("75 days"), y = 0.8*max(ylim), labels = expression(paste(sum(bbp[spikes]),' (',m^-3,')')), cex = 2, adj = 0)

  # -----------------------------------------------------------------------------------------------------
  # Trajectory 1000 m plots
  db2plot <- dbbin # either dbday or dbbin (time-binned onto same tgrid as vertical profiles)

  # 6TH ROW: Trajectory data for CHLA
  par(mar = c(0,5,0,5))
  ytot <- dbbin$mean[[paste0("CHLA_ADJUSTED_",biospike,"despiked")]] + dbbin$mean[[paste0("CHLA_ADJUSTED_",biospike,"spike")]]
  plot(dbbin$mean$DATE, ytot, type="l", col="gray", lwd=0.1, ylog = "F", cex.axis = 1.5, bty = "n", bg = rgb(1,1,1,1),
       ylim = c(0, 0.12), # max 0.08 for 6901480, 0.12 for 6901486
       # ylim=c(0.99, 1.1)*quantile(ytot, c(0.001,0.999), na.rm=T),
       xlim=selperiod, xaxt = "n", xlab = "", ylab = "", xaxs="i")
  lines(dbbin$mean$DATE, dbbin$mean[[paste0("CHLA_ADJUSTED_",biospike,"despiked")]], lwd=2, col=col.chltraj)
  lines(dbbin$mean$DATE, ytot, lwd=0.5, col="black")
  abline(v = dateticks1, lty = 3, lwd = 1, col = "darkgray")
  text(x = dateticks1[1]-as.duration("75 days"), y = 0.15, labels = "Chl (1000 m)", cex = 2, adj = 0)

  # 7TH ROW: Trajectory data for BBP700
  par(mar = c(0,5,0,5))
  ytot <- dbbin$mean[[paste0("BBP700_",biospike,"despiked")]] + dbbin$mean[[paste0("BBP700_",biospike,"spike")]]
  plot(dbbin$mean$DATE, ytot, type="l", col="gray", lwd=0.1, ylog = "F", cex.axis = 1.5, bty = "n", bg = rgb(1,1,1,1),
       ylim = c(0.00000,0.00050),
       # ylim=c(0.8, 0.7)*range(ytot, na.rm=T),
       xlim=selperiod, xaxt = "n", yaxt = "n", xlab = "", ylab = "", xaxs="i")
  axis(side = 4, cex.axis = 1.5)
  lines(dbbin$mean$DATE, dbbin$mean[[paste0("BBP700_",biospike,"despiked")]], lwd=2, col=col.poctraj)
  lines(dbbin$mean$DATE, ytot, lwd=0.5, col="black")
  abline(v = dateticks1, lty = 3, lwd = 1, col = "darkgray")
  text(x = dateticks1[1]-as.duration("75 days"), y = 0.0008, labels = "bbp (1000 m)", cex = 2, adj = 0)

  # 8TH ROW: Trajectory data for TEMP
  par(mar = c(4,4,0,4))
  ytot <- dbbin$mean$CTEMP_despiked + dbbin$mean$CTEMP_spike
  plot(dbbin$mean$DATE, ytot, type="l", col="gray", lwd=0.1, ylog = "T", cex.axis = 1.5, bty = "n", bg = rgb(1,1,1,1),
       ylim=c(0.99, 1.01)*range(ytot, na.rm=T),
       xlim=selperiod, xaxt = "n", xlab = "", ylab = "", xaxs="i")
  lines(dbbin$mean$DATE, dbbin$mean$TEMP_ADJUSTED_despiked, lwd=4, col="lightgray")
  lines(dbbin$mean$DATE, ytot, lwd=0.5, col="black")
  abline(v = dateticks1, lty = 3, lwd = 1, col = "darkgray")
  text(x = dateticks1[1]-as.duration("75 days"), y = 4, labels = "Temp (1000 m)", cex = 2, adj = 0)

  # Common x axis for plot rows 2-8
  axis.POSIXct(side = 1, at = dateticks1, format = "%Y-%m", cex.axis = 1.5, tck=-0.06)
  axis.POSIXct(side = 1, at = dateticks2, labels = "", tck=-0.03)

  # Colorbar for months (row 1)
  v23 <- c( 0.95*vp2[2], 0.96*vp2[2], 0.84*vp2[3], 0.92*vp2[4] )
  v1 <- v23 + c(0.04,0.04,0.15,0.15)
  par( fig=v1, new=TRUE, mar=c(0,0,0,0) )
  cbar.month <- seq(1, 12)
  cbar.tick <- seq(1, 11, 2)
  image(t(cbar.month), cbar.month, col = col.month, xaxt = "n", yaxt = "n", bg = "white")
  axis(4, cex.axis=1.5, mgp = c(0, 1, 0), at = cbar.tick, labels = c("Jan","Mar","May","Jul","Sep","Nov"), las = 1)
  box()

  dev.off()
}


# END OF SCRIPT
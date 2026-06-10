# Figure 2 summarizing difference between weak and strong convection periods in terms of diffusive vs. gravitational fluxes
# Based on test_fig_con_fluxes.R, where additional calculations are available

# Libraries
library(RNetCDF)
library(plyr) # for mapvalues function
library(dplyr) # for select function
library(purrr)
library(tidyr) # for pivot_wider
library(reshape) # for rename function
library(data.table) # for rbindlist
library(oce)
library(ocedata)
library(akima) # for interp function
library(RColorBrewer)
library(lubridate)
library(fields)

# -----------------------------------------------------------------------------------------------------
# Settings and definitions for analysis and plotting

# Select outputs
fig_convfluxes <- T           # Figure 4 of POC convection paper
sfig_sections_map <- T        # Figure 16A of POC convection paper
scheme_budgets <- T           # Figure 5 of POC convection paper
if (scheme_budgets) {
  fig_dcv_budget <- T
  fig_convfluxes <- F
  sfig_sections_map <- F
}

# Region to be analyzed
regname <- "IR_LAB" # Default: "IR_LAB". Other options: "LAB_CB" (restricted to Labrador Central basin), "NASPG" (not limited to Labrador + Irminger Seas)

# Time variables
season <- NULL # NULL, "_FMAM", "_JJASO"
yearsREF <- seq(1958,2019)

# Define experiment and analysis period
# Processing functions prepared to operate on list of experiments, but plotting not yet adapted mostly because of different periods covered by each exeriment
expidS <- list(
  # M1hD1h = "a5xl",   # Restoring (extreme), 1h scale for both  (= M1hD1h), 1998-2019: expid <- "a5xl", expdate <- "v20230421"
  M2hD3h = "a683",     # Restoring (intermediate), 3h/2h scale for M/D, 1998-2019:      expid <- "a683", expdate <- "v20231017"
  REF2 = "a67o"        # New reference, no restoring, 1958-2019:                        expid <- "a67o", expdate <- "v20230630"
)
# Simulation periods
yearS <- list(
  M1hD1h = c(1998,2019),
  M2hD3h = c(1998,2019),
  REF2 = c(1958,2019)
)

# Variables
varnames.flux <- unlist(strsplit(x = "expsdetoc.expldetoc.zafsdetoc.zafldetoc.zafpkt.zdfsdetoc.zdfldetoc.zdfpkt",
                                 split = ".", fixed = T))
varnames.conc <- unlist(strsplit(x = "phymisc.phydiat.zmicro.zmeso.sdetoc.ldetoc",
                                 split = ".", fixed = T))
varnames.conc <- sapply(varnames.conc, function(x) paste0("conc_",x))
varnames.phy2D <- unlist(strsplit(x = "omlda.omldamax.mlotst.mlotstmax.tos.sos.taum",
                                  split = ".", fixed = T))
varnames.bgc2D <- unlist(strsplit(x = "dmpphymisc.dmpphydiat.intpp",
                                  split = ".", fixed = T))
# Base paths
mbasepath <- "~/Desktop/Gali_2026_convectionPOC/input_data/fig_4_5_S14-16_nemo-pisces/"

# -----------------------------------------------------------------------------------------------------
# Define several multiyear periods based on mld, fluxes and sdetoc and doc stocks
col.periods <- c("#90a1a3","#2B6CBE")
names(col.periods) <- c("2009_2013","2014_2018")
jalpha <- 40

# -----------------------------------------------------------------------------------------------------
# Function to load all desired variables from a given experiment (expid)
# all fluxes are defined + downwards

f_load_3D_variables <- function(expid, mbasepath, dirname, varnames, search_pattern) {
  
  # Load preprocessed files with specific season selected?
  ncfileS <- list.files(paste0(mbasepath, expid, "/", dirname), pattern = search_pattern, full.names = T)
  # print(ncfileS) # debugging
  
  # Load into list
  datalist <- lapply(varnames, function(vv) {
    fname <- grep(pattern = vv, x = ncfileS, value = T)
    # print(paste0("Loading ",vv)) # debugging
    ncfile <- open.nc(fname)
    return( var.get.nc(ncfile, variable = gsub("conc_","",vv)) )
    close.nc(ncfile)
  })
  names(datalist) <- varnames
  return(datalist)
}
f_load_2D_variables <- function(expid, mbasepath, dirname, varnames, search_pattern) {
  datalist <- lapply(varnames, function(vv) {
    fname <- grep(pattern = vv,
                  x = list.files(paste0(mbasepath, expid, "/", dirname), pattern = search_pattern, full.names = T),
                  value = T)
    # print(fname) # debugging
    ncfile <- open.nc(fname)
    return( var.get.nc(ncfile, variable = vv) ) # all fluxes are defined + downwards. CHECKED
    close.nc(ncfile)
  })
  names(datalist) <- varnames
  return(datalist)
}

# Function to compute horizontal spatial means for list of 4D (xyzt) or 3D (xyt) arrays, including spatial masks (only ocean cells selected; zeros possible and different from NA)
f_xyzt_xymean_mask <- function(Lxyzt, MASK) { # # Crop volcello array to match experiment duration (restoring experiments start on 19980101)
  
  if (length(dim(Lxyzt[[1]]))==4) {
    amargin <- c(3,4)
    Vxyzt <- Lxyzt[["volcello"]] # assuming volcello has been loaded and is in Lxyzt. TO DO: ADD CHECK AND WARNING
  } else if (length(dim(Lxyzt[[1]]))==3) {
    amargin <- 3
    Vxyzt <- array(areacello, dim = dim(Lxyzt[[1]]))  # assuming areacello has been loaded. TO DO: PUT AREACELLO IN Spre LIST
  }
  if (!is.null(MASK)) {
    Vxyzt <- Vxyzt * array(MASK, dim(Vxyzt))
  }
  lapply(Lxyzt, function(A) {
    AxV <- A * Vxyzt # element by element product
    AxVsum <- apply(AxV, MARGIN = amargin, sum, na.rm=T)
    Vsum <- apply(Vxyzt, MARGIN = amargin, sum, na.rm=T)
    return(AxVsum / Vsum)
  })
}
# Function to compute horizontal spatial quantile for 4D (xyzt) or 3D (xyt) arrays, including spatial masks (only ocean cells selected; zeros possible and different from NA)
# NOTE: does not operate on list of arrays, but on individual arrays
f_xyzt_xyquant_mask <- function(A, MASK, QUANTILES) {
  
  if (length(dim(A))==4) {
    amargin <- c(3,4)
  } else if (length(dim(A))==3) {
    amargin <- 3
  }
  if (!is.null(MASK)) {
    Amask <- A * array(MASK, dim(A))
    A[Amask==0] <- NA
  }
  TMP <- list()
  TMP <- lapply(QUANTILES, function(qq) {
    Aquant <- apply(A, MARGIN = amargin, quantile, as.numeric(gsub("q","",qq))/100, na.rm=T)
  })
  names(TMP) <- QUANTILES
  return(TMP)
}
# Function to compute horizontal sums for list of 3D (xyt) or 4D (xyzt) arrays, including spatial masks
# Piece of code replicated in fig_mlotst_omlda_extent_opera_glorys.R
f_xyzt_xysum_mask <- function(Lxyzt, MASK) {
  
  Sxy <- areacello # assuming areacello has been loaded
  
  if (length(dim(Lxyzt[[1]]))==4) {
    amargin <- c(3,4)
  } else if (length(dim(Lxyzt[[1]]))==3) {
    amargin <- 3
  }
  if (!is.null(MASK)) {
    Sxyt <- array( (MASK * Sxy) , dim(Lxyzt[[1]]))
  }
  lapply(Lxyzt, function(A) {
    AxS <- A * Sxyt # element by element product
    AxSsum <- apply(AxS, MARGIN = amargin, sum, na.rm=T)
    return(AxSsum)
  })
}
# Function to get plot coordinates for posterior overlaying
f_get_pcoords <- function() {
  u <- par("usr")
  v <- c(grconvertX(u[1:2], "user", "ndc"), grconvertY(u[3:4], "user", "ndc"))
  return(v) 
}

# -----------------------------------------------------------------------------------------------------
# Load fixed spatial fields
# Load lists of selected cells (inside polygon), Argo profile counts, etc: produced from "map_profile_counts.NASPGcoriolis.R"
NASPGij <- read.csv("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/NASPGij_orca1.csv")
load(paste0("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/inpolygon_NASPGcells_orca1.Rda"))
load(paste0("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/inpolygon_labCBcells_orca1.Rda"))
load(paste0("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/inpolygon_IRLABcells_orca1.Rda"))

# Load ORCA1 horizontal grid, areacello and deptht
load("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/ORCA1_hgrid.Rda")
anc <- open.nc("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/areacello_Ofx_EC-Earth3_opera-control_r1i1p1f1_gn.nc")
areacello <- var.get.nc(anc, "areacello", start = c(NASPGij$istart, NASPGij$jstart), count = c(NASPGij$icounts, NASPGij$jcounts)) # ORCA1 areacello (m2!)
close.nc(anc)
deptht <- read.csv("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/zgrid_L75.csv", header = T)

# Load global bathymetry mask and subset according to NASPGcells (rectangular lon-lat domain)
if (scheme_budgets) {
  # Use this mask for budgets to avoid bottom effects on diffusive-advective transports that confound convection-driven variability
  maskbathy <- var.get.nc(ncfile = open.nc("~/Desktop/OPERA/bathy_mask/bathymask_2500_0_1.nc"),
                          variable = "bathymask",
                          start = c(NASPGij$istart, NASPGij$jstart),
                          count = c(NASPGij$icounts, NASPGij$jcounts))
} else {
  # Use this mask for 1000-m fluxes (Fig. 4) and tables. Fluxes in >2500-m domain show slightly more pronounced variability, but consistent patterns
  maskbathy <- var.get.nc(ncfile = open.nc("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/a2er2000_bathy_mask.nc"),
                          variable = "bathymask",
                          start = c(NASPGij$istart, NASPGij$jstart),
                          count = c(NASPGij$icounts, NASPGij$jcounts))
}
maskocean <- var.get.nc(ncfile = open.nc("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/a2er800_bathy_mask.nc"),
                        variable = "bathymask",
                        start = c(NASPGij$istart, NASPGij$jstart),
                        count = c(NASPGij$icounts, NASPGij$jcounts))
maskbathy[is.nan(maskbathy)] <- 0
maskocean[is.nan(maskocean)] <- 0

# -----------------------------------------------------------------------------------------------------
# Create rectangular masks of labCB and IRLAB domains based on cell indices **and bathymetry**
basemask <- 0*nav_lat
f_makemask <- function(regbounds, maskcells, bathymask) {
  basemask[maskcells$imatch] <- 1
  maskout <- basemask[regbounds$istart:(regbounds$istart+regbounds$icount-1),regbounds$jstart:(regbounds$jstart+regbounds$jcount-1)]
  return( maskout & bathymask==1 )
}
smask <- list(
  LAB_CB = f_makemask(NASPGij, labCBcells, maskbathy),
  IR_LAB = f_makemask(NASPGij, IRLABcells, maskbathy),
  maskocean_IR_LAB = f_makemask(NASPGij, IRLABcells, maskocean),
  NASPG = f_makemask(NASPGij, NASPGcells, maskocean)
)

# -----------------------------------------------------------------------------------------------------
# Load fluxes in cropped NASPG region and 4 selected depths (100, 500, 1000, 2000; in fact gridW levels 25,40,47,55)
# 4D array with dimensions x,y,z,t is convenient: can reuse functions from hovmoller_pisces_budgetPOC_3D-R
Fpre <- list()
for (nn in names(expidS)) {
  expid <- expidS[[nn]]
  Fpre[[nn]] <- f_load_3D_variables(expid, mbasepath, dirname = "year", varnames.flux, search_pattern = "4z.nc")
}
# Load organic C tracer concentrations. Prepend "conc_" to avoid matching tracer flux variables
Cpre <- list()
for (nn in names(expidS)) {
  expid <- expidS[[nn]]
  Cpre[[nn]] <- f_load_3D_variables(expid, mbasepath, dirname = "year", varnames.conc, search_pattern = "4z.nc")
}

# Load surface physical variables. For simplicity, repeat the same REFERENCE expid for each list item that has that REF expid
Spre <- list()
Spre_bgc <- list()
moSpre <- list() # monthly data for annual maximum mlotstmax and omldamax
for (nn in names(expidS)) {
  
  expid <- expidS[[nn]]
  # print(nn) # debugging
  # print(expid)
  
  # If experiment older than a67o, use older reference a5gj for physical fields. Else, use a67o
  findREF <- (sort(c("a67o",expid)))[1]
  expidREF <- ifelse(findREF=="a67o", "a67o", "a5gj")
  Spre[[nn]] <- f_load_2D_variables(expidREF, mbasepath, dirname = "year", varnames.phy2D, search_pattern = "12.nc")
  Spre_bgc[[nn]] <- f_load_2D_variables(expid, mbasepath, dirname = "year", varnames.bgc2D, search_pattern = "12.nc")
  moSpre[[nn]] <- f_load_2D_variables(expidREF, mbasepath, dirname = "month", c("omldamax","mlotstmax"), search_pattern = "12.nc")
  
  # Replace mlotstmax and omldmax in Spre, which are the annual means of monthly maxima, by the absolute annual mean
  for (vv in c("omldamax","mlotstmax")) {
    TMP <- moSpre[[nn]][[vv]]
    dd <- dim(TMP)
    NTMP <- apply(X = array(TMP, dim = c(dd[1:2], 12, dd[3]/12)), MARGIN = c(1,2,4), function(A) {
      Amax <- max(A, na.rm = T)
      Amax[abs(Amax)==Inf] <- NA
      return(Amax)
    })
    Spre[[nn]][[vv]] <- NTMP
  }
  
  # Remove years in ref experiment for physical variables (eg a5gj) that do not match
  # the sim years in the expid, eg a satellite nudged exp (a5xl)
  if (expid!=expidREF) {
    yexp <- seq(yearS[[nn]][1],yearS[[nn]][2])
    yref <- 1958:2019
    yind <- which(yref %in% yexp)
    Spre[[nn]] <- lapply(Spre[[nn]], function(x) return(x[,,yind]))
  }
  
  # Concatenate lists
  Spre[[nn]] <- c(Spre[[nn]],Spre_bgc[[nn]])
  
}
rm(moSpre)

# -----------------------------------------------------------------------------------------------------
# Add satellite OCCCIv5 data to Spre. If needed, follow same procedure to add phymisc and phydiat from satellite
# Primary production (Sathyendranath et al. 2020, OC-CCI-v5)
load("~/Desktop/Gali_2026_convectionPOC/input_data/fig_4_5_S14-16_nemo-pisces/df_subsampled_esa-pp_opera_NASPG_pp.msk.Rda")
ocyears <- seq(1998,2019)
pp.OCCCIv5 <- array(datadf5$pp.OCCCIv5, dim = c( length(ocyears), 12, dim(areacello) ))
pp.OCCCIv5 <- aperm(pp.OCCCIv5, perm = c(3,4,1,2))
pp.OCCCIv5 <- apply(pp.OCCCIv5, MARGIN = c(1,2,3), weighted.mean, w = c(31,28,31,30,31,30,31,31,30,31,30,31))
for (nn in names(expidS)) {
  expid <- expidS[[nn]]
  dd <- dim(Spre[[nn]][["intpp"]])
  Spre[[nn]][["ppsat"]] <- array(NA, dd)
  Spre[[nn]][["ppsat"]][,,(dd[3]-21):(dd[3])] <- pp.OCCCIv5 / 1000 / 86400
}

# -----------------------------------------------------------------------------------------------------
# Load volcello from the reference experiment of each expid, adjust (crop) time dimension if needed, and add to Fpre and Cpre lists
for (nn in names(expidS)) {
  expid <- expidS[[nn]]
  # If experiment older than a67o, use older reference a5gj for physical fields. Else, use a67o
  findREF <- (sort(c("a67o",expid)))[1]
  expidREF <- ifelse(findREF=="a67o", "a67o", "a5gj")
  TMPvolcello <- f_load_3D_variables(expidREF, mbasepath, dirname = "year", varnames = "volcello", search_pattern = "")[[1]]
  # Add volcello to the corresponding experiment. If needed, crop it to match experiment duration (restoring experiments start on 19980101)
  # Add volcello to vertical fluxes: it is used in posterior array manipulations
  if (dim(Fpre[[nn]][[1]])[4] < length(yearsREF)) {
    Fpre[[nn]][["volcello"]] <- TMPvolcello[,,,yearsREF%in%seq(yearS[[nn]][1],yearS[[nn]][2])]
  } else {
    Fpre[[nn]][["volcello"]] <- TMPvolcello
  }
  # Add volcello to concentrations
  if (dim(Cpre[[nn]][[1]])[4] < length(yearsREF)) {
    Cpre[[nn]][["volcello"]] <- TMPvolcello[,,,yearsREF%in%seq(yearS[[nn]][1],yearS[[nn]][2])]
  } else {
    Cpre[[nn]][["volcello"]] <- TMPvolcello
  }
  rm(TMPvolcello)
}

# -----------------------------------------------------------------------------------------------------
# Compute spatial means over the chosen domain: IR_LAB (default), NASPG, or LAB_CB
F_zt <- lapply(Fpre, f_xyzt_xymean_mask, smask[[regname]])
S_zt <- lapply(Spre, f_xyzt_xymean_mask, smask[[regname]])

# Compute spatial integrals of fluxes over entire domain
F_ztsum <- lapply(Fpre, f_xyzt_xysum_mask, smask[[regname]])

# Subset data, aggregate variables and convert units
for (nn in names(expidS)) {
  lnames <- names(F_ztsum[[nn]])
  lnames <- grep("doc", lnames, value = T, invert = T) # exclude (DOC, volcello)
  lnames <- grep("volcello", lnames, value = T, invert = T) # exclude (DOC, volcello)
  F_ztsum[[nn]][["exppoc"]] <- Reduce('+', F_ztsum[[nn]][grep("exp", lnames, value = T)])
  F_ztsum[[nn]][["zdfpoc"]] <- Reduce('+', F_ztsum[[nn]][grep("zdf", lnames, value = T)])
  F_ztsum[[nn]][["zdfdetoc"]] <- Reduce('+', F_ztsum[[nn]][c("zdfsdetoc","zdfldetoc")])
  F_ztsum[[nn]][["zafpoc"]] <- Reduce('+', F_ztsum[[nn]][grep("zaf", lnames, value = T)])
  
  # Convert units to Tg C yr-1
  F_ztsum[[nn]] <- lapply(F_ztsum[[nn]], function(x) x * 12.011*86400*365/1e12 )
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate difference between fluxes at 4 depths (100-500-100-2000 m) and 500-2000 m for Table S5
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for (nn in names(expidS)) {
  F_ztsum[[nn]] <- lapply(F_ztsum[[nn]], function(x) return(rbind(x, x[1,]-x[2,], x[2,]-x[3,], x[3,]-x[4,], x[2,]-x[4,])))
}


# -----------------------------------------------------------------------------------------------------
# Compute total convected volume, and convected volume between zconvmin-zconvmax (usually 500-2000 m)
zconvmin <- 500
zconvmax <- 2000
S_ztsum <- lapply(Spre, function(x) {
  x1 <- x[c("omldamax","mlotstmax")]
  x2 <- lapply(x1, function(xx) {
    xx[xx>(zconvmax)] <- zconvmax
    xx <- xx-zconvmin
    xx[xx<0 & !is.na(xx)] <- 0
    return(xx)
  })
  names(x2) <- paste0(names(x1),"_500_2000")
  return(c(x1,x2))
})
S_ztsum <- lapply(S_ztsum, f_xyzt_xysum_mask, smask[[regname]])


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define all periods to be used for statistics
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allperiods <- list(sat = seq(1998,2019),
                   all = seq(1958,2019),
                   weak1 = seq(1967,1971),
                   strong1 = seq(1972,1976),
                   weak2 = seq(1977,1981),
                   strong2 = seq(1989,1993),
                   weak3 = seq(2009,2013),
                   strong3 = seq(2014,2018))
allperiods$allweak <- as.vector(do.call(c, allperiods[grep("weak",names(allperiods))]))
allperiods$allstrong <- as.vector(do.call(c, allperiods[grep("strong",names(allperiods))]))

# -----------------------------------------------------------------------------------------------------
# Compute percentiles of annual mlotstmax within smask (usually IR_LAB)
f_smask <- function(A, MASK) {
  dd <- dim(A)
  OUT <- A*array(MASK, dim(A))
  OUT[OUT==0 & !is.na(OUT)] <- NA
  return(OUT)
}
m.weak <- f_smask(Spre$REF2$mlotstmax[,,yearsREF%in%allperiods$weak3], smask$IR_LAB)
m.strong <- f_smask(Spre$REF2$mlotstmax[,,yearsREF%in%allperiods$strong3], smask$IR_LAB)
qq <- c(0.25,0.5,0.75) # c(0.25,0.5,0.75), 1SD is c(0.16,0.5,0.84), 2SD is c(0.025,0.5,0.975)
# quantile(m.weak, qq, na.rm=T)
# quantile(m.strong, qq, na.rm=T)
# mean(m.weak,  na.rm=T)
# mean(m.strong,  na.rm=T)
fweak <- ecdf(x = m.weak)
fstrong <- ecdf(x = m.strong)
# print(fweak(c(500,1000,2000)))
# print(fstrong(c(500,1000,2000)))
# summary(S_ztsum$REF2$mlotstmax_500_2000[yearsREF%in%allperiods$weak3])
# summary(S_ztsum$REF2$mlotstmax_500_2000[yearsREF%in%allperiods$strong3])
# hist(m.weak, 20)
# hist(m.strong, 20)

# -----------------------------------------------------------------------------------------------------
# Convection diagnostics in IR_LAB region in OPERA vs GLORYS
# Load pre-computed MLD metrics: annual means of mlotstmax and omldmax, and area where they exceed a given threshold within the IR_LAB region
load("~/Desktop/Gali_2026_convectionPOC/input_data/fig_4_5_S14-16_nemo-pisces/conv_extent_opera_1958_2019.Rda")         # OPERA only
load("~/Desktop/Gali_2026_convectionPOC/input_data/fig_4_5_S14-16_nemo-pisces/conv_extent_opera_glorys_1998_2019.Rda")  # OPERA and GLORYS overlap period (27 years)

# -----------------------------------------------------------------------------------------------------
# Compute percentage area with annual mlotstmax exceeding 500, 1000, 2000 (data for manuscript text)
# Areas in million m2 (==km2) for consistency with df*yarea data frames
oceanarea.NASPG <- sum(areacello[smask$NASPG==1], na.rm=T)*1e-6
oceanarea.IR_LAB <- sum(areacello[smask$maskocean_IR_LAB==1], na.rm=T)*1e-6
budgetarea.IR_LAB <- sum(areacello[smask$IR_LAB==1], na.rm=T)*1e-6
df.yrarea <- dplyr::select(df.yrarea, -c(omlda_max_500,omlda_max_1000,omlda_max_2000))
avars <- grep("_", names(df.yrarea), value = T)
toappend <- df.yrarea[,..avars]/oceanarea.IR_LAB
names(toappend) <- paste0("frac_",avars)
df.yrarea <- cbind(df.yrarea,toappend)

# -----------------------------------------------------------------------------------------------------
# NAO time series, based on script nao_crudata_uea.R
dfnao <- read.table(file = "~/Desktop/Gali_2026_convectionPOC/input_data/fig_4_5_S14-16_nemo-pisces/nao_crudata_uea.txt")
dfnao[dfnao==-99.99] <- NA
# Calculate DJF and DJFM means
dfnao$djf <- sapply(1:dim(dfnao)[1], function(j) {
  ifelse(j==1,
         return(NA),
         return( mean( unlist( c(dfnao[j-1,12],dfnao[j,1:2]) ) , na.rm=F ) )
  )
})
dfnao$djfm <- sapply(1:dim(dfnao)[1], function(j) {
  ifelse(j==1,
         return(NA),
         return( mean( unlist( c(dfnao[j-1,12],dfnao[j,1:3]) ) , na.rm=F ) )
  )
})

# Calculate 5-year running means
dfnao$djf.runmean <- caTools::runmean(x = dfnao$djf, k = 5, alg = "C", endrule = "NA", align = "right")
dfnao$djfm.runmean <- caTools::runmean(x = dfnao$djfm, k = 5, alg = "C", endrule = "NA", align = "right")
dca.runmean <- caTools::runmean(x = df.opera.yrarea$mlotst_max_1000, k = 5, alg = "C", endrule = "NA", align = "right")
dcv.runmean <- caTools::runmean(x = S_ztsum$REF2$mlotstmax_500_2000, k = 5, alg = "C", endrule = "NA", align = "right")
# Crop to "selperiod"
dfnao <- dfnao[as.numeric(row.names(dfnao)) %in% seq(1958,2019),]

# # NAO-DCA/DCV correlations: print p-value, r and 95% CI only
# print( cor.test(x = dfnao$djfm, y = df.opera.yrarea$mlotst_max_1000)[c(3, 4, 9)] )
# print( cor.test(x = dfnao$djfm, y = S_ztsum$REF2$mlotstmax_500_2000)[c(3, 4, 9)] )
# print( cor.test(x = dfnao$djfm.runmean, y = dca.runmean)[c(3, 4, 9)] )
# print( cor.test(x = dfnao$djfm.runmean, y = dcv.runmean)[c(3, 4, 9)] )

# # Test significant differences by depth and period
# tz <- 2
# tx <- "REF2" # REF2 M1hD1h M2hD3h
# tv <- "exppoc"
# t.test(x = F_ztsum[[tx]][[tv]][tz,seq(yearS[[tx]][1],yearS[[tx]][2])%in%allperiods$weak3],
#        y = F_ztsum[[tx]][[tv]][tz,seq(yearS[[tx]][1],yearS[[tx]][2])%in%allperiods$strong3])

# -----------------------------------------------------------------------------------------------------
# Load monthly fluxes averaged in multiyear periods in cropped NASPG region, all depths
Fperiods <- list()
for (tt in c("ymonmean","yearmean")) {
  for (yy in c("2009_2013","2014_2018")) {
    pname <- paste0(tt,"_",yy)
    Fperiods[[pname]] <- list()
    for (nn in names(expidS)) {
      expid <- expidS[[nn]]
      Fperiods[[pname]][[nn]] <- f_load_3D_variables(expid, mbasepath, dirname = paste0("month_",yy), varnames.flux, search_pattern = tt)
      # Sum grav, diff, adv fluxes of different tracers
      lnames <- names(Fperiods[[pname]][[nn]])
      lnames <- grep("doc", lnames, value = T, invert = T) # exclude DOC
      Fperiods[[pname]][[nn]][["exppoc"]] <- Reduce('+', Fperiods[[pname]][[nn]][grep("exp", lnames, value = T)])
      Fperiods[[pname]][[nn]][["zdfpoc"]] <- Reduce('+', Fperiods[[pname]][[nn]][grep("zdf", lnames, value = T)])
      Fperiods[[pname]][[nn]][["zdfdetoc"]] <- Reduce('+', Fperiods[[pname]][[nn]][c("zdfsdetoc","zdfldetoc")])
      Fperiods[[pname]][[nn]][["zafpoc"]] <- Reduce('+', Fperiods[[pname]][[nn]][grep("zaf", lnames, value = T)])
    }
  }
}

# -----------------------------------------------------------------------------------------------------
# Compute spatial quantiles for each month or year over the chosen domain (NASPG, IR_LAB, LAB_CB)
#     Fquant list contains "ymonmean_2009_2013", "ymonmean_2014_2018", "yearmean_2009_2013", "yearmean_2014_2018"
#     Each of them has nested lists with expids, varnames, quantiles
Fquant <- Fperiods
for (pname in names(Fquant)) {
  for (nn in names(expidS)) {
    Fquant[[pname]][[nn]] <- lapply(Fquant[[pname]][[nn]], function(X) {f_xyzt_xyquant_mask(X, MASK = smask[[regname]], QUANTILES = c("q25","q50","q75")) })
  }
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# -------------------------------------------------------- FIG. 4 -------------------------------------------------------

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (fig_convfluxes) {
  
  # p <- "~/Desktop/Gali_2026_convectionPOC/output/Fig_4_conv_fluxes.png"
  # png(filename = p, width = 15, height = 11, units = 'cm', pointsize = 8, bg = "white", res = 600, type = "quartz")
  
  p <- "~/Desktop/Gali_2026_convectionPOC/output/Fig_4_conv_fluxes.tiff"
  tiff(filename = p, width = 15, height = 11, units = 'cm', pointsize = 8, bg = "white", res = 300, compression = "lzw")
  
  
  # Multipanel setup
  # TOP ROW
  mm <- matrix(data = -1, nrow = 3, ncol = 11)
  m0 <- matrix(data = 0, nrow = 3, ncol = 11)
  m1 <- matrix(data = 1, nrow = 8, ncol = 11)
  m2 <- matrix(data = 2, nrow = 6, ncol = 4)
  m3 <- matrix(data = 3, nrow = 8, ncol = 4)
  mtop <- cbind(rbind(mm,m0,m1), rbind(m2,m3))
  # BOTTOM ROW
  m4 <- matrix(data = 4, nrow = 9, ncol = 4)
  m5 <- matrix(data = 5, nrow = 9, ncol = 3)
  m6 <- matrix(data = 6, nrow = 9, ncol = 4)
  m7 <- matrix(data = 7, nrow = 9, ncol = 4)
  mbot <- cbind(m4, m5, m6, m7)
  
  layout(rbind(mtop+2, mbot+2))
  par(oma = c(1,1,0,0))
  
  # Conversion factor
  convfact <- 12 * 86400 * 1000 # from mol m-2 s-1 to mgC m-2 d-1
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Fig. 4A. Annual time series of spatial means for selected experiments, variables, depths and domains
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # x axis (date) main and secondary ticks
  decdate <- seq(1958,2019,1) + 0.5/12
  date <- floor_date(date_decimal(decdate, tz = "GMT"), "year") + 0.5*as.duration("year")
  selperiod <- force_tz(as.POSIXct(c(paste0("1998-01-01"),paste0("2020-06-01"))), tzone = "GMT")
  jselperiod <- which(date >= selperiod[1] & date <= selperiod[2])
  
  Date <- as.Date(date)
  dateticks1 <- seq(floor_date(min(date), unit = "10 years"), floor_date(max(date+as.duration("1 year")), unit = "10 years"), by = "10 years")
  dateticks2 <- seq(floor_date(min(date), unit = "years"), floor_date(max(date), unit = "years"), by = "years")
  
  # TOP: NAO time series
  par(mar = c(0,5,2,1))
  plot(date, dfnao$djfm, type="l", col = "black", ylim = c(-3.5,3.5), lwd = 1, ylab = "NAO", las = 1,
       xlim = force_tz(as.POSIXct(c(paste0("1958-01-01"),paste0("2020-12-31"))), tzone = "GMT"),
       cex.axis = 1.2, cex.lab = 1.2, xaxs = "i", yaxs = "i", bty = "n", xaxt = "n", yaxt = "n")
  lines(date, dfnao$djfm.runmean, col = "gray", lwd = 2)
  lines(date, dfnao$djfm, col = "black", lwd = 1)
  ineg <- dfnao$djfm<=0
  abline(v = date[54], lwd = 28.5, col = adjustcolor( col.periods[1], alpha.f = 0.3))
  abline(v = date[59], lwd = 28.5, col = adjustcolor( col.periods[2], alpha.f = 0.3))
  abline(v = date[c(12,22)], lwd = 28.5, col = adjustcolor( col.periods[1], alpha.f = 0.10))
  abline(v = date[c(17,34)], lwd = 28.5, col = adjustcolor( col.periods[2], alpha.f = 0.10))
  points(date[ineg], dfnao$djfm[ineg], pch = 16, col = "red", cex = 1)
  points(date[!ineg], dfnao$djfm[!ineg], pch = 16, col = "blue", cex = 1)
  lines(force_tz(as.POSIXct(c(paste0("1958-01-01"),paste0("2019-12-31"))), tzone = "GMT"), y = c(0,0), col="darkgray", lwd=0.2)
  axis.POSIXct(side = 2, at = seq(-4,4,2), labels = c("",seq(-2,4,2)), cex.axis = 1.2, lwd = 0.5, tck=-0.08, las = 1)
  text(labels = c("Wk1","Str1","Wk2","Str2","Wk3","Str3"), font = c(rep(1,4), rep(2,2), cex = 1.1),
       x = date[c(12,17,22,34,54,59)], y = rep(4.5, 6), col = rep(col.periods,3), cex = 1.3, xpd=TRUE)
  mtext(side = 3, adj = 0, "A", cex = 1.4, font = 2)
  
  # MIDDLE: mlotstmax time series  ----
  df.opera.yrarea$mlotst_max_1000.runmean <- caTools::runmean(x = df.opera.yrarea$mlotst_max_1000, k = 5, alg = "C", endrule = "NA", align = "right")
  par(mar = c(0,5,0,1))
  plot(date, df.opera.yrarea$mlotst_max_1000/1e6, type="l", col = "black", ylim = c(-0.2,1.2), lwd = 1, las = 1,
       ylab = expression(paste("DCA (M",km^2,")")),
       xlim = force_tz(as.POSIXct(c(paste0("1958-01-01"),paste0("2020-12-31"))), tzone = "GMT"),
       cex.axis = 1.2, cex.lab = 1.2, xaxs = "i", yaxs = "i", bty = "n", xaxt = "n", yaxt = "n")
  lines(date, df.opera.yrarea$mlotst_max_1000.runmean/1e6, col = "gray", lwd = 2)
  lines(date, df.opera.yrarea$mlotst_max_1000/1e6, col = "black", lwd = 1)
  abline(v = date[54], lwd = 28.5, col = adjustcolor( col.periods[1], alpha.f = 0.3))
  abline(v = date[59], lwd = 28.5, col = adjustcolor( col.periods[2], alpha.f = 0.3))
  abline(v = date[c(12,22)], lwd = 28.5, col = adjustcolor( col.periods[1], alpha.f = 0.10))
  abline(v = date[c(17,34)], lwd = 28.5, col = adjustcolor( col.periods[2], alpha.f = 0.10))
  lines(force_tz(as.POSIXct(c(paste0("1958-01-01"),paste0("2019-12-31"))), tzone = "GMT"), y = c(0,0), col="darkgray", lwd=0.2)
  axis.POSIXct(side = 2, at = seq(0,1,0.5), labels = seq(0,1,0.5), cex.axis = 1.2, lwd = 0.5, tck=-0.06, las = 1)
  
  
  # BOTTOM: Annual Grav and Diff flux time series ----
  # Prepare data for annual time series plotting: subset variables and depths, create aggregated flux variables
  pvarnames <- c("expsdetoc","expldetoc","zdfsdetoc","zdfldetoc","zdfpkt")
  pvarnames <- c("expsdetoc","expldetoc","zdfsdetoc","zdfldetoc","zdfpkt","zafsdetoc","zafldetoc","zafpkt")
  jz <- 3 # 1000 m in 4z files (1,2,3,4) correspond to (100,500,1000,2000)
  
  yplot <- list()
  for (nn in names(expidS)) {
    TMP <- as.data.frame(lapply(F_zt[[nn]][pvarnames], function(x) x[jz,])) * 12 * 86400 * 1000 # mgC m-2 d-1
    TMP$zdfpoc <- Reduce('+', TMP[,grep("zdf", names(TMP), value = T)])
    TMP$zafpoc <- Reduce('+', TMP[,grep("zaf", names(TMP), value = T)])
    TMP$exppoc <- Reduce('+', TMP[,grep("exp", names(TMP), value = T)])
    yplot[[nn]] <- TMP; rm(TMP)
  }
  col.exppoc <- "khaki3"
  col.zafpoc <- "khaki4"
  col.zdfpoc <- "seagreen"
  
  par(mar = c(5,5,0,1))
  plot(date, yplot$REF2$exppoc, type = "l", lwd = 0.1, col = col.exppoc, lty = 1, las = 1, ljoin = 2,
       xlim = force_tz(as.POSIXct(c(paste0("1958-01-01"),paste0("2020-12-31"))), tzone = "GMT"),
       bty = "n", xaxt = "n", yaxt = "n", cex.axis = 1.2, cex.lab = 1.2, xaxs="i", yaxs = "i", xlab = "",
       ylim = c(-2,10.5), ylab = expression(paste("Flux (mg C ",m^-2," ",d^-1,")")) )
  abline(h = 0, col = "darkgray", lty = 1, lwd = 0.5)
  abline(v = date[54], lwd = 28.5, col = adjustcolor( col.periods[1], alpha.f = 0.3))
  abline(v = date[59], lwd = 28.5, col = adjustcolor( col.periods[2], alpha.f = 0.3))
  abline(v = date[c(12,22)], lwd = 28.5, col = adjustcolor( col.periods[1], alpha.f = 0.10))
  abline(v = date[c(17,34)], lwd = 28.5, col = adjustcolor( col.periods[2], alpha.f = 0.10))
  lines(date, yplot$REF2$exppoc, lwd = 1.5, col = col.exppoc, lty = 1, ljoin = 2)
  lines(date, yplot$REF2$zafpoc, lwd = 1.5, col = col.zafpoc, lty = 1, ljoin = 2)
  lines(date, yplot$REF2$zdfpoc, lwd = 2, col = col.zdfpoc, lty = 1, ljoin = 2)
  lines(date[jselperiod], yplot$M2hD3h$exppoc, lwd = 1.5, col = col.exppoc, lty = 3, ljoin = 2)
  lines(date[jselperiod], yplot$M2hD3h$zafpoc, lwd = 1.5, col = col.zafpoc, lty = 3, ljoin = 2)
  lines(date[jselperiod], yplot$M2hD3h$zdfpoc, lwd = 2, col = col.zdfpoc, lty = 3, ljoin = 2)
  axis.POSIXct(side = 1, mgp = c(0, 0.6, 0), at = dateticks1, format = "%Y", cex.axis = 1.2, lwd = 0.5, tck=-0.03)
  axis.POSIXct(side = 1, at = dateticks2, labels = "", lwd = 0.5, tck=-0.01)
  axis(side = 2, at = seq(-2,10,2), labels = seq(-2,10,2), cex.axis = 1.2, lwd = 0.5, tck=-0.03, las = 1)
  mtext(side = 1, line = 2, "Year", cex = 0.8)
  legend(x = date[1], y = 7.5, cex = 1,  bg = "#FFFFFF95", box.col = "#1a1a1a10",
         legend = rep("",3),
         col = c(col.exppoc, col.zdfpoc, col.zafpoc),
         lty = rep(1,3),  lwd = rep(1.5,3), seg.len = rep(2,3))
  legend(x = date[6], y = 7.5, cex = 1,  bg = "#FFFFFF95", box.col = "#1a1a1a10",
         legend = c("Grav",
                    expression(paste("Mix"[Z])),
                    expression(paste("Adv"[Z]))),
         col = c(col.exppoc, col.zdfpoc, col.zafpoc),
         lty = rep(3,3),  lwd = rep(1.5,3), seg.len = rep(2,3))
  text(x = date[c(3,8)], y = 8, labels = c("REF","SAT"), cex = 1)
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Fig. 4B. Maps of POC fluxes
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Load and preprocess grids and bathymetry
  # Functions inherited from map_profile_counts.NASPG.R
  f_crop4map <- function(xmat, ymat, zmat, zmin, zmax, xindvec, yindvec) {
    xmat <- xmat[ xindvec , yindvec ]
    ymat <- ymat[ xindvec , yindvec ]
    zmat <- zmat[ xindvec , yindvec ]
    zmat[zmat<=zmin] <- NA
    zmat[zmat>zmax] <- zmax # very important when adjusting z scale manually! otherwise values >zmax displayed as NA
    return(list(xmat=xmat, ymat=ymat, zmat=zmat))
  }
  f_prep4contour <- function(xin, yin, zin, dlondeg, dlatdeg) {
    iinclude <- !(zin==0 | is.na(zin) | abs(zin)==Inf)
    fld <- akima::interp(x = xin[iinclude], y = yin[iinclude], z = zin[iinclude],
                         xo = seq( min(floor(xin)), max(ceiling(xin)), dlondeg),
                         yo = seq( min(floor(yin)), max(ceiling(yin)), dlatdeg),
                         linear = T, extrap = F)
  }
  load("~/Desktop/OPERA/grids/ORCA1_hgrid.Rda")
  rlon <- c(-66,-19)
  rlat <- c(49,67)
  
  # Bathymetries: 0.5 degree GEBCO2008, 1/12 degree GLORYS12v1
  load("~/Desktop/Gali_2026_convectionPOC/input_data/gebco_08_05degr.Rda")
  bathy.5 <- list(lonvec = gebco$lonvec, latvec = gebco$latvec, zmat = gebco$zmat)
  bathy.5$zmat[bathy.5$zmat>0] <- NA
  bathy.5$zmat[bathy.5$zmat<(-6000)] <- NA
  globathy <- RNetCDF::read.nc( open.nc("~/Desktop/Gali_2026_convectionPOC/input_data/bathymetry_GLORYS12v1_NASPG.nc") ) # GLORYS12v1
  bathy.08 <- list(lonvec = globathy$longitude, latvec = globathy$latitude, zmat = -globathy$deptho)
  bathy.08$zmat[bathy.08$zmat>0] <- NA
  bathy.08$zmat[bathy.08$zmat<(-6000)] <- NA
  
  # Coastline
  data("coastlineWorldMedium"); coastline <- coastlineWorldMedium
  
  # mlotst_max from OPERA: load and reshape array to calculate annual maximum from monthly maxima
  load("~/Desktop/Gali_2026_convectionPOC/input_data/subsampled_opera_a67o_v20230630_omlda.omldamax.mlotst.mlotstmax.tos.sos.taum.Rda")
  dfa <- as.data.frame(edata)
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
  Amlotstmax.2009_2013 <- apply( Amax[,which(vname=="mlotstmax.m"),unique(dfa$year)%in%allperiods$weak3], 1, max, na.rm=T)
  Amlotstmax.2014_2018 <- apply( Amax[,which(vname=="mlotstmax.m"),unique(dfa$year)%in%allperiods$strong3], 1, max, na.rm=T)
  
  # Color palette for maps
  n.diff <- 7
  col.diff <- c("ivory", brewer.pal(9, "YlGn")[3:9])
  
  # Nicer IRLAB polygon for plotting
  plotIRLAB <- IRLABpoly[-2,]
  plotIRLAB[4,] <- c(-55, 50)
  plotIRLAB[plotIRLAB==65] <- 66
  
  # -----------------------------------------------------------------------------------------------------
  # Fig. 4B-top, 2009-2013
  tmin <- as.Date("2009-01-01")
  tmax <- as.Date("2013-12-31")
  jt <- which(Date>tmin & Date<tmax)
  mapdata <- list(
    ZfluxGravSink = apply(Reduce('+', Fpre[["REF2"]][c("expsdetoc","expldetoc")])[,,jz,jt], c(1,2), mean, na.rm=T),
    ZfluxDiffPOC = apply(Reduce('+', Fpre[["REF2"]][c("zdfpkt","zdfsdetoc","zdfldetoc")])[,,jz,jt], c(1,2), mean, na.rm=T),
    ZfluxAdvPOC = apply(Reduce('+', Fpre[["REF2"]][c("zafpkt","zafsdetoc","zafldetoc")])[,,jz,jt], c(1,2), mean, na.rm=T)
  )
  mapdata <- lapply(mapdata, function(x) {x[smask$NASPG==0] <- NA; x[x==0] <- NA; return(x)}) # NOTE some areas with small upwards (negative) Diff. Mask them with x[x<=0] <- NA
  mapdata$DiffContribution <- mapdata$ZfluxDiffPOC / (mapdata$ZfluxGravSink + mapdata$ZfluxDiffPOC)
  tomap <- f_prep4contour(xin = nav_lon[vcell], yin = nav_lat[vcell],
                          zin = mapdata$DiffContribution,
                          dlondeg = 1/12, dlatdeg = 1/12)
  mldmap <- f_prep4contour(xin = nav_lon[vcell], yin = nav_lat[vcell],
                           zin = Amlotstmax.2009_2013,
                           dlondeg = 0.5, dlatdeg = 0.5)
  # Remove interpolated data shallower than 1000 m
  ixb <- tomap$x>=min(bathy.08$lonvec) & tomap$x<=max(bathy.08$lonvec)
  iyb <- tomap$y>=min(bathy.08$latvec) & tomap$y<=max(bathy.08$latvec)
  tomap$x <- tomap$x[ixb]; tomap$y <- tomap$y[iyb]; tomap$z <- tomap$z[ixb,iyb]
  tomap$z[bathy.08$zmat>-deptht[46,1]] <- NA
  
  par(mar = c(1,1,2,1))
  image(tomap$x, tomap$y, tomap$z, xlim = rlon, ylim = rlat, zlim = c(-0.1,n.diff/10),
        xaxt = "n", yaxt = "n", xlab="", ylab="", col = col.diff)
  plot(coastline, clon = 0, clat = 0, col = "black", bg = rgb(1,1,1,0), lwd = 0.5, axes = F, add = T)
  contour(bathy.5$lonvec,bathy.5$latvec,bathy.5$zmat, levels = c(-1000,-2000), lwd = 0.5, lty = 1, labcex = 1, col = "tan4", drawlabels = F, add = T)
  contour(mldmap$x, mldmap$y, mldmap$z, levels = c(1000), lwd = c(2), lty = 1, cex = 1, col = col.periods[1], drawlabels = F, add = T)
  lines(plotIRLAB, lwd = 1, lty = 3, col = "goldenrod")
  text(x = -23, y = 51, labels = "Wk3", col = col.periods[1], cex = 1.4, font = 2)
  mtext(side = 3, adj = 0, "B", line = -0.5, cex = 1.4, font = 2)
  
  # -----------------------------------------------------------------------------------------------------
  # Fig. 4B-bottom. SAME AS ABOVE FOR 2014-2018 PERIOD
  tmin <- as.Date("2014-01-01")
  tmax <- as.Date("2018-12-31")
  jt <- which(Date>tmin & Date<tmax)
  mapdata <- list(
    ZfluxGravSink = apply(Reduce('+', Fpre[["REF2"]][c("expsdetoc","expldetoc")])[,,jz,jt], c(1,2), mean, na.rm=T),
    ZfluxDiffPOC = apply(Reduce('+', Fpre[["REF2"]][c("zdfpkt","zdfsdetoc","zdfldetoc")])[,,jz,jt], c(1,2), mean, na.rm=T),
    ZfluxAdvPOC = apply(Reduce('+', Fpre[["REF2"]][c("zafpkt","zafsdetoc","zafldetoc")])[,,jz,jt], c(1,2), mean, na.rm=T)
  )
  mapdata <- lapply(mapdata, function(x) {x[smask$NASPG==0] <- NA; x[x==0] <- NA; return(x)}) # NOTE some areas with small upwards (negative) Diff. Mask them with x[x<=0] <- NA
  mapdata$DiffContribution <- mapdata$ZfluxDiffPOC / (mapdata$ZfluxGravSink + mapdata$ZfluxDiffPOC)
  tomap <- f_prep4contour(xin = nav_lon[vcell], yin = nav_lat[vcell],
                          zin = mapdata$DiffContribution,
                          dlondeg = 1/12, dlatdeg = 1/12)
  mldmap <- f_prep4contour(xin = nav_lon[vcell], yin = nav_lat[vcell],
                           zin = Amlotstmax.2014_2018,
                           dlondeg = 0.5, dlatdeg = 0.5)
  # Remove interpolated data shallower than 1000 m
  ixb <- tomap$x>=min(bathy.08$lonvec) & tomap$x<=max(bathy.08$lonvec)
  iyb <- tomap$y>=min(bathy.08$latvec) & tomap$y<=max(bathy.08$latvec)
  tomap$x <- tomap$x[ixb]; tomap$y <- tomap$y[iyb]; tomap$z <- tomap$z[ixb,iyb]
  tomap$z[bathy.08$zmat>-deptht[46,1]] <- NA
  
  par(mar = c(7,1,0,1))
  image(tomap$x, tomap$y, tomap$z, xlim = rlon, ylim = rlat, zlim = c(-0.1,n.diff/10),
        xaxt = "n", yaxt = "n", xlab="", ylab="", col = col.diff)
  plot(coastline, clon = 0, clat = 0, col = "black", bg = rgb(1,1,1,0), lwd = 0.5, axes = F, add = T)
  contour(bathy.5$lonvec,bathy.5$latvec,bathy.5$zmat, levels = c(-1000,-2000), lwd = 0.5, lty = 1, labcex = 1, col = "tan4", drawlabels = F, add = T)
  contour(mldmap$x, mldmap$y, mldmap$z, levels = c(1000), lwd = c(2), lty = 1, cex = 1, col = col.periods[2], drawlabels = F, add = T)
  lines(plotIRLAB, lwd = 1, lty = 3, col = "goldenrod")
  text(x = -24, y = 51, labels = "Str3", col = col.periods[2], cex = 1.4, font = 2)
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #               Plot spatial statistics: monthly data at 1000 m and annual vertical profiles
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  expref <- "REF2" # choose here for panels c (x3 subpanels), d
  expalt <- "M2hD3h" # choose here for panels c (x3 subpanels), d
  jz <- 46 # 1000m index
  pz <- 24:55 # 97m to +2000m index
  jalpha <- 30
  x <- seq(1,12)
  zprof <- deptht[pz,1]
  f_rcat <- function(v1, v2) {c(v1, rev(v2))} # concatenate vector v1 and reversed v2 to close polygon (for shaded envelopes)
  
  
  # -----------------------------------------------------------------------------------------------------
  # Fig. 4C-left. POC_pkt diffusion
  vv <- "zdfpkt"
  ydummy <- Fquant$ymonmean_2009_2013[[expref]][[vv]][["q50"]][jz,]
  yrange <- list()
  for (yy in c("2009_2013","2014_2018")) {
    yrange[[paste0("ymonmean_",yy)]] <- f_rcat(Fquant[[paste0("ymonmean_",yy)]][[expref]][[vv]][["q25"]][jz,],
                                               Fquant[[paste0("ymonmean_",yy)]][[expref]][[vv]][["q75"]][jz,])
  }
  par(mar = c(3,6,1,0))
  plot(x, ydummy*convfact, pch=20, cex=0.2, col="black", ylim = c(0,10.1), #ylim = range(unlist(yrange)*convfact, na.rm=T),
       xlab = "", ylab = "",
       cex.lab = 1.3, las=1, lty = 1, bty = "n", xaxt = "n", yaxt = "n", yaxs = "i")
  grid(lwd = 0.5, col = "lightgray")
  for (yy in c("2009_2013","2014_2018")) {
    polygon(y = yrange[[paste0("ymonmean_",yy)]]*convfact, ylab = "",
            x = f_rcat(x,x), border=FALSE, col=paste0(col.periods[yy],jalpha), add = T)
    ymedref <- Fquant[[paste0("ymonmean_",yy)]][[expref]][[vv]][["q50"]][jz,]
    lines(x, ymedref*convfact, lwd=2, col=col.periods[yy])
    ymedalt <- Fquant[[paste0("ymonmean_",yy)]][[expalt]][[vv]][["q50"]][jz,]
    lines(x, ymedalt*convfact, lwd=2, lty = 3, col=col.periods[yy])
  }
  axis(side = 1, at = seq(3,12,3),
       labels = yday(c("2015-03-01", "2015-06-01", "2015-09-01", "2015-12-01")),
       cex.axis = 1.2, lwd = 0.5, tck=-0.03)
  axis(side = 1, at = x, labels = rep("", length(x)), cex.axis = 1.2, lwd = 0.5, tck=-0.01)
  axis(side = 2, at = seq(0,10,2), labels = seq(0,10,2), cex.axis = 1.2, lwd = 0.5, tck=-0.03, las = 1)
  mtext(side = 2, cex = 0.8, line = 3, expression(paste("Flux (mg C ",m^-2," ",d^-1,")")))
  mtext(side = 3, adj = -0.5, "C", line = 1, cex = 1.4, font = 2)
  mtext(side = 3, expression(paste("Mixing, POC"[plankton])), cex = 0.9, line = 0.5, adj = 0.5)
  
  # Legend (common to all C subpanels)
  legend(x = 1, y = 8.5, cex = 1,  bg = "#FFFFFF95", box.col = "#1a1a1a10",
         legend = rep("",2),
         col = c(col.periods[1], col.periods[2]),
         lty = rep(1,3),  lwd = rep(1.5,3), seg.len = rep(2,3))
  legend(x = 4.5, y = 8.5, cex = 1,  bg = "#FFFFFF95", box.col = "#1a1a1a10",
         legend = c("Wk3","Str3"),
         col = c(col.periods[1], col.periods[2]),
         lty = rep(3,3),  lwd = rep(2,2), seg.len = rep(2,2))
  text(x = c(2.5,6), y = 9, labels = c("REF","SAT"), cex = 1)
  
  # -----------------------------------------------------------------------------------------------------
  # Fig. 4C-middle. POC_det diffusion
  vv <- "zdfdetoc"
  ydummy <- Fquant$ymonmean_2009_2013[[expref]][[vv]][["q50"]][jz,]
  yrange <- list()
  for (yy in c("2009_2013","2014_2018")) {
    yrange[[paste0("ymonmean_",yy)]] <- f_rcat(Fquant[[paste0("ymonmean_",yy)]][[expref]][[vv]][["q25"]][jz,],
                                               Fquant[[paste0("ymonmean_",yy)]][[expref]][[vv]][["q75"]][jz,])
  }
  par(mar = c(3,2,1,0))
  plot(x, ydummy*convfact, pch=20, cex=0.2, col="black", ylim = c(0,10.1), #ylim = range(unlist(yrange)*convfact, na.rm=T),
       xlab="", ylab="", cex.lab = 1.3, las=1, lty = 1, bty = "n", xaxt = "n", yaxt = "n", yaxs = "i")
  grid(lwd = 0.5, col = "lightgray")
  for (yy in c("2009_2013","2014_2018")) {
    polygon(y = yrange[[paste0("ymonmean_",yy)]]*convfact,
            x = f_rcat(x,x), border=FALSE, col=paste0(col.periods[yy],jalpha), add = T)
    ymedref <- Fquant[[paste0("ymonmean_",yy)]][[expref]][[vv]][["q50"]][jz,]
    lines(x, ymedref*convfact, lwd=2, col=col.periods[yy])
    ymedalt <- Fquant[[paste0("ymonmean_",yy)]][[expalt]][[vv]][["q50"]][jz,]
    lines(x, ymedalt*convfact, lwd=2, lty = 3, col=col.periods[yy])
  }
  # axis(side = 1, at = seq(3,12,3), labels = seq(3,12,3), cex.axis = 1.2, lwd = 0.5, tck=-0.03)
  axis(side = 1, at = seq(3,12,3),
       labels = yday(c("2015-03-01", "2015-06-01", "2015-09-01", "2015-12-01")),
       cex.axis = 1.2, lwd = 0.5, tck=-0.03)
  axis(side = 1, at = x, labels = rep("", length(x)), cex.axis = 1.2, lwd = 0.5, tck=-0.01)
  axis(side = 2, at = seq(0,10,2), labels = seq(0,10,2), cex.axis = 1.2, lwd = 0.5, tck=-0.03, las = 1)
  mtext(side = 3, expression(paste("Mixing, POC"[detritus])), cex = 0.9, line = 0.5, adj = 0.5)
  mtext(side = 1, cex = 0.8, line = 2, "Day of year")
  
  # -----------------------------------------------------------------------------------------------------
  # Fig. 4C-right. POC_det sinking
  vv <- "exppoc"
  ydummy <- Fquant$ymonmean_2009_2013[[expref]][[vv]][["q50"]][jz,]
  yrange <- list()
  for (yy in c("2009_2013","2014_2018")) {
    yrange[[paste0("ymonmean_",yy)]] <- f_rcat(Fquant[[paste0("ymonmean_",yy)]][[expref]][[vv]][["q25"]][jz,],
                                               Fquant[[paste0("ymonmean_",yy)]][[expref]][[vv]][["q75"]][jz,])
  }
  par(mar = c(3,2,1,4))
  plot(x, ydummy*convfact, pch=20, cex=0.2, col="black", ylim = c(0,25), #ylim = range(unlist(yrange)*convfact, na.rm=T),
       xlab="", ylab="", cex.lab = 1.3, las=1, lty = 3, bty = "n", xaxt = "n", yaxt = "n", yaxs = "i")
  grid(lwd = 0.5, col = "lightgray")
  for (yy in c("2009_2013","2014_2018")) {
    polygon(y = yrange[[paste0("ymonmean_",yy)]]*convfact,
            x = f_rcat(x,x), border=FALSE, col=paste0(col.periods[yy],jalpha), add = T)
    ymedref <- Fquant[[paste0("ymonmean_",yy)]][[expref]][[vv]][["q50"]][jz,]
    lines(x, ymedref*convfact, lwd=2, col=col.periods[yy])
    ymedalt <- Fquant[[paste0("ymonmean_",yy)]][[expalt]][[vv]][["q50"]][jz,]
    lines(x, ymedalt*convfact, lwd=2, lty = 3, col=col.periods[yy])
  }
  axis(side = 1, at = seq(3,12,3),
       labels = yday(c("2015-03-01", "2015-06-01", "2015-09-01", "2015-12-01")),
       cex.axis = 1.2, lwd = 0.5, tck=-0.03)
  axis(side = 1, at = x, labels = rep("",length(x)), cex.axis = 1.2, lwd = 0.5, tck=-0.01)
  axis(side = 2, at = seq(0,25,5), labels = seq(0,25,5), cex.axis = 1.2, lwd = 0.5, tck=-0.03, las = 1)
  mtext(side = 3, expression(paste("Grav. sinking, POC"[detritus])), cex = 0.9, line = 0.5, adj = 0.5)
  
  # -----------------------------------------------------------------------------------------------------
  # Fig. 4D. Vertical profiles of total diffusive fluxes
  vv <- "zdfpoc"
  xdummy <- Fquant$yearmean_2009_2013[[expref]][[vv]][["q50"]][pz]
  xrange <- list()
  for (yy in c("2009_2013","2014_2018")) {
    xrange[[paste0("yearmean_",yy)]] <- f_rcat(Fquant[[paste0("yearmean_",yy)]][[expref]][[vv]][["q25"]][pz],
                                               Fquant[[paste0("yearmean_",yy)]][[expref]][[vv]][["q75"]][pz])
  }
  
  par(mar = c(3,2,1,1))
  xmax <- 25
  plot(xdummy*convfact, zprof, pch=20, cex=0.2, col="black", bty = "n", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i",
       xlim = c(-1,xmax), #range(unlist(xrange)*convfact, na.rm=T),
       xlab="",
       ylim = c(2100,0), # c(max(zprof),min(zprof)),
       ylab="Depth (m)", cex.lab = 1.3, las=1, lty = 1)
  grid(lwd = 0.5, col = "lightgray")
  for (yy in c("2009_2013","2014_2018")) {
    
    # Fluxes
    polygon(x = xrange[[paste0("yearmean_",yy)]]*convfact,
            y = f_rcat(zprof,zprof), border=FALSE, col=paste0(col.periods[yy],jalpha), add = T)
    xmedref <- Fquant[[paste0("yearmean_",yy)]][[expref]][[vv]][["q50"]][pz]
    lines(xmedref*convfact, zprof, lwd=2, col=col.periods[yy])
    xmedalt <- Fquant[[paste0("yearmean_",yy)]][[expalt]][[vv]][["q50"]][pz]
    lines(xmedalt*convfact, zprof, lwd=2, lty=3, col=col.periods[yy])
  }
  # abline(h = c(1000,2000), col = "darkgray", lty = 1, lwd = 0.5) # replaced by visual guides of mean convection depths (above)
  axis(side = 1, at = seq(0,25,5), labels = seq(0,25,5), cex.axis = 1.2, lwd = 0.5, tck=-0.03, las = 1)
  axis(side = 3, at = seq(0,20,5), labels = rep("",5), cex.axis = 1.2, lwd = 0.5, tck=-0.01, las = 1)
  axis(side = 2, at = c(2000,1500,1000,500,100), cex.axis = 1.2, lwd = 0.5, tck=-0.03, las = 1)
  mtext(side = 2, cex = 0.8, line = 4, text = "Depth")
  mtext(side = 1, cex = 0.8, line = 2.2, expression(paste("Flux (mg C ",m^-2," ",d^-1,")")))
  # mtext(side = 3, adj = 0, "d)", line = 1)
  mtext(side = 3, adj = 0, "D ", line = 1, cex = 1.4, font = 2)
  mtext(side = 3, adj = 0.5, expression(paste("Mixing, POC"[total])), cex = 0.9, line = 0.8)
  box(lwd = 0.5)
  
  # -----------------------------------------------------------------------------------------------------
  # Colorbar(s): with this method they have to be called after all plot windows have been filled
  
  # Diffusion contribution to vertical fluxes (panel B)
  par( fig=c(.75,.985,.50,.52), new=TRUE, mar=c(0,0,0,0) ) # with bottom panels nrow=9
  cbar.maps <- seq(0, n.diff/10, length.out = n.diff+1)
  cbar.tick <- seq(-0.1, n.diff/10, 0.1)
  image(cbar.maps, c(0,1), matrix(cbar.maps, c(n.diff+1,2)), col = col.diff, xaxt = "n", yaxt = "n", bg = "white")
  axis(1, cex.axis=1.2, mgp = c(0, 0.6, 0), tck = -0.4, at = cbar.tick+0.05, labels = cbar.tick*100, las = 1, lwd = 0.5)
  mtext(side = 1, las = 1, line = 2, cex = 0.8,
        text = expression(paste("Mix"[Z]," / ","[Mix"[Z]," + Grav"[],"] (%)")) )
  box(lwd = 0.5)
  
  dev.off()
}



# -----------------------------------------------------------------------------------------------------
# Map MLD contours and POC section shown in SM. NOTE: this figure depends on fields computed in fig_convfluxes (just above)

if (sfig_sections_map & fig_convfluxes) {
  
  mldmap.w <- f_prep4contour(xin = nav_lon[vcell], yin = nav_lat[vcell],
                             zin = Amlotstmax.2009_2013,
                             dlondeg = 0.5, dlatdeg = 0.5)
  mldmap.s <- f_prep4contour(xin = nav_lon[vcell], yin = nav_lat[vcell],
                             zin = Amlotstmax.2014_2018,
                             dlondeg = 0.5, dlatdeg = 0.5)
  load("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/section_masks_NASPGcells_orca1.Rda")
  
  png(filename = "~/Desktop/Gali_2026_convectionPOC/output/Fig_S16A_section.png",
      width = 10, height = 12, units = 'cm', pointsize = 16,
      bg = 'white', res = 600, type = 'cairo')
  
  par(mar = c(4,4,1,1))
  image(bathy.08$lonvec, bathy.08$latvec, bathy.08$zmat, xlab = "Lon (ºE)", ylab = "Lat (ºN)",
        xlim = c(-64,-16), ylim = c(40,70),
        col = oce.colorsGebco(n = 100, region = "water")[51:100],
        mgp = c(2, 0.2, 0), tck=-0.01, cex.axis = 1, cex.lab = 1.2, main = "")
  
  plot(coastline, clon = 0, clat = 0, col = "black", bg = rgb(1,1,1,0), lwd = 0.5, axes = F, add = T)
  contour(bathy.08$lonvec,bathy.08$latvec,bathy.08$zmat, levels = c(-1000,-2000), lwd = 0.5, lty = 1, labcex = 1, col = "tan4", drawlabels = F, add = T)
  contour(mldmap.w$x, mldmap.w$y, mldmap.w$z, levels = c(1000), lwd = 3, lty = 1, cex = 1, col = col.periods[1], drawlabels = F, add = T)
  contour(mldmap.s$x, mldmap.s$y, mldmap.s$z, levels = c(1000), lwd = 3, lty = 1, cex = 1, col = col.periods[2], drawlabels = F, add = T)
  points(sectionmask$LABNW$xmat[sectionmask$LABNW$zmat==1],
         sectionmask$LABNW$ymat[sectionmask$LABNW$zmat==1],
         col = "orange1", pch = 20, cex = 1)
  
  dev.off()
  
}



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Calculate some temporal correlations and linear regressions for the text
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

expname <- "M2hD3h" # REF2, M1hD1h, M2hD3h

# Areal fluxes versus convection extent (deep convection area)
zind <- 3 # 1 is 100, 2 is 500, 3 is 1000, 4 is 2000
jy <- which(seq(yearS[[expname]][1],yearS[[expname]][2]) %in% seq(1958,2019))
x <- df.opera.yrarea$mlotst_max_1000[jy + yearS[[expname]][1]-1958] * 1e-6
y1 <- F_ztsum[[expname]][["zdfpoc"]][zind,jy]
y2 <- F_ztsum[[expname]][["zafpoc"]][zind,jy]
y3 <- F_ztsum[[expname]][["exppoc"]][zind,jy]
y <- y1 # + y2 + y3
print( cor.test(x, y1)[c(4,9)] ) # print r and 05% CI only
print( cor.test(x, y2)[c(4,9)] )
print( cor.test(x, y3)[c(4,9)] )
print( cor.test(x, y)[c(4,9)] )
linfit <- lm(y ~ x, data = data.frame(x, y))
cifit <- confint(linfit)
print( linfit$coefficients )
print( (cifit[,2]-cifit[,1])/2 ) # 95% confidence intervals
print( summary.lm(linfit) )      # R2adj


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# -------------------------------------------------------- FIG. 5 -------------------------------------------------------

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (scheme_budgets) {
  
  # Read pre-processed budget data in xls files generated by Andrea (MAOG). Unit conversion to Tg C yr-1 using budgetarea.IR_LAB
  # NOTE: zingest2doc available only for period 2009-2018 and experiments REF2 (a67o) and M2hD3h = SAT restoring (a683)
  xlsrange <- "A1:X23"
  convfactor <- 86400 * 365 * 12 * budgetarea.IR_LAB * 1e6 / 1e12
  budget <- list(
    REF2 = (
      read_excel(
        path = "~/Desktop/Gali_2026_convectionPOC/input_data/fig_4_5_S14-16_nemo-pisces/all_climatologies_1998-2020_a67o.xlsx",
        sheet = "lower_meso_500-1000m",
        range = xlsrange
      ) +
        read_excel(
          path = "~/Desktop/Gali_2026_convectionPOC/input_data/fig_4_5_S14-16_nemo-pisces/all_climatologies_1998-2020_a67o.xlsx",
          sheet = "bathy_1000-2000m",
          range = xlsrange
        )
    ) * convfactor,
    M2hD3h = (
      read_excel(
        path = "~/Desktop/Gali_2026_convectionPOC/input_data/fig_4_5_S14-16_nemo-pisces/all_climatologies_1998-2020_a683.xlsx",
        sheet = "lower_meso_500-1000m",
        range = xlsrange
      ) +
        read_excel(
          path = "~/Desktop/Gali_2026_convectionPOC/input_data/fig_4_5_S14-16_nemo-pisces/all_climatologies_1998-2020_a683.xlsx",
          sheet = "bathy_1000-2000m",
          range = xlsrange
        )
    ) * convfactor
  )
  
  budget <- lapply(budget, function(x) {
    x[["tottrdpoc"]] <- Reduce('+', x[grep("tottrd", names(x), value = T)]) # total trend
    x[["xadpoc"]] <- Reduce('+', x[grep("xad", names(x), value = T)]) # add xadpoc
    x[["yadpoc"]] <- Reduce('+', x[grep("yad", names(x), value = T)]) # add yadpoc
    x[["zadpoc"]] <- Reduce('+', x[grep("zad", names(x), value = T)]) # add zadpoc
    x[["ldfpoc"]] <- Reduce('+', x[grep("ldf", names(x), value = T)]) # add ldfpoc
    x[["lat.trans"]] <- x[["xadpoc"]]+x[["yadpoc"]]+x[["ldfpoc"]] # add lateral transport
    x[["zdfpoc"]] <- Reduce('+', x[grep("zdf", names(x), value = T)]) # add zdfpoc
    x[["exppoc"]] <- Reduce('+', x[grep("trdexp", names(x), value = T)]) # add exppoc
    x[["rempoc"]] <- Reduce('+', x[grep("rem", names(x), value = T)]) # add remin. SINK WITH POSITIVE SIGN (zingest2doc AS WELL)
    x[["partialbalance"]] <- x[["lat.trans"]]+x[["zadpoc"]]+x[["zdfpoc"]]+x[["exppoc"]]-x[["rempoc"]]-x[["zingest2doc"]]
    x[["doc2poc"]] <- x[["tottrdpoc"]] - x[["partialbalance"]] # residual, equal to doc2poc
    x[["poc2doc"]] <- -x[["rempoc"]]-x[["zingest2doc"]]+x[["doc2poc"]]
    return(x)
  })
  for (nn in names(budget)) {
    budget[[nn]][["Year"]] <- seq(1998,2019)
    budget[[nn]][["DCV"]] <- S_ztsum$REF2$mlotstmax_500_2000[(62-21):62]*1e-15
  }
  
  # Test differences for a given variable and experiment across periods
  tx <- "REF2" # REF2 M1hD1h M2hD3h
  tv <- "tottrdpoc"
  t.test(x = budget[[tx]][[tv]][budget[[tx]]$Year%in%allperiods$weak3], y = budget[[tx]][[tv]][budget[[tx]]$Year%in%allperiods$strong3])
  
  # Statistics for Fig. 5 ----
  # POC diffusion, zdfpoc
  # Difference between REF and SAT exp on a given period: only during strong mixing
  t.test(x = budget$REF2$zdfpoc[budget$REF2$Year%in%allperiods$weak3], y = budget$M2hD3h$zdfpoc[budget$M2hD3h$Year%in%allperiods$weak3])
  t.test(x = budget$REF2$zdfpoc[budget$REF2$Year%in%allperiods$strong3], y = budget$M2hD3h$zdfpoc[budget$M2hD3h$Year%in%allperiods$strong3])
  # Difference between periods for REF and SAT, zdfpoc
  t.test(x = budget$REF2$zdfpoc[budget$REF2$Year%in%allperiods$weak3], y = budget$REF2$zdfpoc[budget$REF2$Year%in%allperiods$strong3])
  t.test(x = budget$M2hD3h$zdfpoc[budget$M2hD3h$Year%in%allperiods$weak3], y = budget$M2hD3h$zdfpoc[budget$M2hD3h$Year%in%allperiods$strong3])
  
  # POC grav export, exppoc
  # Difference between REF and SAT exp on a given period: REF 17% larger during all periods
  t.test(x = budget$REF2$exppoc[budget$REF2$Year%in%allperiods$weak3], y = budget$M2hD3h$exppoc[budget$M2hD3h$Year%in%allperiods$weak3])
  t.test(x = budget$REF2$exppoc[budget$REF2$Year%in%allperiods$strong3], y = budget$M2hD3h$exppoc[budget$M2hD3h$Year%in%allperiods$strong3])
  # Difference between periods for REF and SAT
  t.test(x = budget$REF2$exppoc[budget$REF2$Year%in%allperiods$weak3], y = budget$REF2$exppoc[budget$REF2$Year%in%allperiods$strong3])
  t.test(x = budget$M2hD3h$exppoc[budget$M2hD3h$Year%in%allperiods$weak3], y = budget$M2hD3h$exppoc[budget$M2hD3h$Year%in%allperiods$strong3])
  
  # POC to DOC
  # Difference between REF and SAT exp for poc2doc on a given period: SUMMARY
  t.test(x = budget$REF2$poc2doc[budget$REF2$Year%in%allperiods$weak3], y = budget$M2hD3h$poc2doc[budget$M2hD3h$Year%in%allperiods$weak3])
  t.test(x = budget$REF2$poc2doc[budget$REF2$Year%in%allperiods$strong3], y = budget$M2hD3h$poc2doc[budget$M2hD3h$Year%in%allperiods$strong3])
  # Difference between periods for REF and SAT, exppoc
  t.test(x = budget$REF2$poc2doc[budget$REF2$Year%in%allperiods$weak3], y = budget$REF2$poc2doc[budget$REF2$Year%in%allperiods$strong3])
  t.test(x = budget$M2hD3h$poc2doc[budget$M2hD3h$Year%in%allperiods$weak3], y = budget$M2hD3h$poc2doc[budget$M2hD3h$Year%in%allperiods$strong3])
  
  # CO2 flux (a67o only)
  S_fgco2 <- f_load_2D_variables("a67o", mbasepath, dirname = "year", "fgco2", search_pattern = "12.nc")
  S_fgco2 <- f_xyzt_xysum_mask(Lxyzt = S_fgco2, MASK = smask[[regname]])
  fgco2 <- data.frame(Year = 1958:2019,
                      fgco2 = unlist(S_fgco2$fgco2*86400*365*1e-9)) # from kg/s to TgC/yr
  t.test(x = fgco2$fgco2[fgco2$Year%in%allperiods$weak3], y = fgco2$fgco2[fgco2$Year%in%allperiods$strong3])
  
  
  fit_dcv <- lapply(budget, function(x) {
    linfit <- lm(zdfpoc ~ DCV, data = x)
    cifit <- t(confint(linfit))
    names(cifit) <- c("Intercept_025","Intercept_975","Slope_025","Slope_975")
    # print( linfit$coefficients )
    # print( (cifit[,2]-cifit[,1])/2 ) # 95% confidence intervals
    # print( summary.lm(linfit) )      # R2adj
    return(c(linfit, cifit,summary.lm(linfit)))
  })
  # View(fit_dcv)
  lapply(fit_dcv, function(x) {
    print(x$coefficients)
    print((x$Intercept_975 - x$Intercept_025)/2)
    print((x$Slope_975 - x$Slope_025)/2)
  })
  
  # Arrange in long format for plotting and further fitting and hypothesis testing
  dcvplot <- data.table::rbindlist(budget, idcol = "expid") %>% filter(expname != "M1hD1h")
  
  # Alternative fitting and hypothesis testing
  dcvplot$expid <- as.factor(dcvplot$expid) # Ensure expid is a factor
  model <- lm(zdfpoc ~ DCV * expid, data = dcvplot) # Fit model and test for slope differences
  summary(model)
  
  # Compare full model (different slopes & intercepts) vs reduced model (same line)
  reduced_model <- lm(zdfpoc ~ DCV + expid, data = dcvplot)  # same slope, different intercepts
  anova(reduced_model, model)  # Tests if allowing slopes to differ improves fit
  
  # Get the interaction p-value
  interaction_p <- coef(summary(model))[grep(":", rownames(coef(summary(model)))), "Pr(>|t|)"]
  cat("P-value for slope difference:", interaction_p, "\n")
  
  # Additional checks
  shapiro.test(residuals(model)) # Normality of residuals: yes, p almost 1
  bartlett.test(residuals(model) ~ dcvplot$expid) # Homogeneity of variance: yes, p = 0.13
  
  # install.packages("emmeans")
  library(emmeans)
  emtrends(model, ~ expid, var = "DCV")
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # --------------------------------------------------- Fig. 5A ---------------------------------------------------
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if (fig_dcv_budget) {
    
    p <- "~/Desktop/Gali_2026_convectionPOC/output/Fig_5A_DCV_budget.png"
    png(filename = p, width = 6.5, height = 7, units = 'cm', pointsize = 6, bg = "white", res = 600, type = "quartz")
    
    par(mar = c(4,5,4,2))
    par(xpd = T)
    
    plot(x = dcvplot$DCV,
         y = dcvplot$zdfpoc,
         type = "p", pch = 1, cex = 0, las = 1,
         bty = "n",
         xaxs = "i", yaxs = "i",
         cex.axis = 1.2,
         xlim = c(0,1.2), # edit xmax if y boxplot added on margins
         ylim = c(0,8), # edit ymax if x boxplot added on margins 
         cex.lab = 1.4,
         xlab = expression(paste("Deep convection volume (",10^15," ",m^3,")")),
         ylab = expression(paste("Mixing-driven POC supply (Tg C ",yr^-1,")")))
    points(dcvplot$DCV[dcvplot$expid=="REF2"], dcvplot$zdfpoc[dcvplot$expid=="REF2"], pch = 16, col = "black", cex = 2)
    points(dcvplot$DCV[dcvplot$expid=="REF2"&dcvplot$Year%in%allperiods$weak3], dcvplot$zdfpoc[dcvplot$expid=="REF2"&dcvplot$Year%in%allperiods$weak3], pch = 1, col = col.periods[1], cex = 2, lwd=1.6)
    points(dcvplot$DCV[dcvplot$expid=="REF2"&dcvplot$Year%in%allperiods$strong3], dcvplot$zdfpoc[dcvplot$expid=="REF2"&dcvplot$Year%in%allperiods$strong3], pch = 1, col = col.periods[2], cex = 2, lwd=1.6)
    points(dcvplot$DCV[dcvplot$expid=="M2hD3h"], dcvplot$zdfpoc[dcvplot$expid=="M2hD3h"], pch = 16, col = "orange1", cex = 2)
    points(dcvplot$DCV[dcvplot$expid=="M2hD3h"&dcvplot$Year%in%allperiods$weak3], dcvplot$zdfpoc[dcvplot$expid=="M2hD3h"&dcvplot$Year%in%allperiods$weak3], pch = 1, col = col.periods[1], cex = 2, lwd=1.6)
    points(dcvplot$DCV[dcvplot$expid=="M2hD3h"&dcvplot$Year%in%allperiods$strong3], dcvplot$zdfpoc[dcvplot$expid=="M2hD3h"&dcvplot$Year%in%allperiods$strong3], pch = 1, col = col.periods[2], cex = 2, lwd=1.6)
    # points(dcvplot$DCV[dcvplot$expid=="M1hD1h"], dcvplot$zdfpoc[dcvplot$expid=="M1hD1h"], pch = 16, col = "orange3", cex = 2)
    # points(dcvplot$DCV[dcvplot$expid=="M1hD1h"&dcvplot$Year%in%allperiods$weak3], dcvplot$zdfpoc[dcvplot$expid=="M1hD1h"&dcvplot$Year%in%allperiods$weak3], pch = 1, col = col.periods[1], cex = 2, lwd=1.6)
    # points(dcvplot$DCV[dcvplot$expid=="M1hD1h"&dcvplot$Year%in%allperiods$strong3], dcvplot$zdfpoc[dcvplot$expid=="M1hD1h"&dcvplot$Year%in%allperiods$strong3], pch = 1, col = col.periods[2], cex = 2, lwd=1.6)
    
    lines(dcvplot$DCV[dcvplot$expid=="REF2"], fit_dcv$REF2$fitted.values, lty = 1, lwd = 2, col = "black")
    lines(dcvplot$DCV[dcvplot$expid=="M2hD3h"], fit_dcv$M2hD3h$fitted.values, lty = 1, lwd = 2, col = "orange1")
    # lines(dcvplot$DCV[dcvplot$expid=="M1hD1h"], fit_dcv$M1hD1h$fitted.values, lty = 1, lwd = 2, col = "orange3")
    
    legend(0.02, 8.0,
           legend = c("REF","SAT","Wk3 (2009-2013)","Str3 (2014-2018)"),
           bg = "#FFFFFF95", box.col = "#FFFFFF15",
           pch = c(16,16,1,1), cex = rep(1.2,4), col = c("black","orange1",col.periods))
    
    dev.off()
    
  }
  
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # ---------------------------------------- Summarize data for Fig. 5B -----------------------------------------
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # POC budgets
  summary_table <- map_dfr(names(budget), function(xx) {
    
    dfb <- budget[[xx]]
    
    numvars <- setdiff(
      names(dfb)[sapply(dfb, is.numeric)],
      "Year"
    )
    print(numvars)
    
    get_stats <- function(dfsub, period) {
      
      out <- data.frame(
        Simulation = xx,
        Period = period
      )
      for(vv in numvars) {
        out[[paste0(vv, "_mean")]] <- mean(dfsub[[vv]], na.rm = TRUE)
        out[[paste0(vv, "_sd")]]   <- sd(dfsub[[vv]], na.rm = TRUE)
      }
      out
    }
    bind_rows(
      get_stats(filter(dfb, Year %in% allperiods$weak3), "2009-2013"),
      get_stats(filter(dfb, Year %in% allperiods$strong3), "2014-2018")
    )
  })
  
  View(summary_table) # input data for Fig. 5B
  
}

# END OF SCRIPT


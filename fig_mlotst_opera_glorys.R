# Compare (plot) the extent of deep convection in OPERA vs. GLORYS12
# Follows structure of fig_conv_fluxes.R
# Regions defined in map_profile_counts.NASPGcoriolis.R according to convection area, 3000 m isobath and BGC-Argo profile availability
# Functions for processing arrays in list inherited from hovmoller_pisces_budgetPOC_3D.R

# Libraries
library(RNetCDF)
library(ggplot2)
library(plyr) # for mapvalues function
library(dplyr) # for select function
library(tidyr) # for pivot_wider
library(reshape) # for rename function
library(data.table) # for rbindlist
library(ggh4x) # for "stat_difference" in ggplot2
library(oce)
library(ocedata)
library(akima) # for interp function
library(RColorBrewer)
library(lubridate)
library(fields)
library(areaplot)
source("~/Desktop/Gali_2026_convectionPOC/f_skillstats_xvec_yvec.R")

# -----------------------------------------------------------------------------------------------------
# Function to load all desired variables from a given experiment (expid)
# all fluxes are defined + downwards

f_load_2D_variables <- function(expid, mbasepath, varnames) {
  datalist <- lapply(varnames, function(vv) {
    fname <- grep(pattern = vv,
                  x = list.files(paste0(mbasepath, expid, "/", timeres), pattern = ".nc", full.names = T),
                  value = T)
    ncfile <- open.nc(fname)
    return( var.get.nc(ncfile, variable = vv) )
    close.nc(ncfile)
  })
  names(datalist) <- varnames
  return(datalist)
}

# -----------------------------------------------------------------------------------------------------
# Settings, definitions

# Define list of experiments and analysis period
expidS <- list(
  # REF1 = "a5gj",     # Old reference, no restoring, 1958-2019:                       expid <- "a5gj", expdate <- "v20221216"
  REF2 = "a67o"        # New reference, no restoring, 1958-2019:                       expid <- "a67o", expdate <- "v20230630"
)
regname4mean <- "IR_LAB" # NOTE: "LAB_CB" not usable for GLORYS12v1 data processed offline for IR_LAB domain (by María SU)
regname4area <- "NASPG"  # NOTE: IR_LAB currently not available for mlotst_max area above selected thresholds
timeres <- "month" # use month, then process OPERA and GLORYS jointly to get annual metrics
season <- NULL # NULL. Seasonal options "_FMAM", "_JJASO" not enabled here
startyear <- 1958  # first year of OPERA REF simulation
endyear <- 2019 
varnames.phys2D <- c("omlda","omldamax","mlotst","mlotstmax")
newnames.phys2D <- c("omlda_mean","omlda_max","mlotst_mean","mlotst_max")

# Base paths
mbasepath <- "~/Desktop/Gali_2026_convectionPOC/input_data/fig_4_5_S14-16_nemo-pisces/"

# Load lists of selected cells (inside polygon), Argo profile counts, etc: produced from "map_profile_counts.NASPGcoriolis.R"
NASPGij <- read.csv("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/NASPGij_orca1.csv")
load(paste0("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/inpolygon_NASPGcells_orca1.Rda"))
load(paste0("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/inpolygon_labCBcells_orca1.Rda"))
load(paste0("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/inpolygon_IRLABcells_orca1.Rda"))

# Load ORCA1 horizontal grid, areacello, volcello
load("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/ORCA1_hgrid.Rda")
anc <- open.nc("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/areacello_Ofx_EC-Earth3_opera-control_r1i1p1f1_gn.nc")
areacello <- var.get.nc(anc, "areacello", start = c(NASPGij$istart, NASPGij$jstart), count = c(NASPGij$icounts, NASPGij$jcounts)) # ORCA1 areacello (m2!)
close.nc(anc)

# Time axis
ifelse(timeres == "month",
       decdate <- seq(startyear,endyear+11/12,1/12) + 0.5/12,
       decdate <- seq(startyear,endyear,1) + 0.5/12)
date <- floor_date(date_decimal(decdate, tz = "GMT"), timeres) + 0.5*as.duration(timeres)

# Load global bathymetry mask and subset according to NASPGcells (rectangular lon-lat domain)
# Using 2500 m for consistency with POC budget domain. Consistent results with 2000 m mask.
mask2500 <- var.get.nc(ncfile = open.nc("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/bathymask_2500_0_1.nc"),
                       variable = "bathymask",
                       start = c(NASPGij$istart, NASPGij$jstart),
                       count = c(NASPGij$icounts, NASPGij$jcounts))

# Create rectangular masks of labCB and IRLAB domains based on cell indices and bathymetry
basemask <- 0*nav_lat
f_makemask <- function(regbounds, maskcells) {
  basemask[maskcells$imatch] <- 1
  maskout <- basemask[regbounds$istart:(regbounds$istart+regbounds$icount-1),regbounds$jstart:(regbounds$jstart+regbounds$jcount-1)]
  return( maskout & mask2500 )
}
smask <- list(
  LAB_CB = f_makemask(NASPGij, labCBcells),
  IR_LAB = f_makemask(NASPGij, IRLABcells)
)
smask$NASPG <- array(TRUE, dim = dim(smask$IR_LAB))


# -----------------------------------------------------------------------------------------------------
# Load various MLD metrics from OPERA simulation(s)
Spre <- list()
for (nn in names(expidS)) {
  expid <- expidS[[nn]]
  Spre[[nn]] <- f_load_2D_variables(expid, mbasepath, varnames.phys2D)
}

# Rename opera variables by adding "_" before the mean or max stat
Spre <- lapply(Spre, function(x) {names(x) <- newnames.phys2D; return(x)})

# -----------------------------------------------------------------------------------------------------
# Compute spatial means over selected domain (IR_LAB)

# Define function to compute horizontal means for list of 4D (xyzt) or 3D (xyt) arrays, including masks
f_xyzt_xymean_mask <- function(Lxyzt, Vxyzt, MASK) {
  
  if (!is.null(MASK)) {
    Vxyzt <- Vxyzt * array(MASK, dim(Vxyzt))
  }
  if (length(dim(Vxyzt))==4) {
    amargin <- c(3,4)
  } else if (length(dim(Vxyzt))==3) {
    amargin <- 3
  }
  lapply(Lxyzt, function(A) {
    AxV <- A * Vxyzt # element by element product
    AxVsum <- apply(AxV, MARGIN = amargin, sum, na.rm=T) # handling nans properly?
    Vsum <- apply(Vxyzt, MARGIN = amargin, sum, na.rm=T) # handling nans properly?
    return(AxVsum / Vsum)
  })
}
S_zt <- lapply(Spre, f_xyzt_xymean_mask,
               array(areacello, dim = dim(Spre$REF2[[1]])),
               smask[[regname4mean]])

# Function to compute horizontal sums for list of 3D (xyt) or 4D (xyzt) arrays, including spatial masks
f_xyzt_xysum_mask <- function(Lxyzt, MASK) {
  
  Sxy <- areacello # assuming areacello has been loaded. TO DO: PUT AREACELLO IN Spre LIST
  
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


# -----------------------------------------------------------------------------------------------------
# Compute area where MLD metric exceeds a given z (500, 1000, 2000) over selected domain (IR_LAB)

# Compute area (km2) exceeding MLD metric threshold "Zt" for list of 3D (xyt) arrays (monthly)
f_xyt_tarea_mask <- function(Lxyt, Sxy, Zt, MASK) {
  
  # Apply mask if specified
  if (!is.null(MASK)) {
    Sxyt <- array( (MASK * Sxy) , dim(Lxyt[[1]]))
  }
  
  # Rename columns
  names(Lxyt) <- paste0(names(Lxyt),"_",Zt)
  
  # Compute area above Zt threshold
  lapply(Lxyt, function(A) {
    A[A<=Zt] <- NA # check NA treated properly by comparing results with either 0 or NA
    A[A>Zt] <- 1
    AxS <- A * Sxyt # element by element product
    AxSsum <- apply(AxS, MARGIN = 3, sum, na.rm=T) / 1e6
    return(AxSsum)
  })
}
S_tarea_mo_500 <- lapply(Spre, f_xyt_tarea_mask, areacello, 500, smask[[regname4area]])
S_tarea_mo_1000 <- lapply(Spre, f_xyt_tarea_mask, areacello, 1000, smask[[regname4area]])
S_tarea_mo_1500 <- lapply(Spre, f_xyt_tarea_mask, areacello, 1500, smask[[regname4area]])
S_tarea_mo_2000 <- lapply(Spre, f_xyt_tarea_mask, areacello, 2000, smask[[regname4area]])


# -----------------------------------------------------------------------------------------------------
# Define several 5-year periods

# Weak convection periods: 1967..1971, 1977..1981, 2009..2013
# Deep convection periods:    1972..1976, 1989..1993, 2014..2018

highlight_periods <- data.frame(
  xmin = as.Date(c('2009-01-01', '2014-01-01')), # , '1983-01-01' to add high conv period
  xmax = as.Date(c('2013-12-31', '2018-12-31')), # , '1987-12-31' to add high conv period
  ymin = -Inf,
  ymax = Inf,
  period = c("gray80", "gray40", "gray80","gray40", "gray80","gray40") # ,"gray10" to add high conv period
)



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plot multiyear time series
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# x axis (date) main and secondary ticks
dateticks1 <- seq(floor_date(min(date), unit = "10 years"), floor_date(max(date), unit = "10 years"), by = "10 years")
dateticks2 <- seq(floor_date(min(date), unit = "years"), floor_date(max(date), unit = "years"), by = "years")

# mshift <- as.duration("1 month")
mshift <- as.duration(0)
selperiod <- force_tz(as.POSIXct(c("1993-01-01","2019-12-31")), tzone = "GMT") - mshift # Period
jselperiod <- which(date >= selperiod[1] & date <= selperiod[2])

# Format as data frames the OPERA simulation MLD metrics to be compared to GLORYS12v1,
# and subset jselperiod (not indispensable because of long format, but useful to restrict temporal axis when plotting)
df.opera.momean <- cbind(
  date = format(date,"%Y-%m-%d"),
  year = year(date),
  month = month(date),
  as.data.frame(S_zt$REF2)
)[jselperiod,]

df.opera.moarea <- cbind(
  date = format(date,"%Y-%m-%d"),
  year = year(date),
  month = month(date),
  as.data.frame(S_tarea_mo_500$REF2),
  as.data.frame(S_tarea_mo_1000$REF2),
  as.data.frame(S_tarea_mo_1500$REF2),
  as.data.frame(S_tarea_mo_2000$REF2)
)[jselperiod,]


# -----------------------------------------------------------------------------------------------------
# Calculate annual maximum of monthly mlotstmax, its spatial mean over the selected domain (IR_LAB), and the area over 500, 1000 and 2000 m

# Define function to compute annual stat  (currently mean or max) from monthly stat on a list of 3D (xyt) arrays
f_xyt_mo2yr_timestat <- function(Lxyt, timestat) {
  
  # Subset variables with selected stat
  Lxyt <- Lxyt[grep(pattern = timestat, x = names(Lxyt), value=T)]
  
  lapply(Lxyt, function(A) {
    # Reshape from lon-lat-monthS to lon-lat-month-year
    orig_dims <- dim(A)
    Ar <- array(A, dim = c(orig_dims[1:2],12,orig_dims[3]/12))
    if (timestat=="mean") {
      ndaysmonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
      return( apply(Ar, c(1,2,4), weighted.mean, w=ndaysmonth, na.rm=T) )
    } else if (timestat=="max") {
      return( apply(Ar, c(1,2,4), max, na.rm=T) )
    }  else {
      print("Warning: stat not defined! NA array returned")
      return( array(NA, dim = c(orig_dims[1:2],orig_dims[3]/12)) )
    }
  })
}
S_xy_yrmax <- lapply(Spre, f_xyt_mo2yr_timestat, timestat="max")
# S_xy_yrmean <- lapply(Spre, f_xyt_mo2yr_timestat, timestat="mean") # currently unused, but tested and gives consistent picture

# Compute area (km2) exceeding MLD metric threshold "Zt" for list of 3D (xyt) arrays (annual)
S_tarea_yr_500 <- lapply(S_xy_yrmax, f_xyt_tarea_mask, areacello, 500, smask[[regname4area]])
S_tarea_yr_1000 <- lapply(S_xy_yrmax, f_xyt_tarea_mask, areacello, 1000, smask[[regname4area]])
S_tarea_yr_1500 <- lapply(S_xy_yrmax, f_xyt_tarea_mask, areacello, 1500, smask[[regname4area]])
S_tarea_yr_2000 <- lapply(S_xy_yrmax, f_xyt_tarea_mask, areacello, 2000, smask[[regname4area]])

# Compute mean of the annual maximum mlotst in the IR_LAB domain
S_zt <- lapply(S_xy_yrmax, f_xyzt_xymean_mask,
               array(areacello, dim = dim(S_xy_yrmax$REF2[[1]])),
               smask[[regname4mean]])

# Load data frames summarizing GLORYS12v1 data within the NASPG rectangular domain, pre-processed by Maria SU
df.glorys.yrarea <- read.csv(file = "~/Desktop/Gali_2026_convectionPOC/GLORYS12v1_code_data_MSU/output_data/mlotst_yearly_max_area_500_1000_1500_2000_rectangule.csv")
names(df.glorys.yrarea) <- gsub(pattern = "MLD", replacement = "mlotst_", x = names(df.glorys.yrarea))

# Format as data frames the OPERA simulation MLD metrics to be compared to GLORYS12v1 (annual)
operayears <- as.integer(unique(year(date)))
jselyears <- which( operayears %in% df.glorys.yrarea$year )

df.opera.yrarea <- cbind(
  year = operayears,
  as.data.frame(S_tarea_yr_500$REF2),
  as.data.frame(S_tarea_yr_1000$REF2),
  as.data.frame(S_tarea_yr_1500$REF2),
  as.data.frame(S_tarea_yr_2000$REF2)
)
df.opera.yrmaxmean <- cbind(
  year = operayears,
  as.data.frame(S_zt$REF2)
)
# Merge data frames (preparing to put in long format for figures)
df.opera.yrarea <- df.opera.yrarea[jselyears,] # crop to glorys period
df.yrarea <- data.table::rbindlist(list(GLORYS12v1 = df.glorys.yrarea, NEMO4_ORCA1 = df.opera.yrarea), use.names = T, fill = T, idcol = "Dataset")

# Save processed datasets, input for Fig. 4. Keep commented, datasets already available in separate folder to avoid overwriting
# save(df.opera.yrarea, df.opera.yrmaxmean, file = "~/Desktop/Gali_2026_convectionPOC/input_data/test_preprocessing/conv_extent_opera_1958_2019.Rda")
# save(df.yrarea, file = "~/Desktop/Gali_2026_convectionPOC/input_data/test_preprocessing/conv_extent_opera_glorys_1998_2019.Rda")

# Add date corresponding to middle of year (1st July) just for plotting compatibility
# https://stackoverflow.com/questions/30255833/convert-four-digit-year-values-to-class-date
df.yrarea$Date <- as.Date(ISOdate(df.yrarea$year, 7, 1))

# -----------------------------------------------------------------------------------------------------
# Annual skill metrics (remember ref and obs datasets are x and y respectively, opposite to convention)
speriods <- list(whole=c(1993,2019),
                 last=c(2007,2019),
                 weak=c(2009,2013),
                 strong=c(2014,2018))
skill <- lapply(speriods, function(x) {
  jyears <- which(df.glorys.yrarea$year %in% seq(x[1],x[2]))
  data.table::rbindlist(list(
    z500 = f_skillstats_xvec_yvec(yvar = df.glorys.yrarea$mlotst_max_500[jyears], xvar = df.opera.yrarea$mlotst_max_500[jyears]),
    z1000 = f_skillstats_xvec_yvec(yvar = df.glorys.yrarea$mlotst_max_1000[jyears], xvar = df.opera.yrarea$mlotst_max_1000[jyears]),
    z1500 = f_skillstats_xvec_yvec(yvar = df.glorys.yrarea$mlotst_max_1500[jyears], xvar = df.opera.yrarea$mlotst_max_1500[jyears]),
    z2000 = f_skillstats_xvec_yvec(yvar = df.glorys.yrarea$mlotst_max_2000[jyears], xvar = df.opera.yrarea$mlotst_max_2000[jyears])
  ), use.names = T, fill = T, idcol = "threshold")
})
names(skill) <- names(speriods)
skill <- data.table::rbindlist(skill, use.names = T, fill = T, idcol = "period")

# Add sd of x and y
skill$sdy <- skill$cvy * skill$muy
skill$sdx <- skill$sdy * skill$sdstar

# View(skill)


# -----------------------------------------------------------------------------------------------------
# Plot annual time series of MLD metrics exceeding a given Zt threshold, ALL IN ONE PLOT
pvarnames <- c("mlotst_max_500","mlotst_max_1000","mlotst_max_1500","mlotst_max_2000")
yadjust <- 1.3

png(filename = paste0("~/Desktop/Gali_2026_convectionPOC/output/Fig_S14_yrarea_",paste0(pvarnames, collapse = "_"),"_",regname4area,".png"),
    width = 16, height = 10, units = 'cm', pointsize = 8, bg = "white", res = 600, type = "quartz")

toplot <- tidyr::pivot_longer(data = df.yrarea,
                              cols = all_of(pvarnames),
                              cols_vary = "slowest",
                              names_to = "MLD_cutoff",
                              values_to = "value")
toplot$MLD_cutoff  <- plyr::mapvalues(toplot$MLD_cutoff, from = c("mlotst_max_500","mlotst_max_1000","mlotst_max_1500","mlotst_max_2000"), to = c(">500 m",">1000 m",">1500 m",">2000 m"))
toplot$MLD_cutoff <- factor(toplot$MLD_cutoff,
                  levels = c(">500 m",">1000 m",">1500 m",">2000 m"),
                  ordered = T)

p <- ggplot(toplot, aes(x = Date, y = value, group=interaction(Dataset, MLD_cutoff))) +
  geom_line(aes(colour = Dataset, linetype = MLD_cutoff)) +
  xlab("Year") +
  ylab(expression(paste("Annual deep convection area (", km^2, ")"))) +
  scale_colour_manual(values=c(GLORYS12v1="cyan3",NEMO4_ORCA1="#165CAA")) +
  scale_linetype_manual(values=c(2,1,4,3)) +
  geom_rect(data = highlight_periods, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=period),
            alpha = 0.2, inherit.aes = F) +
  scale_fill_manual(values = c("gray40", "gray80"), labels = c("Strong","Weak")) +
  annotate("text", x = toplot$Date[18], y = 2.3e6,
       label = paste0( round(100*skill$rbias[skill$period=="weak"&skill$threshold=="z500"]),"%") ) +
  annotate("text", x = toplot$Date[23], y = 2.3e6,
           label = paste0( round(100*skill$rbias[skill$period=="strong"&skill$threshold=="z500"]),"%") ) +
  annotate("text", x = toplot$Date[18], y = 0.25e6,
           label = paste0( round(100*skill$rbias[skill$period=="weak"&skill$threshold=="z1000"]),"%") ) +
  annotate("text", x = toplot$Date[23], y = 0.95e6,
           label = paste0( round(100*skill$rbias[skill$period=="strong"&skill$threshold=="z1000"]),"%") ) +
  theme_light(base_size = 12)
print(p)
dev.off()


# -----------------------------------------------------------------------------------------------------
# 
#                               Deep Convection Volume
# 
# -----------------------------------------------------------------------------------------------------

# Compute total convected volume, and convected volume between zconvmin-zconvmax (usually 500-2000 m)
zconvmin <- 500
zconvmax <- 2000
S_ztsum <- lapply(S_xy_yrmax, function(x) {
  x1 <- x[c("omlda_max","mlotst_max")]
  x2 <- lapply(x1, function(xx) {
    xx[xx>(zconvmax)] <- zconvmax
    xx <- xx-zconvmin
    xx[xx<0 & !is.na(xx)] <- 0
    return(xx)
  })
  names(x2) <- paste0(names(x1),"_500_2000")
  return(c(x1,x2))
})
S_ztsum <- lapply(S_ztsum, f_xyzt_xysum_mask, smask[["IR_LAB"]])
S_ztsum <- lapply(S_ztsum, function(x) {names(x) <- c("DCV_omlda_max","DCV","DCV_omlda_max_500_2000","DCV_500_2000"); return(x)})

df.opera.yrDCV <- cbind(
  year = operayears,
  as.data.frame(S_ztsum$REF2)
)[jselyears,]

# Load GLORYS (pre-processed) and merge with OPERA
df.glorys.yrDCV_tot <- read.csv(file = "~/Desktop/Gali_2026_convectionPOC/GLORYS12v1_code_data_MSU/output_data/deep_convection_volume_annual_IR_LAB.csv")
df.glorys.yrDCV_500_2000 <- read.csv(file = "~/Desktop/Gali_2026_convectionPOC/GLORYS12v1_code_data_MSU/output_data/deep_convection_volume_annual_IR_LAB_500_2000.csv")
df.glorys.yrDCV <- merge(x = df.glorys.yrDCV_tot, y = df.glorys.yrDCV_500_2000, by = "year")

df.yrDCV <- data.table::rbindlist(list(GLORYS12v1 = df.glorys.yrDCV, NEMO4_ORCA1 = df.opera.yrDCV), use.names = T, fill = T, idcol = "Dataset")
df.yrDCV$Date <- as.Date(ISOdate(df.yrDCV$year, 7, 1))


# -----------------------------------------------------------------------------------------------------
# Annual skill metrics (remember ref and obs datasets are x and y respectively, opposite to convention)
speriods <- list(whole=c(1993,2019),
                 last=c(2007,2019),
                 weak=c(2009,2013),
                 strong=c(2014,2018))
skill.DCV <- lapply(speriods, function(x) {
  jyears <- which(df.glorys.yrarea$year %in% seq(x[1],x[2]))
  data.table::rbindlist(list(
    DCV = f_skillstats_xvec_yvec(yvar = df.glorys.yrDCV$DCV[jyears], xvar = df.opera.yrDCV$DCV[jyears]),
    DCV_500_2000 = f_skillstats_xvec_yvec(yvar = df.glorys.yrDCV$DCV_500_2000[jyears], xvar = df.opera.yrDCV$DCV_500_2000[jyears])
  ), use.names = T, fill = T, idcol = "DCV_zrange")
})
names(skill.DCV) <- names(speriods)
skill.DCV <- data.table::rbindlist(skill.DCV, use.names = T, fill = T, idcol = "period")

# Add sd of x and y
skill.DCV$sdy <- skill.DCV$cvy * skill.DCV$muy
skill.DCV$sdx <- skill.DCV$sdy * skill.DCV$sdstar

# View(skill.DCV)


# -----------------------------------------------------------------------------------------------------
# Plot annual time series of annual (maximum) deep convection volume
pvarnames <- c("DCV_500_2000") # "DCV",

png(filename = paste0("~/Desktop/Gali_2026_convectionPOC/output/Fig_S15_yrDCV_",paste0(pvarnames, collapse = "_"),"_",regname4area,".png"),
    width = 16, height = 10, units = 'cm', pointsize = 8, bg = "white", res = 600, type = "quartz")

toplot <- tidyr::pivot_longer(data = df.yrDCV,
                              cols = all_of(pvarnames),
                              cols_vary = "slowest",
                              names_to = "DCV_zrange",
                              values_to = "value")
toplot$DCV_zrange  <- plyr::mapvalues(toplot$DCV_zrange, from = c("DCV","DCV_500_2000"), to = c("total","500-2000 m"))
toplot$DCV_zrange <- factor(toplot$DCV_zrange,
                            levels = c("total","500-2000 m"),
                            ordered = T)

p <- ggplot(toplot, aes(x = Date, y = value, group=interaction(Dataset, DCV_zrange))) +
  geom_line(aes(colour = Dataset, linetype = DCV_zrange)) +
  xlab("Year") +
  ylab(expression(paste("Annual deep convection volume (", m^3, ")"))) +
  scale_colour_manual(values=c(GLORYS12v1="cyan3",NEMO4_ORCA1="#165CAA")) +
  scale_linetype_manual(values=c(1,2)) +
  geom_rect(data = highlight_periods, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=period),
            alpha = 0.2, inherit.aes = F) +
  scale_fill_manual(values = c("gray40", "gray80"), labels = c("Strong","Weak")) +
  annotate("text", x = toplot$Date[18], y = 6e14,
           label = paste0( round(100*skill.DCV$rbias[skill.DCV$period=="weak"&skill.DCV$DCV_zrange=="DCV_500_2000"]),"%") ) +
  annotate("text", x = toplot$Date[23], y = 1.6e15,
           label = paste0( round(100*skill.DCV$rbias[skill.DCV$period=="strong"&skill.DCV$DCV_zrange=="DCV_500_2000"]),"%") ) +
  theme_light(base_size = 12)
print(p)
dev.off()



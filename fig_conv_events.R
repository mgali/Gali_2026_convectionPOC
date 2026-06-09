# Fig. 1 of convection paper and Sm figures 1, 10, 11 and 12
# Prior steps: batch processing from analyze_plot_mergedProfilesTrajSat.R (merged profile-traj-sat data and convperiod* files)
# Marti Gali Tapias, June 2024, marti.gali.tapias@gmail.com, mgali@icm.csic.es

# -----------------------------------------------------------------------------------------------------
# Input arguments
hgrid <- "orca1"
fwmoS <- c("6901480","6901486","6901523","6901524","6901527")
biospike <- ""          # Spike treatment in trajectory files. Default = "" (Tukey's criterion), Briggs-style = "b"
fex <- 6901486          # float used as example on the right panels. NUMERIC
yex <- 2015             # year used as example on the right panels. NUMERIC

# Set to TRUE or FALSE to select the display item that you want to generate
fig_events <- T           # Fig. 1 of convection paper
fig_compareconv <- T      # Fig. S1 of convection paper
fig_trajproc <- T         # Fig. S10 of convection paper
fig_trajstats <- T        # Fig. S11 of POC convection paper (requires running script twice, with biospike set to either "" or "b")
fig_ESD <- T              # Fig. S12 of convection paper

# -----------------------------------------------------------------------------------------------------
# Define paths and names of input variable files
mpath <- "~/Desktop/Gali_2026_convectionPOC/input_data/"            # merged profile, traj and satellite data with extra processing from analyze_plot_mergedProfilesTrajSat.R

# -----------------------------------------------------------------------------------------------------
# Libraries and custom-made functions
library(RNetCDF)
library(abind)
library(lubridate)
library(caTools)
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(oce)
library(ocedata)
library(fields)
library(globe) # simple Earth's sphere plot
library(multcomp)

source('~/Desktop/Gali_2026_convectionPOC/f_plothovmoller_format_data.R')
source('~/Desktop/Gali_2026_convectionPOC/f_despike_generic.R')

# Interpolate MLD data for contour plot (which cannot take non-rectangular grid: both x and y must be monotonically increasing)----
f_prep4contour1deg <- function(xin, yin, zin, dlondeg, dlatdeg) {
  iinclude <- !(zin==0 | is.na(zin) | abs(zin)==Inf)
  fld <- akima::interp(x = xin[iinclude], y = yin[iinclude], z = zin[iinclude],
                       xo = seq( min(floor(xin)), max(ceiling(xin)), dlondeg),
                       yo = seq( min(floor(yin)), max(ceiling(yin)), dlatdeg),
                       linear = T, extrap = F)
}

# -----------------------------------------------------------------------------------------------------
# Load ystats for floats with events and merge into data frame to make boxplot by periods and/or quantify delta(POC) of each event
# OPTIONS: using ystats (precomputed stats for each period, float and year) or dbdp (direct calculation)
Ystats <- lapply(fwmoS, function(fwmo) {
  
  # Uncomment for using ystats (precomputed stats for each period, float and year)
  load(paste0(mpath,"data_",fwmo,"_noTshift_9999_noreg.Rda"))                                # Load ystats dataset for each float
  return(ystats)
  
  # # Uncomment for using dbdp (direct calculation of boxplot stats from merged annual time series of all floats)
  # load(paste0(mpath,"data_",fwmo,"_noTshift_9999_noreg.Rda"))                                # Load dbdp dataset for each float
  # return(dbdp)
})
names(Ystats) <- fwmoS
Ystats <- data.table::rbindlist(Ystats, use.names = T, fill = T, idcol = "fwmo")
pnames <- c("pre","conv","post")

# -----------------------------------------------------------------------------------------------------
# Load convection events by float and merge into a single data frame
convperiods <- lapply(fwmoS, function(fwmo) {
  convperiod.m <- read.csv(file = paste0(mpath,"v0_convperiod_",fwmo,"_manu.csv")) # prepend v0_ to get previous version of convperiods_manu files (manually curated)
  as.data.frame( lapply(convperiod.m, function(x) {x <- force_tz(as.POSIXct(as.Date(x), tzone = "GMT"))} ) )
})
names(convperiods) <- fwmoS
convperiods <- data.table::rbindlist(convperiods, use.names = T, fill = F, idcol = "fwmo")

# Add variables and convert startdates to doy
convperiods <- convperiods[!is.na(convperiods$startdates)&!is.na(convperiods$enddates),]
convperiods$duration <- round(convperiods$enddates - convperiods$startdates, digits = 0)
convperiods$startdoy <- sapply(convperiods$startdates, function(x) {
  round( julian(x, origin = as.POSIXct(paste0(year(x)-1,"-12-31"))) , digits = 0)
})
convperiods$year <- lubridate::year(convperiods$enddates)

# Subset rows (events) of interest
exclude <- convperiods$fwmo=="6901485" | (convperiods$fwmo=="6901527" & convperiods$year==2015)
convperiods <- convperiods[!exclude,]

# Get index of example event (float_year)
ievent <- convperiods$fwmo==fex & year(convperiods$startdates)==yex

# Add more metrics of convection events ("chlsum" over event duration and over the whole year, etc)
toappend <- data.frame(chlcmean = NA, chl2bbp = NA, chlysum = NA, chlysum_despiked = NA, yeardaycounts = NA,
                       parmld = NA, parmld_chlcmax = NA, Kdpar_exp = NA, mld4par = NA, par0exp_chlcmax = NA,
                       GSM1_chlcmax = NA, GSM8_chlcmax = NA, OC1_chlcmax = NA, OC8_chlcmax = NA,
                       # spoc2bbp = NA, tpoc2bbp = NA,  spoc2bbp_chlcmax = NA, tpoc2bbp_chlcmax = NA,
                       # satpoc = NA, satpoc_chlcmax = NA, satpoc2bbp = NA, satpoc2bbp_chlcmax = NA, satbbp = NA, satbbp_chlcmax = NA
                       ilat = NA, ilon = NA, date = NA, days_prevmax = NA)
convperiods <- cbind(convperiods, toappend)

# Preallocate list with profiles sorted by period
Pprof <- list()

for (j in 1:dim(convperiods)[1]) {
  
  # Subset event data and load dataset for each event
  cp <- convperiods[j,]
  # print(cp$fwmo)
  load(paste0(mpath,"data_",cp$fwmo,"_noTshift_9999_noreg.Rda"))
  
  # Ensure same varnameS in all 3D arrays (exclude dissolved oxygen)
  if (j==1) {
    varnameSpprof <- grep("doxy",Mprof$varnameS,value = T, invert = T)
  }
  
  # Append to convperiods the mean (traj data) chl and chl/bbp700 during convection event: use total signal, taking into account baseline correction during "pre" period
  ys <- Ystats[Ystats$fwmo==cp$fwmo & as.numeric(Ystats$year)==year(cp$startdates) & Ystats$varname=="CHLA_ADJUSTED",]
  bs <- Ystats[Ystats$fwmo==cp$fwmo & as.numeric(Ystats$year)==year(cp$startdates) & Ystats$varname=="CHLA_ADJUSTED_over_BBP700",]
  convperiods$chlcmean[j] <- ys$mean[ys$cperiod=="conv"] - ys$mean[ys$cperiod=="pre"]
  convperiods$chl2bbp[j] <- bs$mean[ys$cperiod=="conv"]
  
  # Append to convperiods chl sum over the whole year (total signal). May be biased by negative values when using "CORRCHL" traj files with baseline correction
  dbyy <- dbdp[as.numeric(dbdp$year)==year(cp$startdates),]
  convperiods$chlysum[j] <- sum(dbyy$CHLA_ADJUSTED, na.rm=T)
  convperiods$chlysum_despiked[j] <- sum(dbyy$CHLA_ADJUSTED_despiked, na.rm=T)
  convperiods$yeardaycounts[j] <- nrow(dbyy)
  
  # Append to convperiods the DOY of CHLA maximum and the corresponding value
  tmpchl <- dbyy$CHLA_ADJUSTED_despiked
  tmpchl[dbyy$DATE < cp$startdates |  dbyy$DATE > cp$enddates] <- NA
  convperiods$chlcmax[j] <- max(tmpchl, na.rm=T)
  jchlmax <- which.max(tmpchl)
  convperiods$doy_chlcmax[j] <- julian(dbyy$DATE[jchlmax], origin = as.Date(paste0(year(cp$startdates)-1,"-12-31")))
  convperiods$bbp_chlcmax[j] <- dbyy$BBP700_despiked[jchlmax]
  
  # Subset profiles by period, excluding doxy to ensure same variables in all floats
  TMP <- Mprof$data[,,Mprof$varnameS %in%varnameSpprof] # equivalent to subsetting with !grepl("doxy",Mprof$varnameS)
  TMP0 <- TMP[,mdf$date>=(cp$startdates - as.duration("3 months")) & mdf$date<=cp$startdates &
                !is.na(mdf$tmld_SIGMAT_0.005) & mdf$tmld_SIGMAT_0.005<200,] # pre-convection period, currently not plotted
  TMP1 <- TMP[,mdf$date>=cp$startdates & mdf$date<=cp$enddates & is.na(mdf$tmld_SIGMAT_0.005),]
  TMP2 <- TMP[,mdf$date>=cp$startdates & mdf$date<=cp$enddates & !is.na(mdf$tmld_SIGMAT_0.005),]
  TMP3 <- TMP[,mdf$date>(cp$enddates + as.duration("1 month")) & mdf$date<=(cp$enddates + as.duration("3 months")) &
                !is.na(mdf$tmld_SIGMAT_0.005) & mdf$tmld_SIGMAT_0.005<40,]
  
  # Create indices for calculations on profiles
  jconv <- which(mdf$date>=cp$startdates & mdf$date<=cp$enddates) # indices of conv period
  jmaxNd <- which(mdf$date>=(dbyy$DATE[jchlmax]- as.duration("5 days")) & mdf$date<=dbyy$DATE[jchlmax]) # indices of profiles within N days (5-10) before chlmax
  jprevmax <- which.max(mdf$date[mdf$date<=dbyy$DATE[jchlmax]]) # index of previous profile
  diffdate <- mdf$date - dbyy$DATE[jchlmax]
  diffdate[diffdate<0] <- 999
  jpostmax <- which.min(diffdate) # index of posterior profile
  convperiods$days_prevmax[j] <- as.numeric(dbyy$DATE[jchlmax] - floor_date(mdf$date[jprevmax], "days")) # days
  
  # Add time and interpolate coordinates (for satellite matchups)
  convperiods$date[j] <- dbyy$DATE[jchlmax]
  convperiods$ilat[j] <- approx(x = mdf$date[c(jprevmax,jpostmax)], y = mdf$lat[c(jprevmax,jpostmax)], xout = dbyy$DATE[jchlmax])[["y"]]
  convperiods$ilon[j] <- approx(x = mdf$date[c(jprevmax,jpostmax)], y = mdf$lon[c(jprevmax,jpostmax)], xout = dbyy$DATE[jchlmax])[["y"]]
  
  # Append to convperiods the mean MLD daily PAR during the convection period and that corresponding to the maximum dCHL
  convperiods$parmld[j] <- mean(mdf$parmld[jconv], na.rm=T)
  convperiods$parmld_chlcmax[j] <- mean(mdf$parmld[jmaxNd], na.rm=T) # mean of profiles in previous N days
  # convperiods$parmld_chlcmax[j] <- mean(mdf$parmld[(jprevmax-1):jprevmax]) # mean of 2 previous profiles regardless of elapsed time
  convperiods$par0exp_chlcmax[j] <- mean(mdf$par0_exp[jmaxNd], na.rm=T)
  convperiods$mld4par[j] <- mean(mdf$tmld_SIGMAT_0.005[jconv], na.rm=T)
  convperiods$mld4par[is.na(convperiods$mld4par)] <- 1000
  convperiods$Kdpar_exp[j] <- mean(mdf$Kdpar_exp[jconv], na.rm=T)
  
  # Chl (in mg/m3) from satellite: daily and 8-day GSM and OC. PAR 8-day
  convperiods$GSM1_chlcmax[j] <- mean(mdf$day_1.CHLGSM[jmaxNd], na.rm=T)
  convperiods$GSM8_chlcmax[j] <- mean(mdf$day8_5.CHLGSM[jmaxNd], na.rm=T)
  convperiods$OC1_chlcmax[j] <- mean(mdf$day_1.CHL1[jmaxNd], na.rm=T)
  convperiods$OC8_chlcmax[j] <- mean(mdf$day8_5.CHL1[jmaxNd], na.rm=T)
  convperiods$PAR8_chlcmax[j] <- mean(mdf$day8_5.PAR[jmaxNd], na.rm=T)
  # convperiods$satchl_chlcmax[j] <- mean(mdf$day8_5.CHLGSM[jmaxNd], na.rm=T)
  
  if (j==1) {
    # create arrays at first iteration
    Pprof$pre <- TMP0       # pre-convection
    Pprof$conv.act <- TMP1  # active convection
    Pprof$conv.int <- TMP2  # interrupted convection
    Pprof$post <- TMP3      # post-convection and statified
  } else {
    # append arrays at posterior iterations
    Pprof$pre <- abind(Pprof$pre, TMP0, along = 2)
    Pprof$conv.act <- abind(Pprof$conv.act, TMP1, along = 2)
    Pprof$conv.int <- abind(Pprof$conv.int, TMP2, along = 2)
    Pprof$post <- abind(Pprof$post, TMP3, along = 2)
  }
  
}
convperiods$chlcsum <- convperiods$chlcmean * as.numeric(convperiods$duration) # daily chl x event duration (days)


# Add the bbp700/chla ratio and compute quartiles from Pprof for all periods (resulting arrays have dimension 3 x Levels x Variables)
Pprof <- lapply(Pprof, function(X) {
  abind(X, X[,,Mprof$varnameS=="chla_adjusted"] / X[,,Mprof$varnameS=="bbp700"], along = 3)
})
varnameSpprof <- c(varnameSpprof,"chla_adjusted_over_bbp700")

QPprof <- lapply(Pprof, function(X) {
  iqr <- apply(X, MARGIN = c(1,3), quantile, c(.25,.5,.75), na.rm=T)
  rng <- apply(X, MARGIN = c(1,3), range, na.rm=T)
  return( abind(iqr, rng, along = 1) )
})
names(QPprof) <- names(Pprof)

# Print IQR of top 50 m
quantile((Pprof$post[1:4,,varnameSpprof=="chla_adjusted_over_bbp700"]), c(0.05, 0.10, 0.25, 0.75), na.rm=T)


# -----------------------------------------------------------------------------------------------------
# delta(POC) during convection event for various pools
dvars <- c("CHLA_ADJUSTED","CHLA_ADJUSTED_despiked","CHLA_ADJUSTED_bdespiked",
           "BBP700","BBP700_despiked","BBP700_bdespiked")
de <- Ystats[Ystats$varname%in%dvars & Ystats$cperiod=="conv",c("fwmo","year","varname")]
de$maxdPOC <- Ystats$Max[Ystats$varname%in%dvars & Ystats$cperiod=="conv"] - Ystats$Min.[Ystats$varname%in%dvars & Ystats$cperiod=="pre"]
de$meandPOC <- Ystats$mean[Ystats$varname%in%dvars & Ystats$cperiod=="conv"] - Ystats$mean[Ystats$varname%in%dvars & Ystats$cperiod=="pre"]
de$premeanPOC <- Ystats$mean[Ystats$varname%in%dvars & Ystats$cperiod=="pre"]
de$preminPOC <- Ystats$Min.[Ystats$varname%in%dvars & Ystats$cperiod=="pre"]

# Add dmaxChl to convperiods and write out if necessary
de2merge <- data.frame(
  fwmo = de$fwmo[de$varname=="CHLA_ADJUSTED"],
  year = as.numeric(de$year[de$varname=="CHLA_ADJUSTED"]),
  dmaxChl = de$maxdPOC[de$varname=="CHLA_ADJUSTED"],
  dmaxChl_despiked = de$maxdPOC[de$varname=="CHLA_ADJUSTED_despiked"],
  dmaxChl_bdespiked = de$maxdPOC[de$varname=="CHLA_ADJUSTED_bdespiked"],
  dmeanChl = de$meandPOC[de$varname=="CHLA_ADJUSTED"],
  dmaxbbp700_x1e3 = de$maxdPOC[de$varname=="BBP700"]*1000,
  dmaxbbp700_despiked_x1e3 = de$maxdPOC[de$varname=="BBP700_despiked"]*1000,
  dmaxbbp700_bdespiked_x1e3 = de$maxdPOC[de$varname=="BBP700_bdespiked"]*1000,
  dmeanbbp700_x1e3 = de$meandPOC[de$varname=="BBP700"]*1000,
  premeanbbp700_x1e3 = de$premeanPOC[de$varname=="BBP700"]*1000,
  preminbbp700_x1e3 = de$preminPOC[de$varname=="BBP700"]*1000
)
convperiods <- merge(x = convperiods, y = de2merge, by = c("fwmo","year"))

# Write file and create symlink (always pointing to the most recent one; NOTE: old files are not overwritten)
fileconv <- paste0("~/Desktop/Gali_2026_convectionPOC/input_data/convperiods_all_",format(Sys.time(),"%Y-%m-%d"),".csv")
lastconv <- "~/Desktop/Gali_2026_convectionPOC/input_data/convperiods_all.csv"
write.csv(convperiods, file = fileconv)
system(command = paste0("ln -s -f ",fileconv," ",lastconv))


# -----------------------------------------------------------------------------------------------------
# Load daily binned and raw trajectory data for example float "fex" and year "yex" (typically 6901486 2015)
tlist <- list.files(path = "~/Desktop/Gali_2026_convectionPOC/input_data", pattern = paste0("Mtraj_",fex,"_1000m_day_stats_"), full.names = T)
tcorrchl <- grepl("CORRCHL", tlist)
load(tlist[tcorrchl])

rlist <- list.files(path = "~/Desktop/Gali_2026_convectionPOC/input_data", pattern = paste0("Mtraj_",fex,"_1000m_raw_"), full.names = T)
rcorrchl <- grepl("CORRCHL", rlist)
load(rlist[rcorrchl])

# -----------------------------------------------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -----------------------------------------------------------------------------------------------------

# Plot parameters and map region
pw <- 16
ph <- 14
rlon <- c(-62,-23)
rlat <- c(50,66)
yearS <- sort(unique(year(convperiods$startdates)))

# Z variable ranges for Hovmoller plots
zr.chl <- 10^c(-2.5,1)		# Epi. Max is 1 µg L-1 for NASTG and 10 for NASPG
zr.chls <- 10^c(-3,0)		# chl_spike
zr.bbp <- 10^c(-4,-2)     # bbp700_despiked			
zr.bbps <- 10^c(-5,-3)    # bbp700_spike
zr.bbp <- 10^c(-4.5,-2); zr.bbps <- zr.bbp
zr.spoc <- 10^c(-1,1)			# Epi Max is 10 µM for NASTG and 100 for NASPG
zr.bpoc <- 10^c(-2,1)

# Color settings
ncolors.cont <- 100
col.years <- c("#08306B", "#2171B5", "#5FC2E6", "#ABD3F4")
# col.conv <- RColorBrewer::brewer.pal(9, "Blues")[6]
names(col.years) <- yearS
pal.chl <- colorRampPalette((brewer.pal(9, "YlGn"))[1:8])
col.chl <- pal.chl(ncolors.cont)
col.chlspikes <- "burlywood3"
col.zeu <- "lightyellow"
col.tmld03 <- "darkgray"
col.tmld01 <- "darkgray"
col.tmld005 <- "black"
col.gmldchl <- "black"
col.month <- brewer.pal(12, "Paired")
col.bathy <- gray.colors(ncolors.cont)[1:75]
col.periods <-  c("#efe536","#165CAA",col.chl[70])
names(col.periods) <- c("pre","conv","post")
# Separate color for inter-convection periods
col.conv.int <- "#986AAF"


# Float symbols
pch.f <- c(15,16,18,4,17) # 24,25
cex.f <- c(2,2,2.3,2,2)
names(pch.f) <- fwmoS
names(cex.f) <- fwmoS


# gofigure
# --------------------------------------------------------------------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMBINED HOVMOELLER AND TRAJ PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------------------------------------------------------------------------------------------------------------------------

if (fig_events) { 
  
  p <- "~/Desktop/Gali_2026_convectionPOC/output/Fig_1_conv_periods.png"
  png(filename = p, width = pw, height = ph, units = 'cm', pointsize = 8, bg = "white", res = 1200, type = "quartz")
  
  # # Not nice, banded patterns appear for unknown reason
  # p <- "~/Desktop/Gali_2026_convectionPOC/output/Fig_1_conv_periods.tiff"
  # tiff(filename = p, width = pw, height = ph, units = 'cm', pointsize = 8, bg = "white", res = 300, compression = "lzw")
  
  # Multipanel setup
  # LEFT COLUMN
  m1 <- matrix(data = 1, nrow = 5, ncol = pw/2)
  m2 <- matrix(data = 2, nrow = 1, ncol = pw/2)
  m3 <- matrix(data = 3, nrow = 3, ncol = pw/2)
  m4 <- matrix(data = 4, nrow = 3, ncol = pw/2)
  m5 <- matrix(data = 5, nrow = 4, ncol = pw/2)
  # RIGHT COLUMN
  m6 <- matrix(data = 6, nrow = 6, ncol = pw/2)
  m7 <- matrix(data = 7, nrow = 3, ncol = pw/2)
  m8 <- matrix(data = 8, nrow = 7, ncol = pw/2)
  
  layout(cbind(rbind(m1, m2, m3, m4, m5), rbind(m6, m7, m8)))
  par(oma = c(1,1,0,0))
  
  
  # -----------------------------------------------------------------------------------------------------
  # ############################################ LEFT COLUMN ############################################
  # -----------------------------------------------------------------------------------------------------
  # (m1) Hovmoeller diagram of Chl
  # Load profiles from 6901486 used to illustrate effect of convection
  # NOTE: each file named data_690XXXX* has the same objects: dbdp, ystats, mdf, Mprof, mdf.orig, Mprof.orig. They will be overwritten whenever a similar file is open
  
  load(paste0(mpath,"data_",fex,"_noTshift_9999_noreg.Rda"))
  
  # # Original version time clipping (natural months)
  # selperiod <- force_tz(as.POSIXct(c(paste0(yex-1,"-12-01"),paste0(yex,"-11-30")), tzone = "GMT"))
  # dateticks1 <- seq(floor_date(min(selperiod), unit = "3 months"), ceiling_date(max(selperiod), unit = "3 months"), by = "3 months")
  # dateticks2 <- seq(floor_date(min(selperiod), unit = "months"), ceiling_date(max(selperiod), unit = "months"), by = "months")
  
  # New version time clipping (DOY)
  selperiod <- force_tz(as.POSIXct(c(paste0(yex-1,"-12-01"), paste0(yex,"-11-30")), tz = "GMT"), tzone = "GMT")
  tickorigin <- as.POSIXct(paste0(yex-1, "-12-01"), tz = "GMT")
  dateticks1 <- seq(tickorigin, max(selperiod), by = "90 days"); dateticks1 <- dateticks1[dateticks1 >= min(selperiod)]
  dateticks1[1] <- dateticks1[1] + as.duration("30 days")
  dateticks2 <- seq(tickorigin, max(selperiod), by = "30 days"); dateticks2 <- dateticks2[dateticks2 >= min(selperiod)]
  
  p2 <- f_plothovmoller_format_data(L = list(date=mdf$date, depth = Mprof$zcenter, chla_adjusted_despiked = Mprof$data[,,Mprof$varnameS=="chla_adjusted_despiked"]),
                                    xn = "date", yn = "depth", zn = "chla_adjusted_despiked",
                                    xr = selperiod, yr = c(1,1000), zr = zr.chl, zlog = T)
  p2c <- f_plothovmoller_format_data(L = list(date=mdf$date, depth = Mprof$zcenter, chla_adjusted_despiked = Mprof$data[,,Mprof$varnameS=="chla_adjusted_spike"]),
                                     xn = "date", yn = "depth", zn = "chla_adjusted_despiked",
                                     xr = selperiod, yr = c(1,1000), zr = zr.chl, zlog = F)
  
  par(mar = c(0,5,3,4))
  image(x = p2$x, y = p2$y, z = p2$z, #log = "y",
        main = "Float 6901486, 2015", cex.main = 1.5, # applies to the whole left column
        xlim = p2$xlim, ylim = rev(c(0,1000)), zlim = p2$zlim, col = col.chl, # ylim = rev(p2$ylim)
        xaxt = "n", xlab = "", yaxt = "n", ylab = "", cex.axis = 1.5, xaxs="i")
  axis.POSIXct(side = 1, at = dateticks1, labels = "", lwd = 0.5, tck=-0.03)
  axis.POSIXct(side = 1, at = dateticks2, labels = "", lwd = 0.5, tck=-0.01)
  # axis(side = 2, at = c(5,10,20,50,100,200,500,1000), lwd = 0.5, tck=-0.03, cex.axis = 1.3, las = 1) # with log y
  axis(side = 2, at = c(0,50,100,200,400,600,800,1000), mgp = c(0, 0.6, 0), lwd = 0.5, tck=-0.03, cex.axis = 1.3, las = 1) # with lin y
  mtext("Depth (m)", side = 2, line = 3, cex = 1)
  
  # Uncomment for verifying/showing that periods are coherent across variables and plots
  abline(v = convperiods$startdates[ievent], lty = 3, lwd = 0.5, col = "darkgray")
  abline(v = convperiods$enddates[ievent], lty = 3, lwd = 0.5, col = "darkgray")
  
  contour(x = p2c$x, y = p2c$y, z = p2c$z, levels = c(0.02, 0.2), lwd = c(1, 2), add = T, drawlabels = F, col = col.chlspikes)
  # contour(x = p2c$x, y = p2c$y, z = p2c$z, levels = 0.05, lwd = 2, add = T, drawlabels = F, col = col.chlspikes)
  # lines(mdf.orig$date, 4.6/mdf.orig[["Kdpar_exp"]], lwd = 2, col = col.zeu, ljoin=1)
  # lines(mdf.orig$date, mdf.orig[["gmld_CHLA_ADJUSTED"]], lty = 1, lwd = 1, col = col.gmldchl, ljoin=1)
  # lines(mdf.orig$date, mdf.orig[["tmld_SIGMAT_0.01"]], lty = 1, lwd = 1, col = col.tmld01, ljoin=1)
  lines(mdf.orig$date, mdf.orig[["tmld_SIGMAT_0.03"]], lty = 1, lwd = 1, col = col.tmld03, ljoin=1)
  lines(mdf.orig$date, mdf.orig[["tmld_SIGMAT_0.005"]], lty = 3, lwd = 1, col = col.tmld005, ljoin=0)
  abline(v = dateticks1, lty = 1, lwd = .5, col = "lightgray")
  # text(x = p2c$x[2], y = 7, labels = "A", cex = 2, bty = "n", font = 2) # with log y
  text(x = p2c$x[2], y = 60, labels = "A", cex = 2, bty = "n", font = 2) # with lin y
  
  legend("right",
         cex = 1.2, title.cex = 1.5, title = "MLD metric",
         legend = c(expression(paste("MLD"[0.005])), expression(paste("MLD"[0.03]))),
         lwd = c(1,1),
         lty = c(3,1),
         col = c(col.tmld005, col.tmld01),
         seg.len = c(1,1),
         bty = "n")
  legend(x = p2c$x[46], y = 650,
         cex = 1.2, title.cex = 1.5, title = expression(paste("FChl",italic("a")[spikes])),
         legend = c(0.02,0.2),
         col = rep(col.chlspikes,2),
         lwd = c(1, 2),
         lty = c(1,1),
         seg.len = c(1,1),
         bty = "n")
  
  # -----------------------------------------------------------------------------------------------------
  # (m2) Show periods as horizontal bars (NOTE: simpler option of plotting thick horizontal lines is not accurate)
  par(mar = c(2.6,5,0.1,4)) # previously c(0,3,3,6))
  widthvec <- as.numeric(diff(c(selperiod[1],floor_date(convperiods$startdates[ievent], "1 year"),convperiods$startdates[ievent], convperiods$enddates[ievent],selperiod[2])))
  barplot(height = rep(0.02, length(widthvec)),
          width = widthvec,
          space = 0, #add = T, offset = 0.18,
          col = c("white",col.periods), border = c("white",col.periods),
          xaxs = "i", axes = F)
  axis(side = 1, at = widthvec + 15, labels = c("","pre","conv","post"), col = rgb(1,1,1,0), cex.axis = 1.2, mgp = c(0, 0, 0))
  
  # -----------------------------------------------------------------------------------------------------
  # (m3) Trajectory data for CHLA, example float-year (see NOTE on bbp700 trajectories)
  par(mar = c(0,5,0,4))
  vname <- paste0("CHLA_ADJUSTED_",biospike,"despiked")
  
  load(paste0(mpath,"data_",fex,"_noTshift_9999_noreg.Rda")) # baseline correction: setting pre period to 0
  yref <- ystats$mean[ystats$year==yex & ystats$cperiod=="pre" & ystats$varname==vname]
  # print(yref)
  # ytot <- df.bin$mean[[paste0("CHLA_ADJUSTED_",biospike,"despiked")]] + df.bin$mean[[paste0("CHLA_ADJUSTED_",biospike,"spike")]] - yref
  ytot <- df.bin$mean$CHLA_ADJUSTED - yref
  plot(df.bin$mean$DATE, ytot, type="l", col="gray", lwd=0.1, cex.axis = 1.3, bty = "n", bg = rgb(1,1,1,1),
       ylim = c(-0.01, 0.20),
       xlim=range(dateticks2), xaxt = "n", xlab = "", yaxt = "n", ylab = "", xaxs="i")
  axis(side = 2, cex.axis=1.3, lwd = 0.5, mgp = c(0, 0.6, 0), at = seq(0,0.2,.05), labels = c(0,.05,.1,.15,.2), las = 1)
  abline(v = dateticks1, lty = 1, lwd = .5, col = "lightgray")
  points(df$DATE, df$CHLA_ADJUSTED, cex = 0.1, col=rgb(0.1,0.8,0,alpha = 0.02))
  lines(df.bin$mean$DATE, df.bin$mean[[paste0("CHLA_ADJUSTED_",biospike,"despiked")]], lwd=2, col="gray", ljoin = 2) # , lwd=2, col="gray", ljoin = 2)
  lines(df.bin$mean$DATE, ytot, lwd=0.5, col="black")
  # if (fex==6901486) abline(h = 0.007909, lty = 1, lwd = 1, col = rgb(0,.5,0,0.8))
  
  # Uncomment for verifying/showing that periods are coherent across variables and plots
  abline(v = convperiods$startdates[ievent], lty = 3, lwd = 0.5, col = "darkgray")
  abline(v = convperiods$enddates[ievent], lty = 3, lwd = 0.5, col = "darkgray")
  
  mtext(text = expression(paste(Delta,"FChl",italic("a")," (mg ",m^-3,")")), side = 2, line = 3, cex = 1)
  text(x = p2c$x[2], y = 0.18, labels = "B", cex = 2, bty = "n", font = 2)
  
  legend("topright",
         cex = 1.2, # title.cex = 1.5, title = "Periods",
         legend = c("raw data","daily","daily (despiked)"),
         col = c(rgb(0.1,0.8,0,alpha = 0.5),"black","gray"),
         lwd = c(NA,1,2),
         lty = c(NA,1,1),
         pch = c(16,NA,NA),
         seg.len = c(0,1,1),
         bty = "n")
  
  # -----------------------------------------------------------------------------------------------------
  # (m4) Trajectory data for BBP700, example float-year
  par(mar = c(0,5,0,4))
  
  # ytot <- df.bin$mean[[paste0("BBP700_",biospike,"despiked")]] + df.bin$mean[[paste0("BBP700_",biospike,"spike")]]
  ytot <- df.bin$mean$BBP700
  plot(df.bin$mean$DATE, ytot, type="l", col="gray", lwd=0.1, cex.axis = 1.3, bty = "n", bg = rgb(1,1,1,1),
       ylim = c(0.00010,0.00040),
       xlim=range(dateticks2), xaxt = "n", xlab = "", yaxt = "n", ylab = "", xaxs="i")
  axis(side = 2, cex.axis=1.3, lwd = 0.5, mgp = c(0, 0.6, 0), at = seq(0.00010,0.00040,1e-4), labels = c(.1,.2,.3,""), las = 1)
  abline(v = dateticks1, lty = 1, lwd = .5, col = "lightgray")
  points(df$DATE, df$BBP700, cex = 0.1, col=rgb(0.8,0.8,0,alpha = 0.02))
  lines(df.bin$mean$DATE, df.bin$mean[[paste0("BBP700_",biospike,"despiked")]], lwd=2, col="gray", ljoin = 2)
  lines(df.bin$mean$DATE, ytot, lwd=0.5, col="black")
  
  # Uncomment for verifying/showing that periods are coherent across variables and plots
  abline(v = convperiods$startdates[ievent], lty = 3, lwd = 0.5, col = "darkgray")
  abline(v = convperiods$enddates[ievent], lty = 3, lwd = 0.5, col = "darkgray")
  
  mtext(text = expression(paste("b"[bp700]," (",km^-1,")")), side = 2, line = 3, cex = 1)
  text(x = p2c$x[2], y = 0.00036, labels = "C", cex = 2, bty = "n", font = 2)
  
  # -----------------------------------------------------------------------------------------------------
  # (m5) Trajectory data for THETA (TEMP), example float-year
  # https://stats.oarc.ucla.edu/r/codefragments/greek_letters/
  par(mar = c(4,5,0,4))
  
  # ytot <- df.bin$mean[[paste0("CTEMP_",biospike,"despiked")]] + df.bin$mean[[paste0("CTEMP_",biospike,"spike")]]
  ytot <- df.bin$mean$CTEMP
  plot(df.bin$mean$DATE, ytot, type="l", col="gray", lwd=0.1, cex.axis = 1.3, bty = "n", bg = rgb(1,1,1,1),
       ylim=c(3.18, 3.52),
       xlim=range(dateticks2), xaxt = "n", xlab = "", yaxt = "n", ylab = "", xaxs="i", yaxs = "i")
  axis(side = 2, cex.axis=1.3, lwd = 0.5, mgp = c(0, 0.6, 0), at = seq(3.2,3.5,.1), labels = c(3.2,3.3,3.4,""), las = 1)
  abline(v = dateticks1, lty = 1, lwd = .5, col = "lightgray")
  points(df$DATE, df$CTEMP, cex = 0.1, col=rgb(0.7,0,0,alpha = 0.1))
  lines(df.bin$mean$DATE, df.bin$mean[[paste0("CTEMP_",biospike,"despiked")]], lwd=2, col="gray", ljoin = 2)
  lines(df.bin$mean$DATE, ytot, lwd=0.5, col="black")
  
  # Uncomment for verifying/showing that periods are coherent across variables and plots
  abline(v = convperiods$startdates[ievent], lty = 3, lwd = 0.5, col = "darkgray")
  abline(v = convperiods$enddates[ievent], lty = 3, lwd = 0.5, col = "darkgray")
  
  # # Common x axis for all left panels—original version
  # axis.POSIXct(side = 1, at = dateticks1, format = "%Y-%m", cex.axis = 1.3, lwd = 0.5, tck=-0.06) #selperiod
  # axis.POSIXct(side = 1, at = dateticks2, lwd = 0.5, labels = "", tck=-0.03)
  
  # Common x axis for all left panels—new version showing DOY
  axis.POSIXct(side = 1, at = dateticks1, cex.axis = 1.3, lwd = 0.5, tck = -0.06,
               labels = as.integer(format(dateticks1, "%j"))) # Major ticks (every 3 months)
  axis.POSIXct(side = 1, at = dateticks2, lwd = 0.5, tck = -0.03,
               labels = "") # Minor ticks (monthly, no labels)
  
  mtext(text = expression(paste("Temp. (",Theta, ", ºC)")), side = 2, line = 3, cex = 1)
  mtext(text = "Day of year", side = 1, line = 3, cex = 1)
  text(x = p2c$x[2], y = 3.46, labels = "D", cex = 2, bty = "n", font = 2)
  
  # -----------------------------------------------------------------------------------------------------
  # ############################################ RIGHT COLUMN ###########################################
  
  # -----------------------------------------------------------------------------------------------------
  # (m6) Profiles (with shading for IQR) of FChla/bbp700 for pre-conv-post periods, distinguishing active convection
  f_rcat <- function(v1, v2) {c(v1, rev(v2))} # concatenate vector v1 and reversed v2 to close polygon
  jvar <- which(varnameSpprof=="chla_adjusted_over_bbp700")
  jalpha <- 40
  
  par(mar = c(5,4,3,5)) # NOTE: adjusted to match the y axis of the panel on the left
  plot(QPprof$post[3,,jvar], Mprof$zcenter, type="l", lwd=0.0001, log="x", xlim = c(10,1000), ylim = c(900,0), # plot just for setting axis limits
       main = "All floats, 2014-2017", cex.main = 1.5, # applies to the whole left column
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", cex.axis = 1.3, mgp = c(0, 0.6, 0), xaxs="i", yaxs="i", bty="n")
  box(lwd=0.5)
  polygon(x = f_rcat(QPprof$post[1,,jvar],QPprof$post[3,,jvar]),
          y = f_rcat(Mprof$zcenter,Mprof$zcenter), border=FALSE, col=paste0(col.periods["post"],jalpha))
  polygon(x = f_rcat(QPprof$pre[1,,jvar],QPprof$pre[3,,jvar]),
          y = f_rcat(Mprof$zcenter,Mprof$zcenter), border=FALSE, col=paste0(col.periods["pre"],jalpha))
  polygon(x = f_rcat(QPprof$conv.int[1,,jvar],QPprof$conv.int[3,,jvar]),
          y = f_rcat(Mprof$zcenter,Mprof$zcenter), border=FALSE, col=paste0(col.conv.int,jalpha)) # rgb(0.6,0.4,.7,jalpha/100)
  polygon(x = f_rcat(QPprof$conv.act[1,,jvar],QPprof$conv.act[3,,jvar]),
          y = f_rcat(Mprof$zcenter,Mprof$zcenter), border=FALSE, col=paste0(col.periods["conv"],jalpha))
  lines(QPprof$post[2,,jvar], Mprof$zcenter, lwd=1.5, col=col.periods["post"])
  lines(QPprof$pre[2,,jvar], Mprof$zcenter, lwd=1.5, col=col.periods["pre"])
  lines(QPprof$conv.int[2,,jvar], Mprof$zcenter, lwd=1.5, col=col.conv.int)
  lines(QPprof$conv.act[2,,jvar], Mprof$zcenter, lwd=1.5, col=col.periods["conv"])
  
  # # Add the median and IQR of trajectory data for each event (requires main plot ymax of 1020)
  # for (j in 1:dim(convperiods)[1]) {
  #   cp <- convperiods[j,]
  #   ji <- Ystats$cperiod=="conv" & Ystats$fwmo==cp$fwmo & Ystats$year==cp$year
  #   points(x = Ystats$Median[ji & Ystats$varname=="CHLA_ADJUSTED_over_BBP700"],
  #          y = Ystats$Median[ji & Ystats$varname=="DEPTH"],
  #          pch = pch.f[cp$fwmo], cex = cex.f[cp$fwmo], col=adjustcolor(col.years[as.character(cp$year)], alpha.f = 0.8))
  #   lines(x = c(Ystats$`1st Qu.`[ji & Ystats$varname=="CHLA_ADJUSTED_over_BBP700"],Ystats$`3rd Qu.`[ji & Ystats$varname=="CHLA_ADJUSTED_over_BBP700"]),
  #         y = rep(Ystats$Median[ji & Ystats$varname=="DEPTH"], 2),
  #         lwd = 1, col=adjustcolor(col.years[as.character(cp$year)], alpha.f = 0.5))
  # }
  
  axis(side = 1, cex.axis=1.3, lwd = 0.5)
  mtext(text = expression(paste("FChl",italic("a"),"/","b"[bp700]," (mg ",m^-2,")")), side = 1, line = 3, cex = 1)
  axis(side = 2, at = c(0,50,100,200,400,600,800,1000), labels = rep("",8), mgp = c(0, 0.6, 0), lwd = 0.5, tck=-0.02, cex.axis = 1.3, las = 1) # with lin y
  axis(side = 4, at = c(0,50,100,200,400,600,800,1000), mgp = c(0, 0.6, 0), lwd = 0.5, tck=-0.03, cex.axis = 1.3, las = 1) # with lin y
  mtext(text = "Depth (m)", side = 4, line = 3, cex = 1, las = 3)
  text(x = 14, y = 90, labels = "E", cex = 2, bty = "n", font = 2)
  
  legend(x = 160, y = 450,
         cex = 1.2, #title.cex = 1.5, title = "Periods",
         legend = c(paste0("pre (n=",dim(Pprof$pre)[2],")"),
                    paste0("conv-act (n=",dim(Pprof$conv.act)[2],")"),
                    paste0("conv-stop (n=",dim(Pprof$conv.int)[2],")"),
                    paste0("post (n=",dim(Pprof$post)[2],")")),
         col = c(col.periods["pre"],col.periods["conv"],col.conv.int,col.periods["post"]),
         lwd = rep(2,4),
         lty = rep(1,4),
         seg.len = rep(1,4),
         bty = "n")
  
  # -----------------------------------------------------------------------------------------------------
  # (m7) Trajectory data for CHLA, all events
  # par(mar = c(1,4,0,5)) # without "boxplot" on the right
  par(mar = c(1,4,0,11)) # with boxplot on the right
  vname <- paste0("CHLA_ADJUSTED_",biospike,"despiked")
  
  for (j in 1:dim(convperiods)[1]) {
    
    cp <- convperiods[j,]                                                                         # subset event data
    load(paste0(mpath,"data_",cp$fwmo,"_noTshift_9999_noreg.Rda"))                                # Load dataset for each float and event
    dbdp <- dbdp[!is.na(dbdp$DATE) & dbdp$year==year(cp$startdates),]
    xref <- year(cp$startdates)                                                                   # Define plot variables
    xdoy <- floor( julian(dbdp$DATE, origin = as.POSIXct(paste0(xref-1,"-12-31"))) )
    yref <- ystats$mean[ystats$year==(xref) & ystats$cperiod=="pre" & ystats$varname==vname]
    ynorm <- dbdp[[vname]] - yref
    if (j==1) {
      plot(xdoy, ynorm, type="l", col=rgb(.5,.5,.5,alpha = 0.2), lwd=0.1, cex.axis = 1.3, bty = "n", bg = rgb(1,1,1,1),
           ylim = c(0, 0.20),#c(-0.01, 0.20),
           xlim=c(0,365), xaxt = "n", yaxt = "n", ylab = "", xaxs="i", yaxs="i")
      abline(h = 0)
      abline(v = dateticks1, lty = 1, lwd = .5, col = "lightgray")
    }
    lines(xdoy, ynorm, lwd=1,                                                                                 # 1000 m trace
          col=rgb(.5,.5,.5,alpha = 0.2))                                                                      # transparent gray
    tmax <- which(ynorm==max(ynorm, na.rm=T))                                                                 # time of maximum
    points(xdoy[tmax], ynorm[tmax], pch=pch.f[cp$fwmo], cex=cex.f[cp$fwmo], col=col.years[as.character(xref)])
  }
  axis(side = 1, cex.axis=1.3, lwd = 0.5, mgp = c(0, 0.6, 0), at = seq(30,330,30), las = 1) # NOTE: if plot ylim is < 0, x axis plotted at min(ylim)
  mtext("Day of year", side = 1, line = 2, cex = 1)
  axis(side = 2, cex.axis=1.3, lwd = 0.5, mgp = c(0, 0.6, 0), at = seq(0,0.2,.05), labels = c(0,.05,.1,.15,.2), las = 1)
  mtext(text = expression(paste(Delta,"FChl",italic("a")," (mg ",m^-3,")")), side = 2, line = 3, cex = 1)
  
  legend(x = 180, y = 0.18, legend = fwmoS, pch = pch.f, bty = "n", cex = 1.2, pt.cex = c(2,2,2.3,2,2), title.cex = 1.5, title = "Float") # TITLE IN SEPARATE TEXT( WITH CEX = 1.5)
  text(x = 20, y = 0.18, labels = "F", cex = 2, bty = "n", font = 2)
  
  # -----------------------------------------------------------------------------------------------------
  # (m8) Map of float trajectory during convection events overlaid on bathymetry. Contour of maximum March MLD (reanalysis?)
  
  # Load and preprocess bathymetry
  globathy <- RNetCDF::read.nc( open.nc("~/Desktop/Gali_2026_convectionPOC/input_data/bathymetry_GLORYS12v1_NASPG.nc") ) # GLORYS12v1
  bathy <- list(lonvec = globathy$longitude, latvec = globathy$latitude, zmat = -globathy$deptho)
  bathy$zmat[bathy$zmat>0] <- NA
  bathy$zmat[bathy$zmat<(-6000)] <- NA
  
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
  Amlotstmax <- Amax[,which(vname=="mlotstmax.m"),unique(dfa$year)%in%yearS]
  
  # %%%%%%%%%%%%%%%%%%%%%%%%% ALL THIS COPIED FROM map_profile_counts.NASPG.R %%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Load and preprocess grid (NOTE: See preprocess_gebco.R)
  if (hgrid=="orca1") {
    load("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/ORCA1_hgrid.Rda")
  } else if (hgrid=="lon1lat1") {
    nav_lat <- t(array(data = rep(seq(-89.5,89.5,1), 360), dim = c(180, 360)))
    nav_lon <- array(data = rep(seq(-179.5,179.5,1), 180), dim = c(360, 180))
  }
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Map
  par(mar = c(3,4,4,3))
  
  # with 1/12 degree bathy
  image(bathy$lonvec, bathy$latvec, bathy$zmat, xlab = "", ylab = "", xlim = rlon, ylim = rlat, col = col.bathy,
        mgp = c(0, 0.2, 0), tck=-0.01, cex.axis = 1, main = "")
  plot(coastline, clon = 0, clat = 0, span = c(length(bathy$lon), length(bathy$lon)),
       col = "black", bg = rgb(1,1,1,0), lwd = 0.5, axes = F, add = T)
  contour(bathy$lonvec,bathy$latvec,bathy$zmat, levels = c(-1000,-2000), lwd = 0.5, lty = 1, labcex = 1, col = col.bathy[1], drawlabels = F, add = T)
  
  # Maximum mixed layer depth (mlotstmax)
  for (jy in 1:length(yearS)) {
    
    # OPERA annual maximum mlotst (plot version 4)
    zoomm <- f_prep4contour1deg(xin = nav_lon[vcell], yin = nav_lat[vcell],
                                zin = Amlotstmax[,jy],
                                dlondeg = 0.5, dlatdeg = 0.5)
    contour(zoomm$x, zoomm$y, zoomm$z, levels = c(1000), lwd = c(1.3, 1.3, 1.2, 1.2), lty = 1, cex = 1, col = col.years[jy], drawlabels = F, add = T)
    
  }
  
  # Map of float-detected convection events
  for (j in 1:dim(convperiods)[1]) {
    
    # Subset event data and load dataset for each event
    cp <- convperiods[j,]
    load(paste0(mpath,"data_",cp$fwmo,"_noTshift_9999_noreg.Rda"))
    
    cpp <- mdf.orig[!is.na(mdf.orig$date) & mdf.orig$date>=cp$startdates & mdf.orig$date<=cp$enddates,]
    cpp <- cpp[!is.na(cpp$longitude) & !is.na(cpp$latitude),]
    if (cp$fwmo==fex & year(cp$startdates[1])==yex) {
      points(cpp$longitude[1], cpp$latitude[1], col = "orange1", pch = 0, cex = 3)
      lines(mdf$longitude[mdf$year==yex],mdf$latitude[mdf$year==yex], col = "orange1", pch = 20, lwd = 0.6)
      # print(c(cpp$longitude[1],cpp$latitude[1]))
    }
    points(cpp$longitude[1], cpp$latitude[1], col = col.years[as.character(year(cp$startdates))],
           pch = pch.f[cp$fwmo], cex = cex.f[cp$fwmo])
    if (cp$fwmo=="6901485") {
      # float with no profile data between cycles 101-125 and sensor malfunctioning. "Interpolated" coordinates
      points(-57.48, 62.47, col = "red", pch = pch.f[cp$fwmo], cex = cex.f[cp$fwmo])
    }
  }
  
  mtext("Longitude (º)", side = 1, line = 2, cex = 1)
  mtext("Latitude (º)", side = 2, line = 1, cex = 1)
  text(x = -60, y = 65, labels = "H", cex = 2, font = 2)
  
  # Legend box used as background for world map (plotted at the bottom of script)
  legend("bottomright", cex = 2.4, box.col = "black", box.lwd = 0.3, bg = "gray95",
         legend = c("    ","    "))
  
  # -----------------------------------------------------------------------------------------------------
  # Colorbars: with this method they have to be called after all plot windows have been filled
  
  # Colorbar for chlorophyll
  par( fig=c(.45,.46,.69,.95), new=TRUE, mar=c(0,0,0,0) )
  cbar.chl <- seq(p2$zlim[1],p2$zlim[2], length.out = ncolors.cont)
  cbar.tick <- c(-3,-2,-1,0,1)
  image(t(cbar.chl), cbar.chl, col = col.chl, xaxt = "n", yaxt = "n", bg = "white")
  axis(4, cex.axis=1.5, mgp = c(0, 0.6, 0), at = cbar.tick, labels = 10^cbar.tick, las = 1, lwd = 0.5)
  mtext(side = 4, text = expression(paste("FChl",italic("a")," (mg ",m^-3,")")), las = 0, line = 3)
  box(lwd = 0.5)
  
  # Colorbar for years
  # par( fig=c(.87,.88,.50,.57), new=TRUE, mar=c(0,0,0,0) ) # location without chl2bbp "boxplot" on the right
  par( fig=c(.80,.81,.50,.57), new=TRUE, mar=c(0,0,0,0) ) # location with chl2bbp "boxplot" on the right
  cbar.yr <- yearS
  cbar.tick <- cbar.yr
  image(t(cbar.yr), cbar.yr, col = rev(col.years), xaxt = "n", yaxt = "n", bg = "white")
  axis(4, cex.axis=1.2, mgp = c(0, 0.6, 0), tck = -0.03, at = cbar.tick, labels = rev(cbar.tick), las = 1, lwd = 0.5)
  mtext(side = 3, text = "Year", las = 0, line = 0.5, adj = 0.2)
  box(lwd = 0.5)
  
  # -----------------------------------------------------------------------------------------------------
  # Add small "boxplot" (manually created) for Chl/bbp700 ratio: median and IQR
  set.seed(42)  # any fixed integer
  xbox <- runif(dim(convperiods)[1], min = 0, max = 0.9)
  par( fig=c(.86,.92,.437,.63), new=TRUE, mar=c(0,0,0,0) )
  
  for (j in 1:dim(convperiods)[1]) {
    
    cp <- convperiods[j,]
    ji <- Ystats$cperiod=="conv" & Ystats$fwmo==cp$fwmo & Ystats$year==cp$year
    if (j==1) {
      plot(x = xbox[j],
           y = Ystats$Median[ji & Ystats$varname=="CHLA_ADJUSTED_over_BBP700"],
           pch = pch.f[cp$fwmo], cex = cex.f[cp$fwmo], col=adjustcolor(col.years[as.character(cp$year)], alpha.f = 1),
           # ylim = c(-68, 680),
           ylim = c(6, 1500), log = "y",
           xlim=c(-0.2,1), xaxt = "n", yaxt = "n", ylab = "", yaxs="i", bty = "n", xpd = T) # xaxs="i",
      abline(h = 10, lwd = 0.5)
    }
    points(x = xbox[j],
           y = Ystats$Median[ji & Ystats$varname=="CHLA_ADJUSTED_over_BBP700"],
           pch = pch.f[cp$fwmo], cex = cex.f[cp$fwmo], col=adjustcolor(col.years[as.character(cp$year)], alpha.f = 1))
    lines(x = rep(xbox[j],2),
          y = c(Ystats$`1st Qu.`[ji & Ystats$varname=="CHLA_ADJUSTED_over_BBP700"],Ystats$`3rd Qu.`[ji & Ystats$varname=="CHLA_ADJUSTED_over_BBP700"]),
          lwd = 1, col=adjustcolor(col.years[as.character(cp$year)], alpha.f = 0.5))
  }
  axis(side = 4, cex.axis=1.3, lwd = 0.5, mgp = c(0, 0.6, 0), at = c(10,100,1000), las = 1)
  mtext(text = expression(paste("FChl",italic("a"),"/","b"[bp700]," (mg ",m^-2,")")), side = 4, line = 4, cex = 1)
  text(x = 0, y = 700, labels = "G", cex = 2, bty = "n", font = 2)
  
  # -----------------------------------------------------------------------------------------------------
  # Overlay Earth's sphere on map (panel m8)
  # par( fig=c(.85,.95,.065,.165), new=TRUE, mar=c(0,0,0,0) ) # with par(mar(c(4,4,4,4))) for the main map
  par( fig=c(.863,.963,.05,.15), new=TRUE, mar=c(0,0,0,0) )
  globeearth(eye=c(-48.17972,56.38395), lty = 1, lwd = 1) # or globeearth(eye=place("titanic"))
  rlon <- rlon+c(-2,2)
  rlat <- rlat+c(-2,2)
  xrlon <- seq(rlon[1],rlon[2],length.out=10)
  xrlat <- seq(rlat[1],rlat[2],length.out=10)
  xmat <- cbind(c(xrlon,rep(rlon[2],10),rev(xrlon),rep(rlon[1],10)),
                c(rep(rlat[1],10),xrlat,rep(rlat[2],10),rev(xrlat)))
  globelines(loc=xmat, col="red", lwd = 2)
  
  # -----------------------------------------------------------------------------------------------------
  dev.off()
}


# --------------------------------------------------------------------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SFIG: EXAMPLE OF TRAJ PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------------------------------------------------------------------------------------------------------------------------

if (fig_trajproc) {
  rm(yref)
  rm(df)
  rm(df.bin)
  load(paste0("~/Desktop/Gali_2026_convectionPOC/input_data/Mtraj_6901486_1000m_raw_2013-06-17_2017-07-20.Rda"))
  load(paste0("~/Desktop/Gali_2026_convectionPOC/input_data/Mtraj_6901486_1000m_day_stats_2013-06-17_2017-07-20.Rda"))
  
  # Uncomment for using ystats (precomputed stats for each period, float and year). NOTE: using uncorrected baseline
  load(paste0(mpath,"data_6901486_noTshift_9999_noreg.NOCORRCHL.Rda"))
  
  # Original version time clipping (natural months)
  selperiod <- force_tz(as.POSIXct(c(paste0(yex-1,"-12-01"),paste0(yex,"-11-30")), tzone = "GMT"))
  dateticks1 <- seq(floor_date(min(selperiod), unit = "3 months"), ceiling_date(max(selperiod), unit = "3 months"), by = "3 months")
  dateticks2 <- seq(floor_date(min(selperiod), unit = "months"), ceiling_date(max(selperiod), unit = "months"), by = "months")
  
  # -----------------------------------------------------------------------------------------------------
  # Re-compute Briggs' spikes on individual measurements for illustrative purposes
  # See "preprocess_trajFile_mergeBC_binT.R"
  # SPIKES B ("Briggs"). Compute spikes with the same approach as in vertical profiles, but setting k to a value
  # that captures the time scale of interest (measurement rate of 6/hour, k = 13 ~ 2h). Only for CHLA and BBP700 variables
  st <- list(BBP700 = list(lower = 2.3e-05, upper = 8e-03 + 6.2e-05),
             CHLA = list(lower = 0.006, upper = 25),
             CDOM = list(lower = 0, upper = diff(range(df$CDOM, na.rm=T))),
             TEMP = list(lower = 0, upper = diff(range(df$TEMP_ADJUSTED, na.rm=T))))
  bspikeappend <- data.frame(DATE = df$DATE)
  for (vs in c('CHLA','CHLA_ADJUSTED','BBP700','BBP700_ADJUSTED')) {
    ifelse(vs %in% c('CHLA','CHLA_ADJUSTED'),
           tmp <- f_despike_generic(df[[vs]], k = 13, st = st$CHLA, vs),
           tmp <- f_despike_generic(df[[vs]], k = 13, st = st$BBP700, vs))
    bspikeappend <- cbind(bspikeappend, tmp)
  }
  # Define total Chla variable and baseline correction (note here we use uncorrected baseline)
  ytot <- df.bin$mean$CHLA_ADJUSTED
  yref <- ystats$mean[ystats$year==yex & ystats$cperiod=="pre" & ystats$varname=="CHLA_ADJUSTED"]
  # Subset data for detailed plots on spike treatment
  selperiod <- force_tz(as.POSIXct(c("2015-03-14 12:00:00 GMT","2015-03-18 12:00:00 GMT")), tzone = "GMT") # mid March
  dt <- as.duration("0.5 day")
  jselperioda <- df.bin$mean$DATE>selperiod[1] & df.bin$mean$DATE<selperiod[2]
  dateticksa <- seq(floor_date(min(df.bin$mean$DATE[jselperioda]), unit = "1 day"), floor_date(max(df.bin$mean$DATE[jselperioda]), unit = "1 day"), by = "1 day")
  jselperiodi <- df$DATE>selperiod[1] & df$DATE<selperiod[2]
  dateticksi <- seq(floor_date(min(df$DATE[jselperiodi]), unit = "1 day"), floor_date(max(df$DATE[jselperiodi]), unit = "1 day"), by = "1 day")
  
  
  # -----------------------------------------------------------------------------------------------------
  p <- "~/Desktop/Gali_2026_convectionPOC/output/Fig_S10_traj_processing.png"
  png(filename = p, width = 17, height = 13, units = 'cm', pointsize = 8, bg = "white", res = 600, type = "quartz")
  
  # Multipanel setup
  # LEFT COLUMN
  sm1 <- matrix(data = 1, nrow = 9, ncol = 9)
  sm2 <- matrix(data = 2, nrow = 9, ncol = 8)
  sm3 <- sm2+1
  # RIGHT COLUMN
  sm4 <- matrix(data = 4, nrow = 8, ncol = 9)
  sm5 <- matrix(data = 5, nrow = 8, ncol = 8)
  sm6 <- sm5+1
  
  layout(rbind( cbind(sm1, sm2, sm3), cbind(sm4, sm5, sm6) ))
  par(oma = c(1,1,0,0))
  
  # (sm1) Trajectory data for CHLA, example float-year, default (Tukey 's criterion) processing
  par(mar = c(1,5,4,1))
  plot(df.bin$mean$DATE, ytot, type="l", col="gray", lwd=0.1, cex.axis = 1.3, bg = rgb(1,1,1,1), bty = "n",
       main = "Tukey filter", cex.main = 1.3,
       ylim = c(0.01, 0.50), log = "y",
       xlim=range(dateticks2), xaxt = "n", xlab = "", yaxt = "n", ylab = "", xaxs="i")
  abline(v = dateticks1, lty = 1, lwd = .5, col = "lightgray")
  abline(v = convperiods$startdates[ievent], lty = 3, lwd = 0.5, col = "darkgray")
  abline(v = convperiods$enddates[ievent], lty = 3, lwd = 0.5, col = "darkgray")
  abline(v = df.bin$mean$DATE[jselperioda][2]+as.duration("12 hour"), lwd = 3, col = adjustcolor( col.periods["conv"], alpha.f = 0.3))
  points(df$DATE, df$CHLA_ADJUSTED, cex = 0.5, pch = 16, col=rgb(0.1,0.8,0,alpha = 0.02))
  lines(df.bin$mean$DATE, df.bin$mean$CHLA_ADJUSTED_despiked, lwd=2, col="gray", ljoin = 2)
  lines(df.bin$mean$DATE, ytot, lwd=0.5, col="black")
  points(x = dateticks2[1], y = yref, col = "red", pch = 1, cex = 2, lwd = 2, xpd = T)
  mtext(text = expression(paste("FChl",italic("a")," (mg ",m^-3,")")), side = 2, line = 4, cex = 1)
  axis(side = 2, cex.axis=1.3, lwd = 0.5, mgp = c(0, 0.6, 0), las = 1)
  text(x = p2c$x[2], y = 0.4, labels = "A", cex = 2, bty = "n", font = 2)
  legend("topright",
         cex = 1.2, # title.cex = 1.5, title = "Periods",
         legend = c("raw data","daily mean","daily mean (despiked)"),
         col = c(rgb(0.1,0.8,0,alpha = 0.5),"black","gray"),
         lwd = c(NA,1,2),
         lty = c(NA,1,1),
         pch = c(16,NA,NA),
         seg.len = c(0,1,1),
         bty = "n")
  
  # (sm2) Trajectory data for CHLA, example float-year, Briggs' processing
  par(mar = c(1,1,4,1))
  plot(df.bin$mean$DATE, ytot, type="l", col = "gray", lwd=0.1, cex.axis = 1.3, bg = rgb(1,1,1,1), bty = "n",
       main = "Briggs filter", cex.main = 1.3,
       ylim = c(0.01, 0.50), log = "y",
       xlim=range(dateticks2), xaxt = "n", xlab = "", yaxt = "n", ylab = "", xaxs="i")
  abline(v = dateticks1, lty = 1, lwd = .5, col="lightgray")
  abline(v = df.bin$mean$DATE[jselperioda][2]+as.duration("12 hour"), lwd = 3, col = adjustcolor( col.periods["conv"], alpha.f = 0.3))
  points(df$DATE, df$CHLA_ADJUSTED, cex = 0.5, pch = 16, col=rgb(0.1,0.8,0,alpha = 0.02))
  lines(df.bin$mean$DATE, df.bin$mean$CHLA_ADJUSTED_bdespiked, lwd=1.5, col="#7851A990", ljoin = 2)
  lines(df.bin$mean$DATE, ytot, lwd=0.5, col="black")
  axis(side = 2, cex.axis=1.3, lwd = 0.5, mgp = c(0, 0.6, 0), las = 1, at = c(.01,.02,.05,.1,.2,.5), labels = rep("",6))
  text(x = p2c$x[2], y = 0.4, labels = "B", cex = 2, bty = "n", font = 2)
  legend("topright",
         cex = 1.2, # title.cex = 1.5, title = "Periods",
         legend = c("raw data","daily mean","daily mean (despiked)"),
         col = c(rgb(0.1,0.8,0,alpha = 0.5),"black","#7851A990"),
         lwd = c(NA,1,2),
         lty = c(NA,1,1),
         pch = c(16,NA,NA),
         seg.len = c(0,1,1),
         bty = "n")
  
  # (sm3) Close up to short period in early March to illustrate processing methods
  par(mar = c(1,5,4,1))
  
  plot(df$DATE[jselperiodi], df$CHLA_ADJUSTED[jselperiodi], pch=20, cex=1.3, col=rgb(0.1,0.8,0,alpha = 0.4), bty = "n",
       main = "Comparison (close up)", cex.main = 1.3,
       ylim = c(0.05, 0.20),
       xlim = c(selperiod[1]+1*dt,selperiod[2]-0.5*dt), xaxt = "n", xlab = "", yaxt = "n", ylab = "")
  abline(v = dateticksi, lty = 1, lwd = .5, col = "lightgray")
  # Briggs treatment
  points(bspikeappend$DATE[jselperiodi],
         (bspikeappend$CHLA_ADJUSTED_bdespiked)[jselperiodi], #  + bspikeappend$CHLA_ADJUSTED_bspike
         pch = 16, cex = 0.5, col="#7851A9")
  lines(df.bin$mean$DATE[jselperioda], df.bin$mean$CHLA_ADJUSTED_bdespiked[jselperioda], col = "#7851A9", type="s", lwd=2)
  # Tukey treatment: despiked mean, mean, daily spike cutoff
  lines(df.bin$mean$DATE[jselperioda], df.bin$mean$CHLA_ADJUSTED_despiked[jselperioda], col = "gray", type="s", lwd=2)
  lines(df.bin$mean$DATE[jselperioda], df.bin$mean$CHLA_ADJUSTED[jselperioda], col = "black", type="s", lwd=1)
  lines(df.bin$qSpike$DATE[jselperioda], df.bin$qSpikeT$CHLA_ADJUSTED[jselperioda], col="blue", type="s", lwd=1, lty=2)
  points(df.bin$med$DATE[jselperioda]+dt, df.bin$med$CHLA_ADJUSTED[jselperioda], col = "black", pch = 1)
  points(df.bin$q25$DATE[jselperioda]+dt, df.bin$q25$CHLA_ADJUSTED[jselperioda], col = "black", pch = 24)
  points(df.bin$q75$DATE[jselperioda]+dt, df.bin$q75$CHLA_ADJUSTED[jselperioda], col = "black", pch = 25)
  axis(side = 2, cex.axis=1.3, lwd = 0.5, mgp = c(0, 0.6, 0), las = 1)
  text(x = df.bin$mean$DATE[jselperioda][2]-1.5*dt, y = 0.19, labels = "C", cex = 2, bty = "n", font = 2)
  legend("top",
         cex = 1, title.cex = 1.2,
         title = "Tukey", title.adj = 0,
         legend = c("cutoff","d. mean"),
         col = c("blue","gray"),
         lwd = c(1,2),
         lty = c(2,1),
         seg.len = c(1,1),
         bty = "n")
  legend("topright",
         cex = 1, title.cex = 1.2,
         title = "Briggs", title.adj = 0,
         legend = c("no spike","d. mean"),
         col = c("purple","#7851A9"),
         lwd = c(NA,1),
         lty = c(NA,1),
         pch = c(20,NA),
         seg.len = c(0,1),
         bty = "n")
  legend(x = df.bin$mean$DATE[jselperioda][3], y = 0.18,
         cex = 1, title.cex = 1.2,
         title = "Daily quartiles",
         title.adj = 0,
         legend = c("q75","median","q25"),
         col = c("black","black","black"),
         pch = c(25,1,24),
         seg.len = c(1,0,0,0),
         bty = "n")
  
  
  # (sm4) Trajectory data for BBP700, example float-year, default (Tukey's criterion) processing
  par(mar = c(4,5,0,1))
  ytot <- df.bin$mean$BBP700
  plot(df.bin$mean$DATE, ytot, type="l", col="gray", lwd=0.1, cex.axis = 1.3, bty = "n", bg = rgb(1,1,1,1),
       ylim = c(0.00010,0.00100), log = "y",
       xlim=range(dateticks2), xaxt = "n", xlab = "", yaxt = "n", ylab = "", xaxs="i")
  abline(v = dateticks1, lty = 1, lwd = .5, col = "lightgray")
  abline(v = convperiods$startdates[ievent], lty = 3, lwd = 0.5, col = "darkgray")
  abline(v = convperiods$enddates[ievent], lty = 3, lwd = 0.5, col = "darkgray")
  abline(v = df.bin$mean$DATE[jselperioda][2]+as.duration("12 hour"), lwd = 3, col = adjustcolor( col.periods["conv"], alpha.f = 0.3))
  points(df$DATE, df$BBP700, cex = 0.5, pch = 16, col=rgb(0.8,0.8,0,alpha = 0.02))
  lines(df.bin$mean$DATE, df.bin$mean$BBP700_despiked, lwd=2, col="gray", ljoin = 2)
  lines(df.bin$mean$DATE, ytot, lwd=0.5, col="black")
  mtext(text = expression(paste("b"[bp700]," (",m^-1,")")), side = 2, line = 4, cex = 1)
  axis(side = 2, cex.axis=1.3, lwd = 0.5, mgp = c(0, 0.6, 0), las = 1)
  axis.POSIXct(side = 1, at = dateticks1, format = "%Y-%m", cex.axis = 1.3, lwd = 0.5, tck=-0.03) #selperiod
  axis.POSIXct(side = 1, at = dateticks2, lwd = 0.5, labels = "", tck=-0.01)
  text(x = p2c$x[2], y = 0.00080, labels = "D", cex = 2, bty = "n", font = 2)
  legend("topright",
         cex = 1.2, # title.cex = 1.5, title = "Periods",
         legend = c("raw data","daily mean","daily mean (despiked)"),
         col = c(rgb(0.8,0.8,0,alpha = 0.5),"black","gray"),
         lwd = c(NA,1,2),
         lty = c(NA,1,1),
         pch = c(16,NA,NA),
         seg.len = c(0,1,1),
         bty = "n")
  
  # (sm5) Trajectory data for CHBBP700, example float-year, Briggs' processing
  par(mar = c(4,1,0,1))
  ytot <- df.bin$mean$BBP700
  plot(df.bin$mean$DATE, ytot, type="l", col="gray", lwd=0.1, cex.axis = 1.3, bty = "n", bg = rgb(1,1,1,1),
       ylim = c(0.00010,0.00100), log = "y",
       xlim=range(dateticks2), xaxt = "n", xlab = "", yaxt = "n", ylab = "", xaxs="i")
  abline(v = dateticks1, lty = 1, lwd = .5, col = "lightgray")
  abline(v = convperiods$startdates[ievent], lty = 3, lwd = 0.5, col = "darkgray")
  abline(v = convperiods$enddates[ievent], lty = 3, lwd = 0.5, col = "darkgray")
  abline(v = df.bin$mean$DATE[jselperioda][2]+as.duration("12 hour"), lwd = 3, col = adjustcolor( col.periods["conv"], alpha.f = 0.3))
  points(df$DATE, df$BBP700, cex = 0.5, pch = 16, col=rgb(0.8,0.8,0,alpha = 0.02))
  lines(df.bin$mean$DATE, df.bin$mean$BBP700_bdespiked, lwd=1.5, col="#7851A990", ljoin = 2)
  lines(df.bin$mean$DATE, ytot, lwd=0.5, col="black")
  axis(side = 2, cex.axis=1.3, lwd = 0.5, mgp = c(0, 0.6, 0), las = 1, at = 1e-3*c(.1,.2,.5,1), labels = rep("",4))
  axis.POSIXct(side = 1, at = dateticks1, format = "%Y-%m", cex.axis = 1.3, lwd = 0.5, tck=-0.03) #selperiod
  axis.POSIXct(side = 1, at = dateticks2, lwd = 0.5, labels = "", tck=-0.01)
  text(x = p2c$x[2], y = 0.00080, labels = "E", cex = 2, bty = "n", font = 2)
  legend("topright",
         cex = 1.2, # title.cex = 1.5, title = "Periods",
         legend = c("raw data","daily mean","daily mean (despiked)"),
         col = c(rgb(0.8,0.8,0,alpha = 0.5),"black","#7851A990"),
         lwd = c(NA,1,2),
         lty = c(NA,1,1),
         pch = c(16,NA,NA),
         seg.len = c(0,1,1),
         bty = "n")
  
  # (sm6)  Close up to short period in early March to illustrate processing methods
  par(mar = c(4,5,0,1))
  
  plot(df$DATE[jselperiodi], df$BBP700[jselperiodi], pch=20, cex=1.3, col=rgb(0.8,0.8,0,alpha = 0.4), bty = "n",
       ylim = c(0.0001, 0.0005),
       xlim = c(selperiod[1]+1*dt,selperiod[2]-0.5*dt), xaxt = "n", xlab = "", yaxt = "n", ylab = "")
  axis(side = 2, cex.axis=1.3, lwd = 0.5, mgp = c(0, 0.6, 0), las = 1)
  # Briggs treatment
  points(bspikeappend$DATE[jselperiodi],
         (bspikeappend$BBP700_bdespiked)[jselperiodi], #  + bspikeappend$BBP700_bspike
         pch = 16, cex = 0.5, col="#7851A9")
  abline(v = dateticksi, lty = 1, lwd = .5, col = "lightgray")
  lines(df.bin$mean$DATE[jselperioda], df.bin$mean$BBP700_bdespiked[jselperioda], col = "#7851A9", type="s", lwd=2)
  # Tukey treatment: despiked mean, mean, daily spike cutoff
  lines(df.bin$mean$DATE[jselperioda], df.bin$mean$BBP700_despiked[jselperioda], col = "gray", type="s", lwd=2)
  lines(df.bin$mean$DATE[jselperioda], df.bin$mean$BBP700[jselperioda], col = "black", type="s", lwd=1)
  lines(df.bin$qSpike$DATE[jselperioda], df.bin$qSpikeT$BBP700[jselperioda], col="blue", type="s", lwd=1, lty=2)
  points(df.bin$med$DATE[jselperioda]+dt, df.bin$med$BBP700[jselperioda], col = "black", pch = 1)
  points(df.bin$q25$DATE[jselperioda]+dt, df.bin$q25$BBP700[jselperioda], col = "black", pch = 24)
  points(df.bin$q75$DATE[jselperioda]+dt, df.bin$q75$BBP700[jselperioda], col = "black", pch = 25)
  axis.POSIXct(side = 1, at = dateticksi, format = "%m-%d", cex.axis = 1.3, lwd = 0.5, tck=-0.03) #selperiod
  axis.POSIXct(side = 1, at = dateticksi, lwd = 0.5, labels = "", tck=-0.01)
  text(x = df.bin$mean$DATE[jselperioda][2]-1.5*dt, y = 0.00048, labels = "F", cex = 2, bty = "n", font = 2)
  legend("top",
         cex = 1, title.cex = 1.2,
         title = "Tukey",
         title.adj = 0,
         legend = c("cutoff","d. mean"),
         col = c("blue","gray"),
         lwd = c(1,2),
         lty = c(2,1),
         seg.len = c(1,1),
         bty = "n")
  legend("topright",
         cex = 1, title.cex = 1.2,
         title = "Briggs",
         title.adj = 0,
         legend = c("no spike","d. mean"),
         col = c("purple","#7851A9"),
         lwd = c(NA,1),
         lty = c(NA,1),
         pch = c(20,NA),
         seg.len = c(0,1),
         bty = "n")
  legend(x = df.bin$mean$DATE[jselperioda][3], y = 0.00043,
         cex = 1, title.cex = 1.2,
         title = "Daily quartiles",
         title.adj = 0,
         legend = c("q75","median","q25"),
         col = c("black","black","black"),
         pch = c(25,1,24),
         seg.len = c(1,0,0,0),
         bty = "n")
  
  dev.off()
}
# stop("Stopping after SFig traj processing")

# --------------------------------------------------------------------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%% SFIG: Boxplots of spike frequency by period and spike method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------------------------------------------------------------------------------------------------------------------------
# Boxplots of spike frequency by period and spike method

if (fig_trajstats) {
  
  # Precomputed stats for each period, float and year
  Ystats$cperiod <- factor(Ystats$cperiod, ordered = T, levels = pnames)
  Ystats$fwmo_year <- paste0(Ystats$fwmo,"_",Ystats$year)
  convperiods$fwmo_year <- paste0(convperiods$fwmo,"_",convperiods$year)
  Ybox <- Ystats[Ystats$fwmo_year %in% convperiods$fwmo_year,]
  
  cstat <- "Mean"                             # either Mean or Median
  spikeproc <- paste0(biospike,"despiked")    # either "despiked" (Tukey) or "bdespiked (Briggs)
  
  # Tukey or Briggs (b-) despiking
  if (spikeproc=="despiked") {
    spikevars1 <- c("CHLA_ADJUSTED_spike","CHLA_ADJUSTED_despiked")
    spikevars2 <- c("BBP700_spike","BBP700_despiked")
  } else if (spikeproc=="bdespiked") {
    spikevars1 <- c("CHLA_ADJUSTED_bspike","CHLA_ADJUSTED_bdespiked")
    spikevars2 <- c("BBP700_bspike","BBP700_bdespiked")
  }
  
  ifelse(spikeproc == "despiked",
         p <- paste0("~/Desktop/Gali_2026_convectionPOC/output/Fig_S11A_traj_relspike_",cstat,".png"),
         p <- paste0("~/Desktop/Gali_2026_convectionPOC/output/Fig_S11B_traj_relspike_",cstat,".png"))
  png(filename = p, width = 8, height = 11, units = 'cm', pointsize = 8, bg = "white", res = 600, type = "quartz")
  
  # Multipanel setup: 2x1
  sm1 <- matrix(data = 1, nrow = 8, ncol = 8)
  sm2 <- matrix(data = 2, nrow = 8, ncol = 8)
  layout(rbind( sm1, sm2 ))
  par(oma = c(1,1,0,0))
  
  # chla spike
  jvars <- Ybox$varname %in% spikevars1
  ybox1 <- Ybox[[cstat]][jvars]
  ygrp <- Ybox$cperiod[jvars]
  yvar <- Ybox$varname[jvars]
  bshift <- median(ybox1[yvar==spikevars1[2] & ygrp=="pre"], na.rm=T)
  # print(bshift)
  ybox1 <- ybox1 - bshift # adjust pre-convection to zero (inconsequential for analysis)
  
  # Kruskal-Wallis tests
  k1d <- kruskal.test(y ~ grp, data = data.frame(y=ybox1[yvar==spikevars1[2]], grp=ygrp[yvar==spikevars1[2]]))
  k1s <- kruskal.test(y ~ grp, data = data.frame(y=ybox1[yvar==spikevars1[1]], grp=ygrp[yvar==spikevars1[1]]))
  # Wilcoxon pairwise tests
  w1d <- pairwise.wilcox.test(x = ybox1[yvar==spikevars1[2]], g = ygrp[yvar==spikevars1[2]], p.adjust.method = "holm", paired = T)
  w1s <- pairwise.wilcox.test(x = ybox1[yvar==spikevars1[1]], g = ygrp[yvar==spikevars1[1]], p.adjust.method = "holm", paired = T)
  
  par(mar = c(1,5,3,1))
  boxplot(ybox1 ~ ygrp:yvar, range = F,
          notch = T, outline = F, plot = TRUE, lwd = 0.5, xaxs="i",
          border = "black", col = col.periods, cex.lab = 2, cex.axis = 1, main = "", cex.main = 1.1,
          axes = F, xlab = "", ylab = "", xlim = c(0.5,6.5), ylim = c(-.005,.075)) # 
  axis(side = 1, cex.axis=1.5, cex.lab = 2, lwd = 0.5, mgp = c(0, 0.6, 0), tck=-0.03, at = seq(1,6), labels = rep("",6), las = 1)
  mtext(text = expression(paste(Delta,"FChl",italic("a")," (mg ",m^-3,")")), side = 2, line = 4, cex = 1)
  axis(side = 2, cex.axis=1.5, cex.lab = 2, lwd = 0.5, mgp = c(0, 0.6, 0), tck=-0.03, at = seq(-.01,.07,0.01), las = 1)
  abline(v=3.5)
  box(lwd = 0.5)
  grid(lty = 3, lwd = 0.3, col = "gray")
  mtext(side = 3, adj = .15, text = "Despiked signal", cex = 1.3)
  mtext(side = 3, adj = .8, text = "Spike signal", cex = 1.3)
  
  # bbp700 spike
  jvars2 <- Ybox$varname %in% spikevars2
  ybox2 <- Ybox[[cstat]][jvars2]
  ygrp <- Ybox$cperiod[jvars2]
  yvar <- Ybox$varname[jvars2]
  
  # Kruskal-Wallis tests
  k2d <- kruskal.test(y ~ grp, data = data.frame(y=ybox2[yvar==spikevars2[2]], grp=ygrp[yvar==spikevars2[2]]))
  k2s <- kruskal.test(y ~ grp, data = data.frame(y=ybox2[yvar==spikevars2[1]], grp=ygrp[yvar==spikevars2[1]]))
  # Wilcoxon pairwise tests
  w2d <- pairwise.wilcox.test(x = ybox2[yvar==spikevars2[2]], g = ygrp[yvar==spikevars2[2]], p.adjust.method = "holm", paired = T)
  w2s <- pairwise.wilcox.test(x = ybox2[yvar==spikevars2[1]], g = ygrp[yvar==spikevars2[1]], p.adjust.method = "holm", paired = T)
  
  par(mar = c(4,5,0,1))
  boxplot(ybox2 ~ ygrp:yvar, range = F,
          notch = T, outline = F, plot = TRUE, lwd = 0.5, xaxs="i",
          border = "black", col = col.periods, cex.lab = 2, cex.axis = 1,
          axes = F, xlab = "", ylab =  "", xlim = c(0.5,6.5), ylim = c(0,0.00027)) # 
  axis(side = 1, cex.axis=1.5, cex.lab = 2, lwd = 0.5, mgp = c(0, 0.6, 0), tck=-0.03, at = seq(1,6), labels = rep(c("pre","conv","post"),2), las = 1)
  axis(side = 2, cex.axis=1.5, cex.lab = 2, lwd = 0.5, mgp = c(0, 0.6, 0), tck=-0.03, at = seq(0,.00030,0.00005), labels = 1e3*seq(0,.00030,0.00005), las = 1) # labels = rep("", 6), 
  mtext(text = expression(paste("b"[bp700]," (",km^-1,")")), side = 2, line = 4, cex = 1)
  abline(v=3.5)
  box(lwd = 0.5)
  grid(lty = 3, lwd = 0.3, col = "gray")
  
  dev.off()
}

if (fig_ESD) {
  # --------------------------------------------------------------------------------------------------------------------------
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SFIG: Profiles (median and IQR) of ESDbbp by period %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # --------------------------------------------------------------------------------------------------------------------------
  # See panel m6 of main Figure (paper Fig. 1)
  f_rcat <- function(v1, v2) {c(v1, rev(v2))} # concatenate vector v1 and reversed v2 to close polygon
  pvar <-  "ESDbbp"
  jvar <- which(varnameSpprof==pvar)
  jalpha <- 40
  
  p <- paste0("~/Desktop/Gali_2026_convectionPOC/output/Fig_S12_profiles_",pvar,".png")
  png(filename = p, width = 8, height = 8, units = 'cm', pointsize = 8, bg = "white", res = 600, type = "quartz")
  
  xrange <- round( range(as.vector(unlist(QPprof$post[,,jvar])), na.rm=T) , digits=1 )
  xrange[1] <- floor(xrange[1]/10)*10
  xrange[2] <- ceiling(xrange[2]/10)*10
  ifelse(pvar=="ESDbbp", xrange <- c(10,80), xrange <- c(0,10))
  
  plot(QPprof$post[3,,jvar], Mprof$zcenter, type="l", lwd=0.0001, xlim = xrange, ylim = c(900,0), # plot just for setting axis limits
       main = "", cex.main = 1.5, # applies to the whole left column
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", cex.axis = 1.3, mgp = c(0, 0.6, 0), xaxs="i", yaxs="i", bty="n")
  box(lwd=0.5)
  grid(col="gray")
  polygon(x = f_rcat(QPprof$post[1,,jvar],QPprof$post[3,,jvar]),
          y = f_rcat(Mprof$zcenter,Mprof$zcenter), border=FALSE, col=paste0(col.periods["post"],jalpha))
  polygon(x = f_rcat(QPprof$pre[1,,jvar],QPprof$pre[3,,jvar]),
          y = f_rcat(Mprof$zcenter,Mprof$zcenter), border=FALSE, col=paste0(col.periods["pre"],jalpha))
  polygon(x = f_rcat(QPprof$conv.int[1,,jvar],QPprof$conv.int[3,,jvar]),
          y = f_rcat(Mprof$zcenter,Mprof$zcenter), border=FALSE, col=paste0(col.conv.int,jalpha)) # rgb(0.6,0.4,.7,jalpha/100)
  polygon(x = f_rcat(QPprof$conv.act[1,,jvar],QPprof$conv.act[3,,jvar]),
          y = f_rcat(Mprof$zcenter,Mprof$zcenter), border=FALSE, col=paste0(col.periods["conv"],jalpha))
  lines(QPprof$post[2,,jvar], Mprof$zcenter, lwd=1.5, col=col.periods["post"])
  lines(QPprof$pre[2,,jvar], Mprof$zcenter, lwd=1.5, col=col.periods["pre"])
  lines(QPprof$conv.int[2,,jvar], Mprof$zcenter, lwd=1.5, col=col.conv.int)
  lines(QPprof$conv.act[2,,jvar], Mprof$zcenter, lwd=1.5, col=col.periods["conv"])
  
  axis(side = 3, cex.axis=1, lwd = 0.5)
  xlab <- ifelse(pvar=="ESDbbp",expression(paste("ESD"[bbp700]," (µm)")), expression(paste("ESD"[Chl]," index (-)")))
  mtext(text = xlab, side = 3, line = 2, cex = 1.2)
  axis(side = 2, at = c(0,50,100,200,400,600,800,1000), mgp = c(0, 0.6, 0), lwd = 0.5, tck=-0.03, cex.axis = 1, las = 1) # with lin y
  mtext(text = "Depth (m)", side = 2, line = 3, cex = 1.2, las = 3)
  
  legend("topright",
         cex = 1,
         legend = c(paste0("pre (n=",dim(Pprof$pre)[2],")"),
                    paste0("conv-act (n=",dim(Pprof$conv.act)[2],")"),
                    paste0("conv-stop (n=",dim(Pprof$conv.int)[2],")"),
                    paste0("post (n=",dim(Pprof$post)[2],")")),
         col = c(col.periods["pre"],col.periods["conv"],col.conv.int,col.periods["post"]),
         lwd = rep(2,4),
         lty = rep(1,4),
         seg.len = rep(1,4),
         bg = "#FFFFFF95", box.col = "#FFFFFF00")
  dev.off()
}


# --------------------------------------------------------------------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SFIG: Convection extent in OPERA vs GLORYS12v1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------------------------------------------------------------------------------------------------------------------------

if (fig_compareconv) {
  # mlotst_max from GLORYS12v1
  # NOTE: This figure here because all datasets have already been loaded and processed for Fig. 1
  # Except for GLORYS data, processed by MSU, code in GLORYS12v1_code_data_MSU
  ncfileglo <- open.nc("~/Desktop/Gali_2026_convectionPOC/input_data/mlotst_max_yearly_2014_2018_GLORYS12v1.nc")
  glorys <- RNetCDF::read.nc(ncfileglo)
  
  pw <- 16 # 17 if 6x4
  ph <- 14 # 10 if 6x4
  
  col.ope <- "#165CAA" # blue3 or "#165CAA"
  col.glo <- "cyan3" # cyan3 or steelblue1
  col.a005 <- "orange2"
  col.a01 <- "indianred2" # Formerly 
  col.a03 <- "red4" # Formerly "goldenrod1"
  
  p <- "~/Desktop/Gali_2026_convectionPOC/output/Fig_S1_conv_extent_opera_glorys.png"
  png(filename = p, width = pw, height = ph, units = 'cm', pointsize = 8, bg = "white", res = 600, type = "quartz")
  
  # Multipanel setup 4x4
  # LEFT COLUMN
  m1 <- matrix(data = 1, nrow = 7, ncol = pw/2)
  m2 <- matrix(data = 2, nrow = 7, ncol = pw/2)
  # RIGHT COLUMN
  m3 <- matrix(data = 3, nrow = 7, ncol = pw/2)
  m4 <- matrix(data = 4, nrow = 7, ncol = pw/2)
  layout(rbind(cbind(m1, m2), cbind(m3, m4)))
  
  par(oma = c(1,1,0,0))
  
  # Plot 4 panels in loop (yearS 2014 to 2017)
  for (jy in 1:length(yearS)) {
    
    
    # Map
    par(mar = c(3,3,1,1))
    
    # with 1/12 degree bathy
    image(bathy$lonvec, bathy$latvec, bathy$zmat, xlab = "", ylab = "", xlim = c(-62,-23), ylim = c(49,66), col = col.bathy,
          mgp = c(0, 0.2, 0), tck=-0.01, cex.axis = 1, main = "")
    plot(coastline, clon = 0, clat = 0, span = c(length(bathy$lon), length(bathy$lon)),
         col = "black", bg = rgb(1,1,1,0), lwd = 0.5, axes = F, add = T)
    contour(bathy$lonvec,bathy$latvec,bathy$zmat, levels = c(-1000,-2000), lwd = 0.5, lty = 1, labcex = 1, col = col.bathy[1], drawlabels = F, add = T)
    
    # OPERA annual maximum mlotst
    zoomm <- f_prep4contour1deg(xin = nav_lon[vcell], yin = nav_lat[vcell],
                                zin = Amlotstmax[,jy],
                                dlondeg = 0.5, dlatdeg = 0.5)
    contour(zoomm$x, zoomm$y, zoomm$z, levels = c(1000), lwd = 2, lty = 1, cex = 1, col = col.ope, drawlabels = F, add = T) # color of converction period in Fig. 1
    
    # GLORYS annual maximum mlotst
    contour(glorys$longitude, glorys$latitude, glorys$mlotst[,,jy], levels = c(1000), lwd = 1, lty = 1, cex = 1, col = col.glo, drawlabels = F, add = T)
    
    # Map of float-detected convection events
    yearconv <- convperiods[convperiods$year==yearS[jy],]
    
    for (j in 1:dim(yearconv)[1]) {
      
      # Subset event data and load dataset for each event
      cp <- yearconv[j,]
      load(paste0(mpath,"data_",cp$fwmo,"_noTshift_9999_noreg.Rda"))
      
      cpp <- mdf.orig[!is.na(mdf.orig$date) & mdf.orig$date>=cp$startdates & mdf.orig$date<=cp$enddates,]
      cpp <- cpp[!is.na(cpp$longitude) & !is.na(cpp$latitude),]
      # Plot all profiles in convperiod with MLD0.03 > 1000
      points(cpp$longitude[is.na(cpp$tmld_SIGMAT_0.03)], cpp$latitude[is.na(cpp$tmld_SIGMAT_0.03)], col = col.a03,
             pch = pch.f[cp$fwmo], cex = 2*cex.f[cp$fwmo]/1.4)
      # Plot profiles in convperiod with MLD0.01 > 1000
      points(cpp$longitude[is.na(cpp$tmld_SIGMAT_0.01)], cpp$latitude[is.na(cpp$tmld_SIGMAT_0.01)], col = col.a01,
             pch = pch.f[cp$fwmo], cex = 1.5*cex.f[cp$fwmo]/1.4)
      # Plot profiles in convperiod with MLD0.01 > 1000
      points(cpp$longitude[is.na(cpp$tmld_SIGMAT_0.005)], cpp$latitude[is.na(cpp$tmld_SIGMAT_0.005)], col = col.a005,
             pch = pch.f[cp$fwmo], cex = cex.f[cp$fwmo]/1.4)
    }
    
    # Annotations
    legend("topleft", legend = paste0("Year ",yearS[jy]), cex = 1.4, bg = "#FFFFFF90", box.col = "#FFFFFF90")
    if (jy %in% c(1,3)) mtext("Latitude (º)", side = 2, line = 2, cex = 1)
    if (jy %in% c(3,4)) mtext("Longitude (º)", side = 1, line = 2, cex = 1)
    if (jy==1) {
      legend("bottomright",
             cex = 1.2, title.cex = 1.3,
             title = "Argo float",
             title.col = "black",
             legend = names(pch.f), 
             text.col = rep("black",length(pch.f)),
             pch = pch.f,
             pt.cex = cex.f/1.4,
             bg = "#FFFFFF97",
             box.col = "#FFFFFF95")
      legend("topright",
             cex = 1.2, title.cex = 1.3,
             title = "Argo MLD",
             title.col = "black",
             legend = c(
               expression(paste("MLD"[0.03]," > 1000")),
               expression(paste("MLD"[0.01]," > 1000")),
               expression(paste("MLD"[0.005]," > 1000"))
             ), 
             text.col = c(col.a03, col.a01, col.a005),
             bg = "#FFFFFF97",
             box.col = "#FFFFFF95")
      text(x = -38, y = 56.5, labels = "NEMO4_ORCA1", col = col.ope, adj = 0, cex = 1.5, font = 2)
      text(x = -38, y = 58.0, labels = "GLORYS12v1", col = col.glo, adj = 0, cex = 1.5, font = 2)
    } 
  }
  dev.off()
}


# -----------------------------------------------------------------------------------------------------
# AGGREGATE TRAJ FILE STATISTICS FOR VISUAL EXAMINATION AND FOR MANUSCRIPT
Ltraj <- list()
lfwmoS <- c(fwmoS) # only floats with convective events
# lfwmoS <- c(fwmoS,"6901472","6901516") # all floats shown in supporting figures, no instrumental issues
Ltraj <- lapply(lfwmoS, function(ff) {
  llist <- list.files(path = "~/Desktop/Gali_2026_convectionPOC/input_data", pattern = paste0("Mtraj_",ff,"_1000m_day_stats_"), full.names = T)
  rcorrchl <- grepl("CORRCHL", llist)
  load(llist[rcorrchl])
  # Ensure same columns prior to rbindlist
  ldf <- df.bin$mean
  jkeep <- which( !( grepl("DOXY", names(ldf)) | grepl("BBP700_ADJUSTED", names(ldf)) | grepl("NITRATE", names(ldf)) ))
  Ltraj[[ff]] <- ldf[, jkeep]
})
Ltraj <- rbindlist(Ltraj, use.names = T, idcol = "fwmo")


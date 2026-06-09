# Analyze and plot convection events in PISCES 1D (Labrador Sea simulation)
# Compare to PISCES 3D (5-day output sims, year 2015) and to BGC-Argo observations (see fig_conv_events.R)
# June 2024-March 2025, marti.gali.tapias@gmail.com, mgali@icm.csic.es

# -----------------------------------------------------------------
# Libraries
library(fields)         # for image.plot
library(lubridate)
library(RNetCDF)
library(oce)            # So far used for oceColors() palettes
library(RColorBrewer)
library(classInt) # for function classIntervals in color scale
library(areaplot)
library(data.table)
library(dplyr)
library(tidyr)

# ---------------------------------------------------------
fig_pisces1D <- T       # Fig. S13 of the POC convection paper
fig_dateChlmax <- T     # Fig. 3 of the POC convection paper

# -----------------------------------------------------------------
# Analysis settings
spinyears <- seq(1995,1997)
simyears <- seq(1998,2019)
nspinyears <- length(spinyears)
nsimyears <- length(simyears)
levmax <- 60 # vertical levels < 2000m
tpocvars <- c("sdetoc","phymisc","phydiat","zmicro","ldetoc","zmeso")
spocvars <- c("sdetoc","phymisc","phydiat","zmicro")
pktvars <- c("phymisc","phydiat","zmicro","zmeso")

# -----------------------------------------------------------------
# Plot settings

# Multiply x10 if using g instead of mol (factor of 12, rounded)

# Z variable ranges for Hovmoller plots
zr.fphy <- c(0,100)
zr.mu <- c(0.1,0.5)
zr.spoc <- 10^c(0,2)
zr.bpoc <- 10^c(-1,1)

# Colors
ncolors.cont <- 100
pal.chl <- colorRampPalette(brewer.pal(9, "YlGn")[2:8])
col.chl <- pal.chl(ncolors.cont)
col.chlspikes <- "darkgreen"
col.fphy <- "darkblue"
col.mu <- viridisLite::viridis(ncolors.cont, alpha = 1, begin = 0, direction = 1)
pal.det <- colorRampPalette((rev(brewer.pal(11, "BrBG")))[7:11])
col.det <- pal.det(ncolors.cont)
col.tmld03 <- "darkgray"
col.tmld01 <- "darkgray"
col.tmld005 <- "#F78606"
col.periods <-  c("#efe536","#165CAA",col.chl[70])
col.periods <-  c("white","white","white")

# Andrea's colors
col.phydiat <- "#1B7837"
col.phymisc <- "#5AAE61"
col.zmeso <- "#c51b7d"
col.zmicro <- "#de77ae"
col.sdetoc <- "#762A83"
col.ldetoc <- '#C2A5CF'
col.spocvars <- c(col.sdetoc, col.phymisc, col.phydiat, col.zmicro)
col.lpocvars <- c(col.ldetoc, col.zmeso)
col.tpocvars <- c(col.spocvars, col.lpocvars)
names(col.tpocvars) <- tpocvars
col.barspoc <- paste0( col.spocvars, "80") # transparency
col.barlpoc <- paste0( col.lpocvars, "50") # transparency
col.bartpoc <- c(col.barspoc, col.barlpoc)

col.areaplot <- paste0( c(col.phymisc,col.det[30],"#FFF200","#FCD12A",col.det[70]), "80") # transparency

# -----------------------------------------------------------------
# Custom functions
source('~/Desktop/Gali_2026_convectionPOC/f_plothovmoller_format_data.R')
source('~/Desktop/Gali_2026_convectionPOC/f_myNAstats.R')
source("~/Desktop/Gali_2026_convectionPOC/f_skillstats_xvec_yvec.R")

f_get_pcoords <- function() {
  u <- par("usr")
  v <- c(grconvertX(u[1:2], "user", "ndc"), grconvertY(u[3:4], "user", "ndc"))
  return(v) 
}

# -----------------------------------------------------------------
# Load PISCES1D DATA AND NEMO-DERIVED DYNA FIELDS
# -----------------------------------------------------------------

# PISCES 1D offline simulation, reference configuration: ORCA1 cell 83138 (Central Labrador Sea)
trcfile1d <- open.nc("~/Desktop/Gali_2026_convectionPOC/input_data/fig_3_S13_nemo-pisces/PISCES_1d_19950101_20191231_pisces_grid_T_3D_trc.nc")
pnames <- c("phymisc","phydiat","zmicro","zmeso","sdetoc","ldetoc")
tdiscard <- 365*nspinyears
P <- lapply(pnames, function(p) {
  # Read vars of interest (pnames), only central grid cell of 3x3, all vertical levels, discard 3 first years (spin-up)
  var.get.nc(trcfile1d, variable = p, start = c(2,2,1,tdiscard+1), count = c(1,1,levmax,9125 - tdiscard))
})
names(P) <- pnames
P$deptht <- var.get.nc(trcfile1d, variable = "deptht", start = 1, count = levmax)
P$deptht_bounds <- var.get.nc(trcfile1d, variable = "deptht_bounds", start = c(1,1), count = c(2,levmax))
close.nc(trcfile1d)

# -----------------------------------------------------------------
# NEMO-PISCES. Ooutput cropped to NASPG region. Year 2015, 5-day output for 1000 m POC tracers and MLD diagnostics
trc1000m5dREF <- open.nc("~/Desktop/Gali_2026_convectionPOC/input_data/fig_3_S13_nemo-pisces/a6qp_5d_20150101_20151231_pisces_grid_T_3D_trc_6pocvars_1000m_NASPG.nc")
REF5days2015 <- lapply(tpocvars, function(p) {var.get.nc(trc1000m5dREF, variable = p, start = c(1,1,1,1))})
names(REF5days2015) <- tpocvars
close.nc(trc1000m5dREF)
REF5days2015$phyc <- REF5days2015$phymisc + REF5days2015$phydiat

trc1000m5dRES <- open.nc("~/Desktop/Gali_2026_convectionPOC/input_data/fig_3_S13_nemo-pisces/a6qw_5d_20150101_20151231_pisces_grid_T_3D_trc_6pocvars_1000m_NASPG.nc")
RES5days2015 <- lapply(tpocvars, function(p) {var.get.nc(trc1000m5dRES, variable = p, start = c(1,1,1,1))})
names(RES5days2015) <- tpocvars
close.nc(trc1000m5dRES)
RES5days2015$phyc <- RES5days2015$phymisc + RES5days2015$phydiat

mldvars <- c("mlotst","mlotstmax","omlda","omldamax")
mld5dREF <- open.nc("~/Desktop/Gali_2026_convectionPOC/input_data/fig_3_S13_nemo-pisces/a6qp_5d_20150101_20151231_grid_T_2D_4mldvars_NASPG.nc")
MLD5days2015 <- lapply(mldvars, function(p) {var.get.nc(mld5dREF, variable = p, start = c(1,1,1,1))})
names(MLD5days2015) <- mldvars
close.nc(mld5dREF)

load(file = "~/Desktop/Gali_2026_convectionPOC/input_data/fig_3_S13_nemo-pisces/smask_NASPG_IRLAB_LABCB_bathy2000m.Rda")
smask <- lapply(smask, function(x) {x[x==0] <- NA; return(x)} )
REF5days2015 <- lapply(REF5days2015, function(x) return(x*array(smask$IR_LAB, dim = c(dim(smask$IR_LAB),73))))
RES5days2015 <- lapply(RES5days2015, function(x) return(x*array(smask$IR_LAB, dim = c(dim(smask$IR_LAB),73))))
MLD5days2015 <- lapply(MLD5days2015, function(x) return(x*array(smask$IR_LAB, dim = c(dim(smask$IR_LAB),73))))

# -----------------------------------------------------------------
# Diagnostics that are not part of detrital POC budgets, 1D
diafile1d <- open.nc("~/Desktop/Gali_2026_convectionPOC/input_data/fig_3_S13_nemo-pisces/PISCES_1d_19950101_20191231_pisces_grid_T_3D_diagnostics.nc")
mud <- var.get.nc(diafile1d, variable = "mud", start = c(2,2,1,tdiscard+1), count = c(1,1,levmax,9125 - tdiscard))
mun <- var.get.nc(diafile1d, variable = "mun", start = c(2,2,1,tdiscard+1), count = c(1,1,levmax,9125 - tdiscard))
P$muphyto <- 86400*( mud*P$phydiat + mun*P$phymisc ) / (P$phydiat + P$phymisc)
close.nc(diafile1d)

# -----------------------------------------------------------------
# Dyna fields for 1D simulation (append to the "P" list)
dyna_T_2D <- open.nc("~/Desktop/PISCES_1D/LAB_83138_1998_2019_OPERA1D_run31802099/ORCA1_OFF_PISCES_83138_LAB_ncrcat/a6u8_1d_19980101_20191231_grid_T_2D.nc")
P$mldkz5 <- var.get.nc(dyna_T_2D, variable = "mldkz5", start = c(2,2,1), count = c(1,1,8030))
P$mldr10_3 <- var.get.nc(dyna_T_2D, variable = "mldr10_3", start = c(2,2,1), count = c(1,1,8030))
close.nc(dyna_T_2D)

# -----------------------------------------------------------------
# Add date in POSIXct format
tvec_1y <- seq(1, 365)
tvec_doy <- rep(tvec_1y, nsimyears)
P$decdate <- sort(rep(simyears, 365)) + (tvec_doy - 0.5)/365 # center at noon, otherwise some dates shifted to previous day
P$date <- floor_date(date_decimal(P$decdate, tz = "GMT"), "days")

# -----------------------------------------------------------------
# Temporal binning of profile time series: daily to 5 days
# Create vector b4vec for binning data frame columns, and df b4mat for 2D binning
# Use DOY for binning (see note 4 at bottom of mergeByFloat_binZ_selected_floats.R)
ndays <- 5
ddate <- decimal_date(P$date) - year(P$date[1])
b4vec <- floor(ddate*365/ndays)*ndays
b4mat <- data.frame(DATEBIN = sort(rep(b4vec, length(P$deptht))), DEPTH = rep(P$deptht, length(P$date)))
vars2Dtobin <- c(tpocvars)

P5days <- lapply(vars2Dtobin, function(v) {
  A <- P[[v]]
  arrd <- dim(A)
  if (arrd[1]==length(P$deptht) & arrd[2]==length(P$date)) {
    TOBIN <- array(A, dim = arrd[1]*arrd[2]) # re-arrange in long format prior to binning
    TMP <- aggregate(TOBIN, by = list(DEPTH = b4mat$DEPTH, DATEBIN = b4mat$DATEBIN), function(x) nanmean(x), simplify = TRUE)
    TMP <- subset.data.frame(TMP, select = -c(DEPTH, DATEBIN))
    return( array(unlist(TMP), dim = c(arrd[1], arrd[2]/ndays) ) )
  } else {
    return(NULL)
  }
})
names(P5days) <- vars2Dtobin
P5days$date <- P$date[!(tvec_doy%%5)] - as.duration("2 days")

# -----------------------------------------------------------------
# Convert units
for ( vv in tpocvars ) {
  P[[vv]] <- P[[vv]] * 1e3 * 12 # mol/m3 to mg/m3
  P5days[[vv]] <- P5days[[vv]] * 1e3 * 12
  # REF[[vv]] <- REF[[vv]] * 1e3 * 12
  # RES[[vv]] <- RES[[vv]] * 1e3 * 12
}

# -----------------------------------------------------------------
# Calculate new variables from addition of existing ones (ensure unit conversions are previously applied to all variables)
P$phyc <- P$phymisc + P$phydiat
P$tpoc <- Reduce('+', P[tpocvars])
P$spoc <- Reduce('+', P[spocvars])
# P$c2chl <- P$phyc*12 / P$chl
P$fdiat <- 100*P$phydiat / P$phyc
P$fphy <- 100*P$phyc / P$spoc

P5days$phyc <- P5days$phymisc + P5days$phydiat
P5days$phyc <- P5days$phymisc + P5days$phydiat
P5days$tpoc <- Reduce("+", P5days[tpocvars])


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAKE FIGURES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Axis settings
dateticks1 <- seq(floor_date(min(P$date), unit = "3 months"), floor_date(max(P$date), unit = "3 months"), by = "3 months")
dateticks2 <- seq(floor_date(min(P$date), unit = "months"), floor_date(max(P$date), unit = "months"), by = "months")


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ------------------------------------------------------ FIG. S13 -------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (fig_pisces1D) {
  
  # Define period
  yy1 <- 2014; yy2 <- 2017
  
  p <- "~/Desktop/Gali_2026_convectionPOC/output/Fig_S13_pisces1D.png"
  png(filename = p, width = 15, height = 13, units = 'cm', pointsize = 8, bg = "white", res = 600, type = "quartz")
  
  # Multipanel setup
  m1 <- matrix(data = 1, nrow = 2, ncol = 15)
  m2 <- matrix(data = 2, nrow = 3, ncol = 15)
  m3 <- matrix(data = 3, nrow = 3, ncol = 15)
  m4 <- matrix(data = 4, nrow = 5, ncol = 15)
  layout(rbind(m1, m2, m3, m4))
  par(oma = c(1,1,0,0))
  
  # Time axes and related labelling
  mshift <- as.duration("1 month")
  selperiod <- force_tz(as.POSIXct(c(paste0(yy1-1,"-12-31"),paste0(yy2+1,"-01-01"))), tzone = "GMT") - mshift # Period
  jselperiod <- which(P$date > selperiod[1] & P$date < selperiod[2])
  
  if (yy1==2008 & yy2==2009) {
    xlett <- P$date[jselperiod][15]
  } else if (yy1==2014 & yy2==2017) {
    xlett <- P$date[jselperiod][30]
  }
  
  # -----------------------------------------------------------------
  # (m1) Phyto growth rates (biomass-weighted nano + diatoms)
  p11 <- f_plothovmoller_format_data(L = P,
                                     xn = "date", yn = "deptht",
                                     zn = "muphyto", zlog = F, zr = zr.mu,
                                     xr = selperiod, yr = c(8,200))
  par(mar = c(0,6,2,8))
  image(x = p11$x, y = p11$y, z = p11$z, log = "y",
        xlim = p11$xlim, ylim = rev(p11$ylim), zlim = p11$zlim, col = col.mu,
        xaxt = "n", yaxt = "n", xlab = "", ylab = "", cex.axis = 1.5)
  lines(P$date, P$mldr10_3, lty = 1, lwd = 1.5, col = col.tmld03, ljoin=1)
  lines(P$date, P$mldkz5, lty = 3, lwd = 1, col = col.tmld005, ljoin=0)
  axis.POSIXct(side = 1, at = dateticks2, labels = "", lwd = 0.5, tck=1, lty = 3, col = "gray")
  axis.POSIXct(side = 1, at = dateticks1, labels = "", lwd = 0.5, tck=-0.03)
  axis.POSIXct(side = 1, at = dateticks2, labels = "", lwd = 0.5, tck=-0.01)
  axis(side = 2, at = c(10,100,1000), lwd = 0.5, tck=-0.03, cex.axis = 1.4, las = 1)
  abline(h = c(100), col = "darkgray", lty = 3, lwd = 1)
  mtext("Depth (m)", side = 2, line = 4, cex = 1)
  text(x = xlett, y = 12, labels = "A", cex = 2, col = "white", font = 2)
  
  # -----------------------------------------------------------------
  # (m2) SPOC overlaid with contours of LPOC
  p21 <- f_plothovmoller_format_data(L = P,
                                     xn = "date", yn = "deptht", zn = "sdetoc",
                                     xr = selperiod, yr = c(8,3000), zr = zr.spoc, zlog = T)
  par(mar = c(0,6,1,8))
  image(x = p21$x, y = p21$y, z = p21$z, log = "y",
        xlim = p21$xlim, ylim = rev(p21$ylim), zlim = p21$zlim, col = col.chl,
        xaxt = "n", yaxt = "n", xlab = "", ylab = "", cex.axis = 1.5)
  # phytoplankton fraction (fdiat, %, log scale)
  p22 <- f_plothovmoller_format_data(L = P,
                                     xn = "date", yn = "deptht",
                                     zn = "fphy", zlog = F, zr = zr.fphy,
                                     xr = selperiod, yr = c(8,3000))
  contour(p22$x, p22$y, p22$z, levels = c(50,75), lwd = c(0.6,1.5), lty = rep(1,2), col = rep(col.chlspikes,2),
          drawlabels = F, add = T)
  lines(P$date, P$mldr10_3, lty = 1, lwd = 1.5, col = col.tmld03, ljoin=1)
  lines(P$date, P$mldkz5, lty = 3, lwd = 1, col = col.tmld005, ljoin=0)
  axis.POSIXct(side = 1, at = dateticks2, labels = "", lwd = 0.5, tck=1, lty = 3, col = "gray")
  axis.POSIXct(side = 1, at = dateticks1, labels = "", lwd = 0.5, tck=-0.03)
  axis.POSIXct(side = 1, at = dateticks2, labels = "", lwd = 0.5, tck=-0.01)
  axis(side = 2, at = c(10,100,1000), lwd = 0.5, tck=-0.03, cex.axis = 1.4, las = 1)
  abline(h = 1000, col = "darkgray", lty = 3, lwd = 1)
  mtext("Depth (m)", side = 2, line = 4, cex = 1)
  text(x = xlett, y = 15, labels = "B", cex = 2, col = "black", font = 2)
  
  if (yy1==2008 & yy2==2009) {
    leg1xy <- c(p21$x[580],450); leg2xy <- c(p21$x[665],450)
  } else if (yy1==2014 & yy2==2017) {
    leg1xy <- c(p21$x[1300],55); leg2xy <- c(p21$x[1320],450)
  }
  legend(x = leg1xy[1], y = leg1xy[2],
         cex = 1, title.cex = 1.2, title = "MLD metric",
         legend = c(expression(paste("MLD"[Kz])), expression(paste("MLD"[0.03]))),
         lwd = c(1,1.5),
         lty = c(3,1),
         col = c(col.tmld005, col.tmld01),
         seg.len = c(1,1),
         bg = "#FFFFFF80",
         box.col = "#FFFFFF80")
  legend(x = leg2xy[1], y = leg2xy[2],
         cex = 1, title.cex = 1.2, title = "Phyto %",
         legend = c(50,75),
         col = rep(col.chlspikes,2),
         lwd = c(1, 2),
         lty = c(1,1),
         seg.len = c(1,1),
         bg = "#FFFFFF80",
         box.col = "#FFFFFF80")
  
  # -----------------------------------------------------------------
  # (m3) ldetoc
  p3 <- f_plothovmoller_format_data(L = P,
                                    xn = "date", yn = "deptht", zn = "ldetoc",
                                    xr = selperiod, yr = c(8,3000), zr = zr.bpoc, zlog = T)
  par(mar = c(0,6,1,8))
  image(x = p3$x, y = p3$y, z = p3$z, log = "y",
        xlim = p3$xlim, ylim = rev(p3$ylim), zlim = p3$zlim, col = col.det,
        xaxt = "n", yaxt = "n", xlab = "", ylab = "", cex.axis = 1.5)
  axis.POSIXct(side = 1, at = dateticks2, labels = "", lwd = 0.5, tck=1, lty = 3, col = "gray")
  axis.POSIXct(side = 1, at = dateticks1, labels = "", lwd = 0.5, tck=-0.03)
  axis.POSIXct(side = 1, at = dateticks2, labels = "", lwd = 0.5, tck=-0.01)
  axis(side = 2, at = c(10,100,1000), lwd = 0.5, tck=-0.03, cex.axis = 1.4, las = 1)
  lines(P$date, P$mldr10_3, lty = 1, lwd = 1.5, col = col.tmld03, ljoin=1)
  lines(P$date, P$mldkz5, lty = 3, lwd = 1, col = col.tmld005, ljoin=0)
  abline(h = 1000, col = "darkgray", lty = 3, lwd = 1)
  mtext("Depth (m)", side = 2, line = 4, cex = 1)
  text(x = xlett, y = 15, labels = "C", cex = 2, col = "black", font = 2)
  
  # -----------------------------------------------------------------
  # (m4) 1000 m time series. Using areaplot package. Note need to reverse order of y variables and their colors
  # P$pkt <- P$phyc + P$zmicro + P$zmeso
  areavars <- c("phyc","sdetoc","zmicro","zmeso","ldetoc")
  
  par(mar = c(6,6,2,8))
  yarea <- as.data.frame(lapply(P[rev(areavars)], function(x) x[46,jselperiod]))
  ylines <- as.data.frame(lapply(P[rev(tpocvars)], function(x) x[46,jselperiod]))
  yleftmin <- 0
  yleftmax <- 4.5
  areaplot(P$date[jselperiod], yarea, col = rev(col.areaplot), bty = "n", lwd = 0.1, ylim = c(yleftmin,yleftmax), # 1.1*max(Reduce('+',yarea), na.rm=T)
           xaxt = "n", yaxt = "n", xlab = "Date", ylab = "", cex.axis = 1.5, cex.lab = 1.5, xaxs = "i", yaxs = "i")
  axis.POSIXct(side = 1, at = dateticks2, labels = "", lwd = 0.5, tck=1, lty = 3, col = "gray")
  axis.POSIXct(side = 1, at = dateticks1, format = "%Y-%m", cex.axis = 1.3, lwd = 0.5, tck=-0.04)
  axis.POSIXct(side = 1, at = dateticks2, labels = "", lwd = 0.5, tck=-0.01)
  lines(P$date[jselperiod], yarea$phyc, lty = 1, lwd = 2, col = col.phymisc, ljoin=1)
  lines(P$date[jselperiod], Reduce('+',ylines[spocvars]), lty = 1, lwd = 2, col = "black", ljoin=1)
  lines(P$date[jselperiod], Reduce('+',ylines), lty = 1, lwd = 1, col = "black", ljoin=1)
  
  axis(side = 2, cex.axis=1.4, lwd = 0.5, mgp = c(0, 0.6, 0), at = seq(0,5), las = 1)
  mtext(text = expression(paste("POC (mg C ",m^-3,")")), side = 2, line = 3, cex = 1)
  text(x = xlett, y = 4.1, labels = "D", cex = 2, col = "black", font = 2)
  box(lwd = 0.5)
  
  # Annotations
  # dshift <- 0 # 0 months
  dshift <- 31 # 1 month
  # Arrows to indicate deltaPOCphy and deltaPOC
  if (yy1==2008 & yy2==2009) {
    # Annotate year 2008
    arrows(x0 = P$date[jselperiod][52+dshift], x1 = P$date[jselperiod][52+dshift], y0 = 2.1, y1 = 3.9, code = 3, lwd = 0.5, length = 0.05, angle = 15)
    text(x = P$date[jselperiod][8+dshift], y = 3, labels = expression(paste(Delta,"POC")), adj = 0, cex = 1.3)
    arrows(x0 = P$date[jselperiod][62+dshift], x1 = P$date[jselperiod][62+dshift], y0 = 0.15, y1 = 1.15, code = 3, lwd = 0.5, length = 0.05, angle = 15)
    text(x = P$date[jselperiod][65+dshift], y = 0.7, labels = expression(paste(Delta,"POC"[phy])), adj = 0, cex = 1.3)
    # Arrows to indicate small and large POC
    arrows(x0 = P$date[jselperiod][213+dshift], x1 = P$date[jselperiod][213+dshift], y0 = 0, y1 = 2.3, code = 3, lwd = 0.5, length = 0.05, angle = 15)
    text(x = P$date[jselperiod][215+dshift], y = 1.5, labels = "Small POC", adj = 0, cex = 1.3)
    arrows(x0 = P$date[jselperiod][213+dshift], x1 = P$date[jselperiod][213+dshift], y0 = 2.3, y1 = 3.8, code = 3, lwd = 0.5, length = 0.05, angle = 15)
    text(x = P$date[jselperiod][215+dshift], y = 2.9, labels = "Large POC", adj = 0, cex = 1.3)
  } else  if (yy1==2014 & yy2==2017) {
    # Annotate year 2016
    arrows(x0 = P$date[jselperiod][730+15+dshift], x1 = P$date[jselperiod][730+15+dshift], y0 = 2.6, y1 = 4.2, code = 3, lwd = 0.5, length = 0.05, angle = 15)
    text(x = P$date[jselperiod][730-80+dshift], y = 3.5, labels = expression(paste(Delta,"POC")), adj = 0, cex = 1.3)
    arrows(x0 = P$date[jselperiod][730+15+dshift], x1 = P$date[jselperiod][730+15+dshift], y0 = 0.45, y1 = 1.05, code = 3, lwd = 0.5, length = 0.05, angle = 15)
    text(x = P$date[jselperiod][730-90+dshift], y = 1.3, labels = expression(paste(Delta,"POC"[phy])), adj = 0, cex = 1.3)
    # Arrows to indicate small and large POC
    arrows(x0 = P$date[jselperiod][933+dshift], x1 = P$date[jselperiod][933+dshift], y0 = 0, y1 = 1.9, code = 3, lwd = 0.5, length = 0.05, angle = 15)
    text(x = P$date[jselperiod][937+dshift], y = 1.3, labels = "Small POC", adj = 0, cex = 1.3)
    arrows(x0 = P$date[jselperiod][933+dshift], x1 = P$date[jselperiod][933+dshift], y0 = 1.9, y1 = 3.5, code = 3, lwd = 0.5, length = 0.05, angle = 15)
    text(x = P$date[jselperiod][937+dshift], y = 2.4, labels = "Large POC", adj = 0, cex = 1.3)
  }
  # -----------------------------------------------------------------------------------------------------
  # Colorbar(s): with this method they have to be called after all plot windows have been filled
  
  # Colorbar for panel m1
  par( fig=c(.90,.91,.85,.96), new=TRUE, mar=c(0,0,0,0) )
  cbar <- seq(p11$zlim[1],p11$zlim[2], length.out = ncolors.cont)
  cbar.tick <- seq(0.1,0.5,0.2)
  image(t(cbar), cbar, col = col.mu, xaxt = "n", yaxt = "n", bg = "white")
  axis(4, cex.axis=1.5, mgp = c(0, 0.6, 0), at = (cbar.tick), labels = cbar.tick, las = 1, lwd = 0.5)
  mtext(
    text = expression(paste(mu[phy] ~ " (" ~ d^-1 ~ ")")),
    side = 4, las = 0, line = 4)
  box(lwd = 0.5)
  
  # Colorbar for panel m2
  par( fig=c(.90,.91,.62,.82), new=TRUE, mar=c(0,0,0,0) )
  cbar <- seq(p21$zlim[1],p21$zlim[2], length.out = ncolors.cont)
  cbar.tick <- c(1,3,10,30,100)
  image(t(cbar), cbar, col = col.chl, xaxt = "n", yaxt = "n", bg = "white")
  axis(4, cex.axis=1.5, mgp = c(0, 0.6, 0), at = log10(cbar.tick), labels = cbar.tick, las = 1, lwd = 0.5)
  mtext(
    text = expression(paste(SPOC ~ " (mg C" ~ m^-3 ~ ")")),
    side = 4, las = 0, line = 4)
  box(lwd = 0.5)
  
  # Colorbar for panel m3
  par( fig=c(.90,.91,.39,.59), new=TRUE, mar=c(0,0,0,0) )
  cbar <- seq(p3$zlim[1],p3$zlim[2], length.out = ncolors.cont)
  cbar.tick <- c(0.1, 0.3, 1, 3, 10)
  image(t(cbar), cbar, col = col.det, xaxt = "n", yaxt = "n", bg = "white")
  axis(4, cex.axis=1.5, mgp = c(0, 0.6, 0), at = log10(cbar.tick), labels = cbar.tick, las = 1, lwd = 0.5)
  mtext(
    text = expression(paste(LPOC[det] ~ " (mg C" ~ m^-3 ~ ")")),
    side = 4, las = 0, line = 4)
  box(lwd = 0.5)
  
  # Colorbar for panel m4: Can adjust position depending on whether temperature right axis is added
  if (yy1==2008 & yy2==2009) {
    par( fig=c(.90,.91,.11,.31), new=TRUE, mar=c(0,0,0,0) )
  } else if (yy1==2014 & yy2==2017) {
    par( fig=c(.90,.91,.11,.31), new=TRUE, mar=c(0,0,0,0) )
  }
  cbar <- seq(0,4, length.out = length(areavars))
  cbar.tick <- cbar
  image(t(cbar), cbar, col = col.areaplot, xaxt = "n", yaxt = "n", bg = "white")
  axis(4, cex.axis=1.5, mgp = c(0, 0.6, 0), at = (cbar.tick), las = 1, lwd = 0.5,
       labels = c("Phy","Sdet","Micro","Meso","Ldet"))
  box(lwd = 0.5)
  
  dev.off()
  
}


# ------------------------------------------------------------------------
# Clean data for Table S3: estimates of dmaxPOCphy and dmaxPOC
# CALCULATIONS BASED ON CONVPERIODS FILE AND SATELLITE OBSERVATIONS
# ------------------------------------------------------------------------

# Load "merged data "convperiods" table: convection event metrics (generated in fig_conv_events.R)
convperiods <- read.csv("~/Desktop/Gali_2026_convectionPOC/input_data/convperiods_all.csv")

# Scale convperiods parmld to satellite 8-day data
scfact <- convperiods$PAR8_chlcmax/convperiods$par0exp_chlcmax
scfact[is.na(scfact)] <- 1
convperiods$parmld_chlcmax <- convperiods$parmld_chlcmax * scfact

# Photoacclimation model (C:Chl ratio). Jackson et al. 2017 FMARS
thetamax <- 10      # maximum Chl:C (minimum C:Chl ratio)
c2mol <- 12         # C molar mass
f_theta <- function(thetamax, I, Ik) {theta <- thetamax*(Ik/I)*(1-exp(-(I/Ik))); return(theta)}
convperiods$c2chl_prof <- 1 / f_theta(1/thetamax, convperiods$parmld_chlcmax, 3)

convperiods$dmaxPOCchl <- convperiods$dmaxChl * convperiods$c2chl_prof


# Merge data from convperiods with GlobColour satellite matchups (several variables) and OC-CCI matchups (Chl only)
dfd <- convperiods[,c("ilat","ilon","fwmo","year","startdoy","doy_chlcmax","parmld_chlcmax","par0exp_chlcmax",
                      "chlcmean","chl2bbp","dmaxChl","GSM1_chlcmax","GSM8_chlcmax","OC1_chlcmax","OC8_chlcmax")]
load("~/Desktop/Gali_2026_convectionPOC/input_data/globcolour_matchups/match_convperiods_L3bin.Rda")
dfd <- merge(x = mm, y = dfd, intersect(names(mm), names(dfd)), all = T, sort = F) %>% .[with(., order(fwmo, year)), ]
dfd$day5_1.CHLCCI <- c(0.135, NA, NA, 0.167, NA, 0.125, 0.192, 0.133, 0.130, NA, NA, 0.166) # Add OC-CCIv6 manually (data compiled in Table S3)

# Calculate C:Chl ratio for instantaneous photoacclimation case (mixed-layer PAR for MLD = 1000 m)
mld_inst <- 1000
dfd$par1000day <- dfd$day8_5.PAR * (1 / (convperiods$Kdpar_exp * mld_inst)) * (1 - exp(-convperiods$Kdpar_exp * mld_inst))
dfd$c2chl_inst <- 1 / f_theta(1/thetamax, dfd$par1000day, 1)
dfd$dmaxPOCphy_inst <- convperiods$dmaxChl * dfd$c2chl_inst

# Claculate C:Chl ratio for photoacclimation based on prior profile MLD
dfd$dmaxPOCphy_prof <- convperiods$dmaxPOCchl


dfd$prePOC <- convperiods$preminbbp700_x1e3*c2mol

# Alternative estimates of dPOC
# 0. Add delta bbp700 to be used in all BGC-Argo-based dPOC estimates. NOTES:
# a) very similar results when using absolute convperiod bbp max, rather than bbp on the day of max dCHL
# b) very similar results with full or despiked signal, except for event 11 (6901527_2014) with large spike signal
dfd$dmaxbbp <- convperiods$bbp_chlcmax - convperiods$preminbbp700_x1e3/1e3  # using bbp on the day of maximum delta chl 
# dfd$dmaxbbp_bis <- convperiods$dmaxbbp700_x1e3/1e3                          # using maximum full bbp during convperiod
# dfd$dmaxbbp_bis <- convperiods$dmaxbbp700_despiked_x1e3/1e3                 # using maximum despiked bbp during convperiod

# a) Delta bbp700 and constant POC2bbp700 = 12000
dfd$dmaxPOCbbp_a <- dfd$dmaxbbp * 1e3 * c2mol

# b) and c). Intermediate scenarios with POC2bbp profile prescribed from Galí et al. (2022) algorithm. Discarded.

# d) Delta bbp700 and a less conservative POC2bbp700: the one that ensures bbp-based POC does not exceed satellite POC
dfd$dmaxPOCbbp_d <- dfd$dmaxbbp * 40000

# Stats to estimate the POC/bbp of convection-supplied particles and consistency between BGC-Argo and satellite Chl estimates
stats <- f_skillstats_xvec_yvec(dfd$day8_5.CHL1, dfd$dmaxChl)     # GlobColour OC
stats <- f_skillstats_xvec_yvec(dfd$day8_5.CHLGSM, dfd$dmaxChl)   # GlobColour GSM
stats <- f_skillstats_xvec_yvec(dfd$day5_1.CHLCCI, dfd$dmaxChl)   # OC-CCI

# Print stats out
print(stats[c("bias","rbias","rmse","mapd")])

# ========================================================================
# Stats for main text
# ========================================================================
summary(100*dfd$dmaxPOCphy_prof/dfd$dmaxPOCbbp_d)
cor.test(x = dfd$doy_chlcmax, y = dfd$dmaxPOCphy_prof/dfd$dmaxPOCbbp_d, method = "pearson")
summary(dfd$dmaxPOCbbp_d)
summary(dfd$dmaxPOCbbp_d + dfd$prePOC)
summary(100*dfd$dmaxPOCbbp_d/dfd$prePOC)


# ------------------------------------------------------------------------
# Compare patterns of dPOCphy versus DOY in observations and simulation
# ------------------------------------------------------------------------

# Process 1000-m simulation data (REF and satellite restoring RES)
# Entire NASPG region, year 2015, 5-day output

# Create a mask (lon*lat*time) that includes only omldmax > 1000 m
# NOTE: POC jumps match date of turbocline depth (mldkz5) crossing 1000 m
mask_mldkz5_1000 <- MLD5days2015$omldamax > 1000  & !is.na(MLD5days2015$omldamax)
mask_mldkz5_1000[mask_mldkz5_1000==0] <- NA
# Index and first day mldkz5 exceeds 1000 m
i1_mldkz5_1000 <- apply(mask_mldkz5_1000, c(1,2), function(x) {
  if (sum(!is.na(x))) return(min( which(!is.na(x)) )) else return(NA)
} )
i0_mldkz5_1000 <- i1_mldkz5_1000 - 1
doy1_mldkz5_1000 <- i1_mldkz5_1000*5 - 2

# Extract the min POCphy for each dataset (regardless of MLD)
minPOCphy.REF <- apply(REF5days2015$phyc, c(1,2), function(x) {xmin <- min(x, na.rm=T); xmin[abs(xmin)==Inf] <- NA; return(xmin)} )
minPOCphy.RES <- apply(RES5days2015$phyc, c(1,2), function(x) {xmin <- min(x, na.rm=T); xmin[abs(xmin)==Inf] <- NA; return(xmin)} )

# Extract the max POCphy for each dataset during the MLD-masked period
maxPOCphy.REF <- apply(REF5days2015$phyc*mask_mldkz5_1000, c(1,2), function(x) {xmax <- max(x, na.rm=T); xmax[abs(xmax)==Inf] <- NA; return(xmax)} )
maxPOCphy.RES <- apply(RES5days2015$phyc*mask_mldkz5_1000, c(1,2), function(x) {xmax <- max(x, na.rm=T); xmax[abs(xmax)==Inf] <- NA; return(xmax)} )

# Extract the day of max POCphy for each dataset during the MLD-masked period
doymaxPOCphy.REF <- apply(REF5days2015$phyc*mask_mldkz5_1000, c(1,2), function(x) {if (sum(is.na(x))==73) return(NA) else return(which.max(x))} )
doymaxPOCphy.RES <- apply(RES5days2015$phyc*mask_mldkz5_1000, c(1,2), function(x) {if (sum(is.na(x))==73) return(NA) else return(which.max(x))} )

# Calculate dPOCphy
dPOCphy.REF <- maxPOCphy.REF - minPOCphy.REF
dPOCphy.RES <- maxPOCphy.RES - minPOCphy.RES


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------------------------------- Fig. 3: deltaChlmax vs. doy when it occurs ------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (fig_dateChlmax) {
  
  # Prepare 3D simulation data for plotting in long-format data frame
  toplot <- list(REF = data.frame(doy = as.vector(doymaxPOCphy.REF)*5-2, period = as.vector(doymaxPOCphy.REF), dPOCphy = as.vector(dPOCphy.REF)*12000),
                 SAT = data.frame(doy = as.vector(doymaxPOCphy.RES)*5-2, period = as.vector(doymaxPOCphy.RES), dPOCphy = as.vector(dPOCphy.RES)*12000))
  toplot <- data.table::rbindlist(toplot, idcol = "sim")
  toplot$sim <- factor(toplot$sim, levels = c("REF","SAT"), ordered = T)
  toplot <- toplot[!is.na(toplot$doy),]
  # set binning period (days)
  simbin <- 15 # # 10-, 15- or 30-day bins
  toplot$period <- as.character(floor(toplot$period/(simbin/5))+1)
  toplot$period <- factor(toplot$period, levels = as.character(seq(1,8)), ordered = T)
  toplot <- toplot[with(toplot, order(sim, period)), ]
  xaxis <- as.numeric(levels(toplot$period))
  xticks <- c(0,xaxis)
  # View(toplot)
  
  # Count events by 15-day period in sim data
  hist.sim <- list(
    REF = hist(as.numeric(toplot$period[toplot$sim=="REF"]), breaks = seq(0,length(xaxis))+0.5, plot = F),
    SAT = hist(as.numeric(toplot$period[toplot$sim=="SAT"]), breaks = seq(0,length(xaxis))+0.5, plot = F)
  )
  
  # Version with PAR_MLD irradiance
  zcol.max <- 4
  ncol.par <- 2*zcol.max
  parpal <- rev(brewer.pal(ncol.par, 'YlGnBu'))
  zcol <- convperiods$parmld_chlcmax # Variable used for coloring
  zcol[zcol>zcol.max] <- zcol.max
  zbreaks <- seq(0,zcol.max,0.5)
  zlabels <- zbreaks
  zlabels[seq(2,8,2)] <- ""
  
  # Common
  class <- classIntervals(zcol, ncol.par, style="fixed", fixedBreaks = zbreaks) # fixed breaks
  parcol <- findColours(class, parpal)
  
  # p <- "~/Desktop/Gali_2026_convectionPOC/output/Fig_3_dateChlmax.png"
  # png(filename = p, width = 14, height = 10, units = 'cm', pointsize = 6, bg = "white", res = 600, type = "quartz")
  
  p <- "~/Desktop/Gali_2026_convectionPOC/output/Fig_3_dateChlmax.tiff"
  tiff(filename = p, width = 14, height = 10, units = "cm", pointsize = 6, bg = "white", res = 300, compression = "lzw")
  
  # Multipanel setup
  m1 <- matrix(data = 1, nrow = 8, ncol = 10)
  m2 <- matrix(data = 2, nrow = 8, ncol = 9)
  layout(cbind(m1, m2))
  par(oma = c(1,1,0,0))
  
  # dPOCphy in BGC-Argo Observations
  # NOTE: see old code at the bottom
  par(mar = c(6,7,5,0.5))
  plot(convperiods$doy_chlcmax, convperiods$dmaxPOCchl,
       # xlim = c(25,125),
       xlim = c(26,124),
       ylim=c(0,6), pch=16, cex=6, col="black", 
       xlab = "", ylab = "", main = "Observations", cex.main = 3,
       las = 1, cex.lab = 2, cex.axis = 1.5, xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
  mtext("A", cex = 2, line = 2, adj = 0, font = 2)
  # Manual grid
  axis(side = 1, at = xticks*simbin, labels = rep("",length(xticks)), lwd = 0.5, tck=1, lty = 3, col = "gray")
  axis(side = 2, at = seq(1,5,1), labels = rep("",5), lwd = 0.5, tck=1, lty = 3, col = "gray")
  for (j in 1:dim(convperiods)[1]) {
    points(rep(convperiods$doy_chlcmax[j],2), c(10,20)*convperiods$dmaxChl[j], col="black", pch = "-", cex=5)
    lines(rep(convperiods$doy_chlcmax[j],2), c(10,20)*convperiods$dmaxChl[j], col="black", lwd = 2)
  }
  points(convperiods$doy_chlcmax, convperiods$dmaxPOCchl, pch=16, cex=5, col=parcol)
  # Manual axes
  axis(side = 1, at = xticks*simbin, labels = xticks*simbin, tck=-0.01, cex.axis = 2, mgp = c(0,1.2,0))
  axis(side = 2, at = seq(0,6,1), labels = seq(0,6,1), las = 1, tck=-0.01, cex.axis = 2)
  # Manual axis labels
  mtext("Day of year", side = 1, line = 4, cex = 2)
  mtext(expression(paste("Maximum ",Delta,"POC"[phy]," (mg ",m^-3,")")), side = 2, line = 4, cex = 2)
  
  # ----
  # dPOCphy in 3D simulations for year 2015
  par(mar = c(6,1,5,2))
  boxplot(dPOCphy ~ period, toplot, border = "white", col = "white",
          xlim = c(2.5,120/simbin+0.5), # limited x axis 30-day bins
          ylim = c(0,6),
          at = xaxis, log = "", xaxt = "n", yaxt = "n",
          xlab = "", ylab = "", main = "Simulations", cex.main = 3, adj.main = 0,
          cex.axis = 1.5, las = 1, cex.lab = 2, xaxs = "i", yaxs = "i")
  mtext("B", cex = 2, line = 2, adj = 0, font = 2)
  # Manual grid
  axis(side = 1, at = xticks+0.5, labels = rep("",length(xticks)), lwd = 0.5, tck=1, lty = 3, col = "gray")
  axis(side = 2, at = seq(1,5,1), labels = rep("",5), lwd = 0.5, tck=1, lty = 3, col = "gray")
  # Manual axes
  axis(side = 1, at = xticks+0.5, labels = xticks*simbin, lwd = 1, tck=-0.01, cex.axis = 2, mgp = c(0,1.2,0))
  # Manual axis labels
  mtext("Day of year", side = 1, line = 4, cex = 2)
  
  boxplot(dPOCphy ~ period, toplot[toplot$sim=="REF",], outline = F, border = "black", col = "gray50",
          xaxt = "n", at = xaxis - 0.1, boxwex = 0.3, yaxt = "n", add=T)
  boxplot(dPOCphy ~ period, toplot[toplot$sim=="SAT",], outline = F, border = "black", col = "orange1", # previously, "lightgreen"
          xaxt = "n", at = xaxis + 0.1, boxwex = 0.3, yaxt = "n", add=T)
  text(x = xaxis-0.2, y = rep(4.7, length(xaxis)), labels = hist.sim$REF$counts, cex = 2, col = "gray50")
  text(x = xaxis-0.2, y = rep(4.4, length(xaxis)), labels = hist.sim$SAT$counts, cex = 2, col = "orange1")
  legend("topleft",
         legend = c("REF","SAT"),
         pch = c(15,15),
         col = c("gray50","orange1"),
         cex = 3,
         pt.cex = c(3,3),
         bg = "#FFFFFF95",
         box.col = "#FFFFFF00")
  
  # Color bar for PAR
  par( fig=c(.12,.14,.55,.85), new=TRUE, mar=c(0,0,0,0) )
  cbar <- zbreaks[1:length(zbreaks)-1]
  image(t(cbar), cbar, col = parpal, xaxt = "n", yaxt = "n")
  axis(4, cex.axis=2, mgp = c(0, 1, 0), at = zbreaks-0.25, labels = zlabels, las=1)
  mtext(text = expression(paste('mol photons ',m^-2,' ',d^-1)), side = 4, cex = 1.5, line = 5) # PAR version
  box(lwd = 0.5)
  
  dev.off()
  
}


# END OF SCRIPT


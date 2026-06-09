# Compute and plot anomalies across sections in the subpolar North Atlantic
# Anomalies defined as difference between weak or strong convection periods and climatology
# Ocean sections from map_profile_counts.NASPGcoriolis.R

# Libraries
library(RNetCDF)
# library(ggplot2)
library(collapse) # install.packages("collapse")
library(tidyr) # for pivot_wider
# library(ggh4x) # for "stat_difference" in ggplot2
library(oce)
library(ocedata)
library(RColorBrewer)
library(fields)
library(cowplot)

fig_sections <- T

# -----------------------------------------------------------------------------------------------------

expid <- "a67o" # a67o, a5kq, a683, a5xl, a5gj (reference to a5* series)

mbasepath <- "~/Desktop/Gali_2026_convectionPOC/input_data/fig_4_5_S14-16_nemo-pisces/"
opath <- "~/Desktop/Gali_2026_convectionPOC/output/"
if (!dir.exists(opath)) {dir.create(opath, recursive = T)}

varnames <- c("sdetoc","ldetoc","phymisc","phydiat","zmicro","zmeso")
# varnames <- c("sdetoc","ldetoc","phymisc","phydiat","zmicro","zmeso","dissoc","dissic","thetao","so",
#               "expsdetoc","expldetoc","zdfsdetoc","zdfldetoc","zdfpkt","zdfdoc","zafsdetoc","zafldetoc","zafpkt","zafdoc",
#               "xafsdetoc","xafldetoc","xafpkt","xafdoc","yafsdetoc","yafldetoc","yafpkt","yafdoc")

varnames.phy2D <- c("omlda","omldamax","mlotst","mlotstmax")

periods <- list(
  weak = c("2009","2013"),
  strong = c("2014","2018"),
  # clim = c("1998","2019")  # a683
  clim = c("1958","2019")    # a67o
)

# -----------------------------------------------------------------------------------------------------
# Functions

# Function to load all desired variables from a given experiment (expid)
# all fluxes are defined + downwards.

f_load_3D_variables <- function(periodpath, varnames, timestat) {
  
  ncfileS <- list.files(paste0(periodpath), pattern = timestat, full.names = T)
  # print(ncfileS)
  
  # Load into list
  datalist <- lapply(varnames, function(vv) {
    
    # print(vv)
    ifelse(vv %in% c("zdfsdetoc","zdfldetoc","zdfpkt","zdfdoc"),
           vname <- paste0("/f_",vv),
           vname <- paste0("/",vv))
    # print(vname)
    filepath <- grep(pattern = vname, ncfileS, value = T)
    print(filepath)
    ncfile <- open.nc(filepath)
    return( var.get.nc(ncfile, variable = vv) )
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
    # print(fname) # for debugging
    ncfile <- open.nc(fname)
    return( var.get.nc(ncfile, variable = vv) ) # all fluxes are defined + downwards. CHECKED
    close.nc(ncfile)
  })
  names(datalist) <- varnames
  return(datalist)
}


# -----------------------------------------------------------------------------------------------------
# Load data

# Tracers and fluxes
Fpre <- list()
for (nn in names(periods)) {
  periodpath <- paste0(mbasepath,expid,"/month_",periods[[nn]][1],"_",periods[[nn]][2])
  Fpre[[nn]] <- f_load_3D_variables(periodpath, varnames, timestat = "yearmean")
}

# MLD metrics: load monthly data and calculate annual maximum
findREF <- (sort(c("a67o",expid)))[1]
expidREF <- ifelse(findREF=="a67o", "a67o", "a5gj")
Spre <- f_load_2D_variables(expidREF, mbasepath, dirname = "month", varnames.phy2D, search_pattern = "12.nc")
Spre <- lapply(Spre, function(A) {
  dima <- dim(A)
  A <- array(A, dim = c(dima[1],dima[2],dima[3]/62,dima[3]/12))
  return(apply(A, c(1,2,4), max, na.rm=T))
})

# Load lists of selected cells (inside polygon), Argo profile counts, etc: produced from "map_profile_counts.NASPGcoriolis.R"
NASPGij <- read.csv("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/NASPGij_orca1.csv")
load("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/inpolygon_NASPGcells_orca1.Rda")
load("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/section_masks_NASPGcells_orca1.Rda")

# Load ORCA1 horizontal grid, areacello
load("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/ORCA1_hgrid.Rda")
anc <- open.nc("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/areacello_Ofx_EC-Earth3_opera-control_r1i1p1f1_gn.nc")
areacello <- var.get.nc(anc, "areacello", start = c(NASPGij$istart, NASPGij$jstart), count = c(NASPGij$icounts, NASPGij$jcounts)) # ORCA1 areacello (m2!)
close.nc(anc)

# Load L75 vertical grid
deptht <- read.csv("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/zgrid_L75.csv", header = T)

# Load global bathymetry mask and subset according to NASPGcells (rectangular lon-lat domain)
maskbathy <- var.get.nc(ncfile = open.nc("~/Desktop/Gali_2026_convectionPOC/input_data/nemo_grid/a2er2000_bathy_mask.nc"),
                        variable = "bathymask",
                        start = c(NASPGij$istart, NASPGij$jstart),
                        count = c(NASPGij$icounts, NASPGij$jcounts))


# Create rectangular masks of labCB and IRLAB domains based on cell indices and bathymetry
basemask <- 0*nav_lat
f_makemask <- function(regbounds, maskcells) {
  basemask[maskcells$imatch] <- 1
  maskout <- basemask[regbounds$istart:(regbounds$istart+regbounds$icount-1),regbounds$jstart:(regbounds$jstart+regbounds$jcount-1)]
  return( maskout & maskbathy )
}
smask <- list(
  # LAB_CB = f_makemask(NASPGij, labCBcells),
  # IR_LAB = f_makemask(NASPGij, IRLABcells),
  NASPG = f_makemask(NASPGij, NASPGcells)
)

# Function to compute horizontal means for list of 4D (xyzt) or 3D (xyz) arrays, including masks for entire domain (NASPG), gyre area (IR_LAB) and central Labrador basin (LAB_CB)
# NOTE: using areacello because spatial 3D integrals will be computed from variables already multiplied by volcello
f_xyzt_xymean_mask <- function(Lxyzt, Sxy, MASK) {
  
  if (!is.null(MASK)) {
    Vxyzt <- array( (MASK * Sxy) , dim(Lxyzt[[1]]))
  }
  if (length(dim(Vxyzt))==4) {
    amargin <- c(3,4)
  } else if (length(dim(Vxyzt))==3) {
    amargin <- 3
  }
  lapply(Lxyzt, function(A) {
    AxV <- A * Vxyzt # element by element product
    AxVsum <- apply(AxV, MARGIN = amargin, sum, na.rm=T)
    Vsum <- apply(Vxyzt, MARGIN = amargin, sum, na.rm=T)
    return(AxVsum / Vsum)
  })
}

# -----------------------------------------------------------------------------------------------------
# Function to apply horizontal mask to subset a section, or regional data, without further averaging
# NOTE: can use variables already multiplied by volcello to obtain absolute mass and transports 
f_xyzt_smask <- function(Lxyzt, MASK) { # # Crop volcello array to match experiment duration (restoring experiments start on 19980101)
  
  if (!is.null(MASK)) {
    Vxyzt <- array(MASK$zmat, dim = dim(Lxyzt[[1]])) # currently not multiplying by volcello
  }
  lapply(Lxyzt, function(A) {
    SECT <- array(A[which(Vxyzt==1)], dim = c( sum(MASK$zmat==1, na.rm=T), dim(A)[length(dim(A))] ))
    SECT[SECT==0 & !is.na(SECT)] <- NA
    return(SECT)
  })
}

# -----------------------------------------------------------------------------------------------------
# LUMP TOGETHER ALL POC TRACERS AND CORRESPONDING FLUXES. UPDATE VARNAMES

Fpre <- lapply(Fpre, function(X) {
  OUT <- list()
  OUT$poc <- Reduce('+', X[c("sdetoc","ldetoc","phymisc","phydiat","zmicro","zmeso")])
  # OUT$zdfpoc <- Reduce('+', X[c("zdfsdetoc","zdfldetoc","zdfpkt")])
  # OUT$zafpoc <- Reduce('+', X[c("zafsdetoc","zafldetoc","zafpkt")])
  # OUT$xafpoc <- Reduce('+', X[c("xafsdetoc","xafldetoc","xafpkt")])
  # OUT$yafpoc <- Reduce('+', X[c("yafsdetoc","yafldetoc","yafpkt")])
  # OUT$exppoc <- Reduce('+', X[c("expsdetoc","expldetoc")])
  # if (expid=="a67o") {
  #   return( c(OUT,X[c("dissoc","dissic","thetao","so","zdfdoc","zafdoc","xafdoc","yafdoc","uo","vo","wo")]) )
  # } else {
  #   return( c(OUT,X[c("dissoc","dissic","thetao","so","zdfdoc","zafdoc","xafdoc","yafdoc")]) )
  # }
  OUT
})
varnames <- names(Fpre[[1]])

# -----------------------------------------------------------------------------------------------------
# Calculate MLDmax for each period
# This step can be omitted if code lines Smean.section and Smax.section are uncommented.
# The function f_xyzt_smask works well on the xyt arrays of MLD variables, and the output from the
# ocean section subsetting (xt or yt dimensions) has to be further processed to compute period means or maxs before plotting
Smean <- list()
Smax <- list()
for (pname in names(periods)) {
  year1 <- periods[[pname]][1]
  year2 <- periods[[pname]][2]
  iperiod <- which(seq(1958,2019) %in% seq(year1,year2)) # note: only whole period, applied to ref exps only
  Smean[[pname]] <- lapply(Spre, function(A) {
    apply(A[,,iperiod], c(1,2), mean, na.rm=T)
  })
  Smax[[pname]] <- lapply(Spre, function(A) {
    apply(A[,,iperiod], c(1,2), max, na.rm=T)
  })
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Figure S16 B and C: changes in fluxes and concentrations across sections
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (fig_sections) {
  
  # Extract sections
  Fsection <- list()
  for (jsect in "LABNW") { # names(sectionmask)
    for (pname in names(Fpre)) {
      Fsection[[jsect]] <- lapply(Fpre, f_xyzt_smask, sectionmask[[jsect]])
    }
    
    # X axis for section plot: use longitude for zonal sections, else use latitude
    if (grepl("X", jsect)) {
      xs <- NASPG.lon[which(sectionmask[[jsect]][["zmat"]]==1)]
      xlab <- "Longitude"
    } else {
      xs <- NASPG.lat[which(sectionmask[[jsect]][["zmat"]]==1)]
      xlab <- "Latitude (ºN)"
    }
    iy <- which(deptht$zcenter < 3000)
    
    for (vv in c("poc")) { # varnames, c("poc","dissoc")
      
      # Calculate anomalies
      if (grepl("df",vv) | grepl("af",vv) | vv%in%c("thetao","so","uo","vo","wo")) {
        z1 <- Fsection[[jsect]][["weak"]][[vv]][,iy] - Fsection[[jsect]][["clim"]][[vv]][,iy]
        z2 <- Fsection[[jsect]][["strong"]][[vv]][,iy] - Fsection[[jsect]][["clim"]][[vv]][,iy]
        absmax <- max(abs(quantile(c(z1,z2),c(.01,.99),na.rm=T)))
        ll <- "Anomaly"
      } else {
        z1 <- 100 * (Fsection[[jsect]][["weak"]][[vv]][,iy] / Fsection[[jsect]][["clim"]][[vv]][,iy] - 1)
        z2 <- 100 * (Fsection[[jsect]][["strong"]][[vv]][,iy] / Fsection[[jsect]][["clim"]][[vv]][,iy] - 1)
        absmax <- max(abs(c(z1,z2)),na.rm=T)
        ll <- "% anomaly"
      }
      if (is.na(absmax) | abs(absmax)==Inf) {
        zlim <- c(-1,1) # safety net to avoid errors
      } else {
        zlim <- 35 * c(-1,1)
      }
      # Basic settings
      zback <- matrix(0, nrow = dim(z1)[1], ncol = dim(z1)[2]) # bathymetry background
      ncolors <- 14
      colors <- oce.colorsTwo(ncolors)
      breaks <- seq(zlim[1], zlim[2], length.out=ncolors+1)
      
      # Middle/Top: weak convection period
      png(filename = paste0(opath,"Fig_S16BC_",vv,"_section",jsect,"_yearmean_",
                            paste0(unlist(periods), collapse = "-"),"_",
                            expid,"_clim_",
                            paste0(periods$clim, collapse = "-"),"_weak.png"),
          width = 12, height = 7, units = 'cm', pointsize = 8, bg = "white", res = 600, type = "quartz")
      
      filled.contour(x = xs, y = deptht$zcenter[iy], z1, ylim=rev(range(deptht$zcenter[iy])), zlim=zlim, cex.lab=1.5, cex.axis=1.5,
                     levels = breaks, col = colors, xlab=xlab, ylab="Depth (m)")
      
      # https://stackoverflow.com/questions/28257533/r-plotting-a-scatterplot-over-a-filled-contour-plot
      nlevels <- ncolors
      zlim <- range(z1, finite = TRUE)
      las <- 1 
      levels <- pretty(zlim, nlevels)
      xlim <- range(xs, finite = TRUE)
      ylim <- rev(range(deptht$zcenter[iy], finite = TRUE))
      xaxs <- "i"
      yaxs <- "i"
      asp <- NA
      mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
      w <- (3 + mar.orig[2L]) * par("csi") * 2.54
      layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
      plot.window(ylim=ylim,xlim=xlim)
      
      lines(xs, Smean$weak$mlotstmax[which(sectionmask[[jsect]][["zmat"]]==1)], lwd=2, lty=2)
      lines(xs, Smax$weak$mlotstmax[which(sectionmask[[jsect]][["zmat"]]==1)], lwd=2)
      mtext("Weak convection (2009-2013)", line = 1, adj=0, cex = 1.5, col = "#90a1a3", font = 2)
      dev.off()
      
      # ----
      # Right/bottom: strong convection period
      png(filename = paste0(opath,"Fig_S16BC_",vv,"_section",jsect,"_yearmean_",
                            paste0(unlist(periods), collapse = "-"),"_",
                            expid,"_clim_",
                            paste0(periods$clim, collapse = "-"),"_strong.png"),
          width = 12, height = 7, units = 'cm', pointsize = 8, bg = "white", res = 600, type = "quartz")
      
      filled.contour(x = xs, y = deptht$zcenter[iy], z2, ylim=rev(range(deptht$zcenter[iy])), zlim=zlim, cex.lab=1.5, cex.axis=1.5,
                     levels = breaks, col = colors, xlab=xlab, ylab="Depth (m)")
      
      # https://stackoverflow.com/questions/28257533/r-plotting-a-scatterplot-over-a-filled-contour-plot
      nlevels <- ncolors
      zlim <- range(z2, finite = TRUE)
      las <- 1 
      levels <- pretty(zlim, nlevels)
      xlim <- range(xs, finite = TRUE)
      ylim <- rev(range(deptht$zcenter[iy], finite = TRUE))
      xaxs <- "i"
      yaxs <- "i"
      asp <- NA
      mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
      w <- (3 + mar.orig[2L]) * par("csi") * 2.54
      layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
      plot.window(ylim=ylim,xlim=xlim)
      
      lines(xs, Smean$strong$mlotstmax[which(sectionmask[[jsect]][["zmat"]]==1)], lwd=2, lty=2)
      lines(xs, Smax$strong$mlotstmax[which(sectionmask[[jsect]][["zmat"]]==1)], lwd=2)
      mtext("Strong convection (2014-2018)", line = 1, adj=0, cex = 1.5, col = "#2B6CBE", font = 2)
      dev.off()
      
    }
  }
  
}



# Partly based on analyze_float_matchups.R and analyze_plot_mergedProfilesTrajSat.R
# Dependencies: batch processing from analyze_plot_mergedProfilesTrajSat.R (merged profile/traj/sat data and convperiod.m files)
# See fig_conv_events.R
# Marti Gali Tapias, Dec 2025, marti.gali.tapias@gmail.com, mgali@icm.csic.es

# -----------------------------------------------------------------------------------------------------
# Input arguments
fwmoS <- c("6901480","6901486","6901523","6901524","6901527")

# Set to TRUE or FALSE to select the display item that you want to generate
paper_figure <- T        # Fig. 2 of convection paper
animation <- F           # Movie S1 of convection paper

# -----------------------------------------------------------------------------------------------------
# Define paths and names of input variable files
mpath <- "~/Desktop/Gali_2026_convectionPOC/input_data/"
spath <- "~/Desktop/Gali_2026_convectionPOC/input_data/"
gpath <- "~/Desktop/Gali_2026_convectionPOC/input_data/glorys12v1_daily/"
opath <- "~/Desktop/Gali_2026_convectionPOC/output/"             # merged profile, traj and satellite data with extra processing from analyze_plot_mergedProfilesTrajSat.R

if (animation) system( paste0("mkdir -p ",opath,"frames") )

# -----------------------------------------------------------------------------------------------------
# Libraries and functions
library(RNetCDF)
library(abind)
library(lubridate)
library(caTools)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(oce)
library(ocedata)
library(fields)
library(multcomp)
library(magick)
library(av)

# Custom functions
clamp <- function(x, minv, maxv) {
  pmin(pmax(x, minv), maxv)
}
bin_5day <- function(d) {
  y     <- lubridate::year(d)
  doy   <- data.table::yday(d) - 1           # 0-based day of year
  b_doy <- (doy %/% 5) * 5        # floored to 5-day bin
  lubridate::make_date(y, 1, 1) + b_doy
}
f_points_traj_chl <- function(dbdp, points_date, cbreaks, point_alpha, conv_red, halo) {
  jdates <- which(dbdp$DATE==points_date)
  if(sum(jdates)) {
    # print("Data found")
    bdata <- dbdp[jdates, c("fwmo","cperiods","DATE","DATEBIN","LONGITUDE","LATITUDE","CHLA_ADJUSTED")]
    bdata$CHLA_ADJUSTED <- clamp(bdata$CHLA_ADJUSTED, min(cbreaks), max(cbreaks))
    col_index <- cut(bdata$CHLA_ADJUSTED, breaks = cbreaks, include.lowest = TRUE, labels = FALSE)
    point_cols <- adjustcolor(cols[col_index], point_alpha)
    circle_cols <- rep("#FFFFFF30",length(bdata$cperiods))
    if (conv_red) circle_cols[bdata$cperiods=="conv"] <- "red"
    base_cex <- 1000^bdata$CHLA_ADJUSTED + 1
    base_cex[is.na(base_cex)] <- 2 # ;print(base_cex)
    if (halo) points(bdata$LONGITUDE, bdata$LATITUDE, pch=pch.f[bdata$fwmo], cex=base_cex+1, bg = "#FFFFFF30", col = NA)     # white halo
    points(bdata$LONGITUDE, bdata$LATITUDE, pch=pch.f[bdata$fwmo], cex=base_cex+0.1, bg = point_cols, col="white", lwd = 0.3)
    points(bdata$LONGITUDE, bdata$LATITUDE, pch=pch.f[bdata$fwmo], cex=base_cex+0.2, bg = NA, col = circle_cols)       # outer red circle
  }
}

# -----------------------------------------------------------------------------------------------------
# Load dbdp of each float (daily-binned traj data with additional period classification) and merge
dball <- lapply(fwmoS, function(fwmo) {
  load(paste0(mpath,"data_",fwmo,"_noTshift_9999_noreg.Rda"))
  return(dbdp)
})
names(dball) <- fwmoS

# Turn list into data frame
dbdp <- data.table::rbindlist(dball, use.names = T, fill = T, idcol = "fwmo")

# Add 5-day binned date centered on first day (matching OC-CCI 5D image dates)
dbdp$DATEBIN <- bin_5day(dbdp$DATE)
dbdp$DOY <- data.table::yday(dbdp$DATE)


# -----------------------------------------------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%% Map satellite chl along with BGC-Argo trajectory data %%%%%%%%%%%%%%%%%%%%%%%
# -----------------------------------------------------------------------------------------------------

# Satellite resolution
sres <- c("4.6km")

# Symbols
pch.f <- rep(21, 5) # c(15,16,18,4,17)
cex.f <- rep(3, 5) # c(2,2,2.3,1.8,2) - 0.2
cexout.f <- cex.f + 0.3 # c(2,2,2.3,1.8,2) + 0.4
cexhalo.f <- cex.f + 1
names(pch.f) <- fwmoS
names(cex.f) <- fwmoS
names(cexout.f) <- fwmoS
names(cexhalo.f) <- fwmoS

# Colors
cbreaks <- c(seq(0.02,0.2,0.02),0.30,0.40)
cols   <- colorRampPalette(viridis(9, direction=1))(length(cbreaks)-1)

# Define bounds and crop bathymetry to slightly smaller area to avoid overplotting plot box with the contours
rlon <- c(-61,-29) # v0 -62, -32
rlat <- c(53,65) # v0 53, 63
eps <- 0.05   # tweak (e.g. 0.02–0.1 depending on resolution)

# Bathymetry
globathy <- RNetCDF::read.nc( open.nc("~/Desktop/Gali_2026_convectionPOC/input_data/bathymetry_GLORYS12v1_NASPG.nc") ) # GLORYS12v1
bathy <- list(lonvec = globathy$longitude, latvec = globathy$latitude, zmat = -globathy$deptho)
bathy$zmat[bathy$zmat>0] <- NA
bathy$zmat[bathy$zmat<(-6000)] <- NA
eps <- 0.09   # tweak (e.g. 0.02–0.1 depending on resolution)
ilon <- which(bathy$lonvec >= (rlon[1] + eps) &
                bathy$lonvec <= (rlon[2] - eps))
ilat <- which(bathy$latvec >= (rlat[1] + eps) &
                bathy$latvec <= (rlat[2] - eps))
lon_crop <- bathy$lonvec[ilon]
lat_crop <- bathy$latvec[ilat]
zmat_crop <- bathy$zmat[ilon, ilat]

# Coastline
data("coastlineWorldMedium"); coastline <- coastlineWorldMedium


# ======================================================================
# ANIMATION: DAILY FLOAT DATA AND GLORYS MLD WITH 5D SATELLITE
# ======================================================================
if (animation) {
  
  # Plot parameters and map region
  pw <- 3400 # original 3400
  ph <- 1900 # original 1900
  
  # Selected year (ONLY 2015 INCLUDED IN REPOSITORY)
  for (sel_year in 2015) {
    sdateS <- seq(as.Date(paste0(sel_year,"-02-05")), as.Date(paste0(sel_year,"-04-21")), by="5 days") # 04-21 (better) or 26
    
    # GLORYS12v1 mixed layer depth (mlotst). Load Feb-Apr for selected year, keep in separate arrays
    months <- sprintf("%02d", 2:4)
    paths  <- file.path(gpath, paste0("mlotst_", sel_year, months, "_NASPG.nc"))
    ncfiles <- lapply(paths, open.nc)
    on.exit(lapply(ncfiles, close.nc), add = TRUE)
    glo <- c(
      list(
        lonvec = var.get.nc(ncfiles[[1]], "lon"),
        latvec = var.get.nc(ncfiles[[1]], "lat")
      ),
      setNames(
        lapply(ncfiles, var.get.nc, variable = "mlotst"),
        paste0("mlotst", months)
      )
    )
    
    for (sdate in sdateS) {
      
      sdate <- as.Date(sdate)
      # Satellite image
      fpath <- paste0(spath,"occciv6_chlos_5D_",sres,"/chlos_",gsub(pattern = "-", replacement = "", x = sdate),"_NASPG.nc")
      if (!file.exists(fpath)) stop("File not found")
      ncfile <- open.nc(fpath)
      cci <- list(
        lonvec = var.get.nc(ncfile, variable = "lon"),
        latvec = var.get.nc(ncfile, variable = "lat"),
        chlos = var.get.nc(ncfile, variable = "chlos")
      )
      close.nc(ncfile)
      # Corrections
      if (sres=="0.25deg") {
        cci$lon <- cci$lon - 360                      # shift longitude to degrees east, needed only for 0.25 deg data
      } else if (sres=="4.6km") {
        cci$lat <- rev(cci$lat)                       # make latitude increasing
        cci$chlos <- cci$chlos[, ncol(cci$chlos):1]   # flip lon-lat matrix along the latitude dimension
      }
      cci$chlos <- clamp(cci$chlos, min(cbreaks), max(cbreaks))
      
      # New loop on daily data within 5D satellite image
      ddateS <- seq(sdate, sdate+as.duration("4 days"), by="1 day")
      
      for (ddate in ddateS) {
        
        ddate <- as.Date(ddate)
        # Month (string) and day index for glorys
        mm <- as.character(format(ddate, "%m")) 
        iday <- as.integer(format(ddate, "%d"))
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        p <- paste0(opath,"frames/","fig_chl_sat_traj_mld_",ddate,"_",sres,".png")
        png(filename = p, width = pw, height = ph, pointsize = 8, bg = "white", res = 600, type = "cairo") # note that for animation plot size must be divisible by 2
        
        par(mar = c(3.3,4,2.7,8))
        
        # Satellite chl, coastline and contours
        image(x=cci$lon, y=cci$lat, z=cci$chlos, bg,
              xlab = "Longitude (º)", ylab = "Latitude (º)",
              main = paste0("DOY = ",data.table::yday(ddate),"  (",yday(ddateS[1]),"-",yday(ddateS[5]),")"),
              breaks = cbreaks, col = cols,
              xlim = rlon, ylim = rlat,
              mgp = c(2, 0.2, 0), tck=-0.01, cex.axis = 1)
        box()
        usr <- par("usr")
        clip(usr[1], usr[2], usr[3], usr[4])
        plot(coastline, clon = 0, clat = 0, span = c(length(bathy$lon), length(bathy$lon)),
             col = "gray50", bg = rgb(1,1,1,0), lwd = 0.5, axes = F, add = T)
        contour(lon_crop, lat_crop, zmat_crop, levels = c(-1000), lwd = c(2), lty = 1, labcex = 1, col = "gray70", drawlabels = F, add = T)
        contour(lon_crop, lat_crop, zmat_crop, levels = c(-2000), lwd = c(1), lty = 1, labcex = 1, col = "gray70", drawlabels = F, add = T)
        contour(glo$lonvec, glo$latvec, glo[[paste0("mlotst",mm)]][,,iday], levels = c(900), lwd = 0.5, lty = 1, labcex = 1, col = c("black"), drawlabels = F, add = T)
        contour(glo$lonvec, glo$latvec, glo[[paste0("mlotst",mm)]][,,iday], levels = 1000, lwd = 1.2, lty = 1, labcex = 1, col = "black", drawlabels = F, add = T)
        
        # Float trajectories with concentration as color scale, past days with increasing transparency to show trajectory streaks
        # f_points_traj_chl(dbdp, points_date = ddate - as.duration("1 day"), cbreaks, point_alpha = 0.6, conv_red = F, halo = F) # previous day points
        f_points_traj_chl(dbdp, points_date = ddate, cbreaks, point_alpha = 1.0, conv_red = T, halo = F) # same day points
        
        legend("topleft",
               inset = 0.008,
               cex = 1.2, title.cex = 1.5,
               title = expression(paste("Daily MLD"[0.03])),
               legend = c(
                 "950 m",
                 "1000 m"),
               lwd = c(1,1),
               lty = c(1,3),
               col = c("black","black"),
               seg.len = rep(1,4),
               bg = "#FFFFFFE6",
               box.col = "#FFFFFFE6"
        )
        
        # Colorbar for chlorophyll
        par( fig=c(.86,.88,.15,.85), new=TRUE, mar=c(0,0,0,0) )
        cbar.int <- diff(range(cbreaks))/(length(cbreaks)-1)
        cbar.chl <- seq(min(cbreaks), max(cbreaks), cbar.int)
        cbar.tick <- cbar.chl - cbar.int/2
        cbar.labels <- as.character(cbreaks)
        cbar.labels[length(cbar.labels)] <- paste0("> ", cbreaks[length(cbreaks)])
        image(t(cbar.chl), cbar.chl, col = cols, xaxt = "n", yaxt = "n", bg = "white")
        axis(4, cex.axis=1.5, mgp = c(0, 0.6, 0), at = cbar.tick, labels = cbar.labels, las = 1, lwd = 0.5, cex.axis = 1)
        mtext(side = 4, text = expression(paste("Chl",italic("a")," (mg ",m^-3,")")), las = 0, line = 3.5, cex = 1.5)
        box(lwd = 0.5)
        
        dev.off()
        
      } # loop on days within 5D period
    } # loop on satellite 5D periods
    
    files <- list.files(
      paste0(opath, "frames"),
      pattern = "\\.png$",
      full.names = TRUE
    )
    
    # files <- files[order(as.numeric(gsub("\\D", "", files)))]
    
    av_encode_video(files,
                    output=paste0(opath, "Movie_S1_", sel_year, ".mp4"),
                    framerate=2)
    
  } # loop on years
}


# ======================================================================
# PAPER FIGURE
# ======================================================================
if (paper_figure) {
  
  
  # p <- paste0(opath,"Fig_2_chl_sat_traj_mld_",sres,".png")
  # png(filename = p, width = 17, height = 10, units = 'cm', pointsize = 8, bg = "white", res = 600, type = "quartz")
  
  p <- paste0(opath,"Fig_2_chl_sat_traj_mld_",sres,".tiff")
  tiff(filename = p, width = 17, height = 10, units = 'cm', pointsize = 8, bg = "white", res = 300, compression = "lzw")
  
  # Layout: 2x2 panels and right colorbar
  sm1 <- matrix(data = 1, nrow = 6, ncol = 9)
  sm2 <- matrix(data = 2, nrow = 6, ncol = 8)
  sm3 <- matrix(data = 3, nrow = 6, ncol = 9)
  sm4 <- matrix(data = 4, nrow = 6, ncol = 8)
  sm5 <- matrix(data = 5, nrow = 12, ncol = 2) # right column for colorbar
  
  layout( cbind(rbind(cbind(sm1, sm2), cbind(sm3, sm4)), sm5) )
  par(oma = c(1,1,0,0))
  
  sdateS <- as.Date(c("2015-03-12",
                      "2015-03-17",
                      "2015-04-01",
                      # "2015-04-06",
                      "2015-04-11")) # v0 just this panel
  
  # v0: c(4,4,2,7))
  pmargins <- list( c(2.3,4,1.9,2),
                    c(2.3,0,1.9,2),
                    c(2.9,4,1.3,2),
                    c(2.9,0,1.3,2))
  
  j <- 0
  
  for (sdate in sdateS) {
    
    sdate <- as.Date(sdate)
    j <- j + 1
    # Satellite image
    fpath <- paste0(spath,"occciv6_chlos_5D_",sres,"/chlos_",gsub(pattern = "-", replacement = "", x = sdate),"_NASPG.nc")
    ncfile <- open.nc(fpath)
    if (!file.exists(fpath)) stop("File not found")
    cci <- list(
      lonvec = var.get.nc(ncfile, variable = "lon"),
      latvec = var.get.nc(ncfile, variable = "lat"),
      chlos = var.get.nc(ncfile, variable = "chlos")
    )
    close.nc(ncfile)
    # Corrections
    if (sres=="0.25deg") {
      cci$lon <- cci$lon - 360                      # shift longitude to degrees east, needed only for 0.25 deg data
    } else if (sres=="4.6km") {
      cci$lat <- rev(cci$lat)                       # make latitude increasing
      cci$chlos <- cci$chlos[, ncol(cci$chlos):1]   # flip lon-lat matrix along the latitude dimension
    }
    cci$chlos <- clamp(cci$chlos, min(cbreaks), max(cbreaks))
    
    # GLORYS12v1 mixed layer depth (mlotst). Concatenate current and posterior month to ensure 5-days subset possible and previous few days as well
    gdate <- unlist(strsplit(as.character(sdate), split = "-"))
    rpath1 <- paste0(gpath,"mlotst_",gdate[1],gdate[2],"_NASPG.nc")
    rpath2 <- paste0(gpath,"mlotst_",gdate[1],sprintf("%02i", as.numeric(gdate[2])+1),"_NASPG.nc")
    if (!file.exists(rpath1) | !file.exists(rpath2)) stop("File not found")
    ncfile1 <- open.nc(rpath1); ncfile2 <- open.nc(rpath2)
    glo <- list(
      lonvec = var.get.nc(ncfile1, variable = "lon"),
      latvec = var.get.nc(ncfile1, variable = "lat"),
      mlotst1 = var.get.nc(ncfile1, variable = "mlotst"),
      mlotst2 = var.get.nc(ncfile2, variable = "mlotst")
    )
    close.nc(ncfile1); close.nc(ncfile2)
    
    ind_5D <- seq(as.numeric(gdate[3]), as.numeric(gdate[3])+4)
    glo$mlotst <- abind(glo$mlotst1, glo$mlotst2, along = 3)[,,ind_5D]
    glo$mlotst1 <- NULL
    glo$mlotst2 <- NULL
    glo$mlotst_max <- apply(glo$mlotst, MARGIN = c(1,2), max)
    glo$mlotst_min <- apply(glo$mlotst, MARGIN = c(1,2), min)
    
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT PANEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    par(mar = pmargins[[j]])
    
    # Satellite chl
    image(x=cci$lon, y=cci$lat, z=cci$chlos, xlab = "", ylab = "", breaks = cbreaks, col = cols,
          xlim = rlon, ylim = rlat,
          axes = F,
          cex.main = 1.4, adj = 0.5)
    title(main = paste0("DOY ",data.table::yday(sdate),"-",data.table::yday(sdate)+4),
          line = 0.3,
          cex.main = 1.4, adj = 0.5)
    # Panel-dependent axes and ticks
    if (j == 1) {
      axis(side = 1, mgp = c(2, 0.4, 0), tck=-0.01, cex.axis = 1, labels = F)
      axis(side = 2, mgp = c(2, 0.4, 0), tck=-0.01, cex.axis = 1)
      mtext(side = 2, text = "Latitude (º)", cex = 1, adj = -0.5, line = 2.2)
    }
    if (j == 2) {
      axis(side = 1, mgp = c(2, 0.4, 0), tck=-0.01, cex.axis = 1, labels = F)
      axis(side = 2, mgp = c(2, 0.4, 0), tck=-0.01, cex.axis = 1, labels = F)
    }
    if (j == 3) {
      axis(side = 1, mgp = c(2, 0.4, 0), tck=-0.01, cex.axis = 1)
      axis(side = 2, mgp = c(2, 0.4, 0), tck=-0.01, cex.axis = 1)
      mtext(side = 1, text = "Longitude (º)", cex = 1, adj = 1.2, line = 2.5)
    }
    if (j == 4) {
      axis(side = 1, mgp = c(2, 0.4, 0), tck=-0.01, cex.axis = 1)
      axis(side = 2, mgp = c(2, 0.4, 0), tck=-0.01, cex.axis = 1, labels = F)
    }
    box()
    
    # Coastline and contours
    usr <- par("usr")
    clip(usr[1], usr[2], usr[3], usr[4])
    plot(coastline, clon = 0, clat = 0, span = c(length(bathy$lon), length(bathy$lon)),
         col = "gray50", bg = rgb(1,1,1,0), lwd = 0.5, axes = F, add = T)
    contour(lon_crop, lat_crop, zmat_crop, levels = c(-1000), lwd = c(1), lty = 1, labcex = 0.8, col = "gray70", drawlabels = T, add = T)
    contour(lon_crop, lat_crop, zmat_crop, levels = c(-2000), lwd = c(0.5), lty = 1, labcex = 0.8, col = "gray70", drawlabels = F, add = T)
    contour(glo$lonvec, glo$latvec, glo$mlotst_max, levels = 900, lwd = 0.5, lty = 1, labcex = 1, col = c("black"), drawlabels = F, add = T)
    contour(glo$lonvec, glo$latvec, glo$mlotst_min, levels = 1000, lwd = 1.2, lty = 1, labcex = 1, col = "black", drawlabels = F, add = T)
    text(x = rlon[2]-2, y = rlat[2]-0.7, LETTERS[j], cex = 2.5, font = 2, adj = 0)
    
    # Float trajectories with concentration as color scale
    jdates <- which(dbdp$DATEBIN==sdate)
    if(sum(jdates)) {
      bdata <- dbdp[jdates, c("fwmo","cperiods","DATE","DATEBIN","LONGITUDE","LATITUDE","CHLA_ADJUSTED")]
      bdata$CHLA_ADJUSTED <- clamp(bdata$CHLA_ADJUSTED, min(cbreaks), max(cbreaks))
      col_index <- cut(bdata$CHLA_ADJUSTED, breaks = cbreaks, include.lowest = TRUE, labels = FALSE)
      point_cols <- cols[col_index]
      circle_cols <- rep("#FFFFFF30",length(bdata$cperiods))
      circle_cols[bdata$cperiods=="conv"] <- "red"
      base_cex <- 1000^bdata$CHLA_ADJUSTED + 1
      base_cex[is.na(base_cex)] <- 2 # ;print(base_cex)
      # points(bdata$LONGITUDE, bdata$LATITUDE, pch=pch.f[bdata$fwmo], cex=base_cex+0.5, bg = "#FFFFFF80", col=NA)     # white halo
      points(bdata$LONGITUDE, bdata$LATITUDE, pch=pch.f[bdata$fwmo], cex=base_cex+0.5, col = circle_cols)       # outer red circle
      points(bdata$LONGITUDE, bdata$LATITUDE, pch=pch.f[bdata$fwmo], cex=base_cex+0.1, bg = cols[col_index], col="white", lwd = 0.3) # 0.4
      
      # Annotate floats and MLD contours (legend)
      if (j %in% c(1,4)) {
        loff <- list(
          c(0,-0.8),
          c(3,1),
          c(2.5,-0.8),
          c(3.5,-0.5),
          c(0,0)
        )
        names(loff) <- fwmoS
        
        for (fwmo in fwmoS) {
          jfd <- which(bdata$fwmo == fwmo &
                         bdata$DATE == max(bdata$DATE[bdata$fwmo == fwmo]))
          
          if (length(jfd) > 0) {
            
            x <- bdata$LONGITUDE[jfd] + loff[[fwmo]][1]
            y <- bdata$LATITUDE[jfd] + loff[[fwmo]][2]
            
            lab <- fwmo
            
            # text size
            tcex <- 1
            w <- strwidth(lab, cex = tcex)
            h <- strheight(lab, cex = tcex)
            
            # padding
            pad_x <- 0.1 * w
            pad_y <- 0.4 * h
            
            # draw semi-transparent white box
            rect(x - w/2 - pad_x, y - h/2 - pad_y,
                 x + w/2 + pad_x, y + h/2 + pad_y,
                 col = adjustcolor("grey99", alpha.f = 0.9), border = "grey90")
            
            # draw text on top
            text(x, y, labels = lab, cex = tcex, font = 1)
          }
        }
      }
    }
    if (j==1) {
      legend("topleft",
             inset = 0.008,
             # x = rlon[1] + diff(rlon)*0.01,
             # y = rlat[2] -  diff(rlat)*0.01,
             cex = 1.2, title.cex = 1.5,
             title = expression(paste("MLD"[0.03])),
             legend = c(
               " < 900 m",
               " > 1000 m"),
             lwd = c(0.5,1),
             lty = c(1.2,1),
             col = c("black","black"),
             seg.len = rep(1,4),
             bg = "#FFFFFFE6",
             box.col = "grey90")
    }
    
  }
  
  # Right column: Colorbar for chlorophyll
  plot.new()
  
  par( fig=c(.90,.92,.10,.90), new=TRUE, mar=c(0,0,0,0) )
  cbar.int <- diff(range(cbreaks))/(length(cbreaks)-1)
  cbar.chl <- seq(min(cbreaks), max(cbreaks), cbar.int)
  cbar.tick <- cbar.chl - cbar.int/2
  cbar.labels <- as.character(cbreaks)
  cbar.labels[length(cbar.labels)] <- paste0("> ", cbreaks[length(cbreaks)])
  image(t(cbar.chl), cbar.chl, col = cols, xaxt = "n", yaxt = "n", bg = "white")
  axis(4, cex.axis=1.5, mgp = c(0, 0.6, 0), at = cbar.tick, labels = cbar.labels, las = 1, lwd = 0.5)
  mtext(side = 4, text = expression(paste("Chl",italic("a")," (mg ",m^-3,")")), las = 0, line = 5, cex = 1.4)
  box(lwd = 0.5)
  
  dev.off()
}


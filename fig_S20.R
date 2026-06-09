# DISCLAIMER -------------------------------------------------------------------
# This script reproduces Figure S20 from the manuscript:
# "Convection injects labile particulate organic carbon to the deep ocean", by
#
# Martí Galí, María Andrea Orihuela-García, Yohan Ruprich-Robert,
# Vladimir Lapin, María Sánchez-Urrea, Marcos Fontela,
# Joan Llort, Valentina Sicardi, Raffaele Bernardello.
#
# IMPORTANT:
# The datasets used in this script are derived from the original data sources:
# https://doi.org/10.25921/s4f4-ye35
# http://doi.org/10.20350/DIGITALCSIC/8513
#
# Users MUST cite and acknowledge these original datasets when using this code
# or any derived results in publications or presentations.
#
# The data provided here represent only a geographically filtered subset of the
# original datasets and do not replace or supersede them.
#
# Script author:
# Marcos Fontela (mfontela@iim.csic.es)
#
# LICENSE (MIT) ----------------------------------------------------------------
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# -----------------------------------------------------------------------------

pkgs <- c("tidyverse", "oce", "metR", "patchwork", "marmap", "plyr", "lubridate")
install.packages(pkgs[!pkgs %in% rownames(installed.packages())]); rm(pkgs)

library(tidyverse)
library(oce)
library(metR)
library(patchwork)

# data loading ---------------------------------------------------------

d <- readRDS(
  "~/Desktop/Gali_2026_convectionPOC/input_data/Gali_POC_convection_FigS20.rds") %>%
  mutate(period = ifelse(
    between(lubridate::year(date), 2014, 2018), "Strong", "Weak"
  ))

d_bin <- d %>%
  filter(latitude < 63.5, !is.na(period)) %>%
  mutate(
    swSigmaTheta = oce::swSigmaTheta(ctd_salinity, ctd_temperature, ctd_pressure),
    sal_bin      = plyr::round_any(ctd_salinity, .01),
    theta_bin    = plyr::round_any(theta, .1),
    period     = factor(period, levels = c("Weak", "Strong"))
  ) %>%
  group_by(sal_bin, theta_bin, period) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# map -------------------------------------------------------------------

coordinates <- data.frame(
  lon = c(-65, -30, -20, -42, -65, -65),
  lat = c(65, 70, 65, 50, 50, 65)
)

m <- marmap::getNOAA.bathy(
  lon1 = -65,
  lon2 = -15,
  lat1 = 45,
  lat2 = 70,
  resolution = 5
)

basemap <- ggplot(m, aes(x = x, y = y)) +
  coord_quickmap() +
  geom_raster(aes(fill = z)) +
  scale_fill_gradient2(
    low = "cadetblue1",
    mid = "white",
    high = "darkgreen",
    midpoint = 0,
    guide = FALSE
  ) +
  geom_contour(
    aes(z = z),
    breaks = c(-1000, -2000, -3000, -4000, -5000),
    colour = "gray",
    size = 0.5
  ) +
  geom_contour(
    aes(z = z),
    breaks = c(-2500),
    colour = "gray63",
    size = 0.7
  ) +
  geom_contour(
    aes(z = z),
    breaks = c(0),
    colour = "black",
    size = 0.3
  ) +
  geom_point(
    data = coordinates,
    aes(x = lon + 0.25, y = lat - 0.25),
    color = "black",
    size = 4,
    shape = 13
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "", y = "") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 2),
    text = element_text(size = 20, colour = "black"),
    legend.title = element_blank()
  )

map_plot <- basemap +
  geom_point(
    data = d %>%
      mutate(period = factor(.data$period, levels = c("Weak", "Strong"))) %>%
      filter(latitude > 52, .data$period == "Strong"),
    aes(longitude, latitude, colour = as.character(period)),
    size = 4
  ) +
  geom_point(
    data = d %>%
      mutate(period = factor(.data$period, levels = c("Weak", "Strong"))) %>%
      filter(latitude > 52, .data$period != "Strong"),
    aes(longitude, latitude, colour = as.character(period)),
    size = 3
  ) +
  scale_color_manual(values = c("green4", "lightgreen")) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 14),
    plot.margin = margin(0, 0, 0, 0),
    axis.title.x = element_text(margin = margin(t = 0)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )

rm(coordinates, basemap, m)
# TS plot -------------------------------------------------------------

TSgrid <- expand.grid(ctd_salinity = seq(34.5, 35.3, 0.01), theta = seq(0, 8.5, 0.1)) %>%
  mutate(swSigmaTheta = oce::swSigmaTheta(ctd_salinity, theta, p = 0))

SWT_table <- data.frame(
  wm    = c("IrSPMW", "LSW",   "ISOW",  "DSOW"),
  theta = c(5.0,      3.4,     2.7,     1.3),
  Sal   = c(35.01,    34.86,   35.00,   34.905)
)

TS_plot <- ggplot(d_bin, aes(sal_bin, theta_bin)) +
  geom_contour(
    data = TSgrid,
    aes(x = ctd_salinity, y = theta, z = swSigmaTheta),
    colour = "grey40",
    size = 0.4,
    breaks = seq(27.4, 28.6, 0.05)
  ) +
  metR::geom_text_contour(
    data = TSgrid,
    aes(
      x = ctd_salinity,
      y = theta,
      z = swSigmaTheta,
      label = after_stat(level)
    ),
    breaks = seq(27.4, 28.6, 0.05),
    stroke = 0.2,
    size = 2.5,
    colour = "grey20"
  ) +
  geom_point(
    aes(colour = ctd_pressure),
    shape = 15,
    size = 3,
    alpha = 0.8
  ) +
  geom_label(
    data = SWT_table,
    aes(x = Sal, y = theta, label = wm),
    colour = "gray33",
    inherit.aes = FALSE,
    size = 2
  ) +
  facet_grid(. ~ period) +
  scale_x_continuous(limits = c(34.65, 35.05)) +
  scale_y_continuous(limits = c(1, 8.5)) +
  scale_colour_viridis_c(direction = -1) +
  labs(
    colour = "",
    x = "Salinity",
    y = expression(theta ~ (degree * C))
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 2),
    text = element_text(size = 13, colour = "black"),
    legend.title = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    axis.title.x = element_text(margin = margin(t = 0)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )

rm(SWT_table, TSgrid)
# vertical profile --------------------------------------------------------


profile_plot <- ggplot(
  d_bin %>%
    mutate(new_press = plyr::round_any(ctd_pressure, 200)) %>%
    group_by(new_press, station, period) %>%
    summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)), .groups = "drop"),
  aes(doc, -ctd_pressure)
) +
  geom_point(
    aes(colour = period),
    size = 4,
    alpha = 0.6
  ) +
  geom_smooth(
    aes(group = period, colour = period),
    show.legend = FALSE,
    size = 2,
    method = "loess",
    se = FALSE,
    span = 0.5,
    orientation = "y"
  ) +
  scale_x_continuous(
    position = "top",
    limits = c(43, 70),
    breaks = seq(45, 70, 5)
  ) +
  scale_y_continuous(
    breaks = c(-4000, -3000, -2000, -1000, 0),
    labels = c(4000, 3000, 2000, 1000, 0)
  ) +
  scale_color_manual(values = c("lightgreen", "green4")) +
  labs(
    colour = "",
    x = "DOC (μmol·kg⁻¹)",
    y = "Pressure (dbar)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 19),
    plot.margin = margin(0, 0, 0, 0),
    axis.title.x = element_text(margin = margin(t = 0)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )

# final composition ------------------------------------------------------

FigS20 <- ((map_plot / TS_plot) | profile_plot) +
  plot_layout(widths = c(3, 1))

ggsave(
  filename = "~/Desktop/Gali_2026_convectionPOC/output/Fig_S20b.png",
  plot = FigS20,
  width = 12,
  height = 8,
  dpi = 300
)


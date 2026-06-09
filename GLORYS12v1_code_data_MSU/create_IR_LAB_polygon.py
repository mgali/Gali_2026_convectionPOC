###############################################################################
#  Create IR_LAB Mask
###############################################################################
# ---------------------------------------------------------------------------
# 1. Load template and bathymetry
# 2. Create regional mask (IR_LAB)
# 3. Save results
# ---------------------------------------------------------------------------
# Author    : Maria Sanchez Urrea
# Email     : mariasu@icm.csic.es
# Created   : November 2024
###############################################################################

# User options ---------------------------------------------------------------
show_plots = True    # Set to False to disable plotting

# Libraries ------------------------------------------------------------------
import xarray as xr
import numpy as np
import regionmask
import geopandas as gpd

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cmocean import cm

import warnings
warnings.filterwarnings('ignore')

###############################################################################
# 1. Load Data
###############################################################################

# Paths ----------------------------------------------------------------------
path_data = 'D:/Arctic/projects/deep_convection_NA/data/'
path_save = 'D:/Arctic/projects/deep_convection_NA/data/'

# GLORYS12v1 template --------------------------------------------------------
template_file = path_data + 'template_rectangule_glorys12v1.nc'
dset = xr.open_dataset(template_file)

mask_raw = dset.mlotst
mask = xr.where(np.isnan(mask_raw), 1, 1)  # uniform mask (all ones)
dset.close()

# IR_LAB shapefile -----------------------------------------------------------
shp_path = path_data + "IR_LAB_shapefile/"
shp_file = "IR_LAB_budgets.shp"
IR_LAB_shape = gpd.read_file(shp_path + shp_file)

# GLORYS12v1 bathymetry ------------------------------------------------------
bathy_file = path_data + 'GLO-MFC_001_030_mask_bathy.nc'
dset = xr.open_dataset(bathy_file)

lonW, lonE, latS, latN = -65, -15, 40, 70
bathy = dset.deptho.isel(
    longitude=(dset.longitude >= lonW) & (dset.longitude <= lonE),
    latitude=(dset.latitude >= latS) & (dset.latitude <= latN)
)

# Keep only deep ocean (>2500 m)
bathy = xr.where(bathy >= 2500, bathy, np.nan)

###############################################################################
# 2. Create IR_LAB Mask
###############################################################################

# Bathymetry mask ------------------------------------------------------------
mask_bathy = xr.where(~np.isnan(bathy), mask, np.nan)

# IR_LAB polygon (manual definition) -----------------------------------------
# Coordinates:
# lon: -65, -30, -20, -42, -65, -65
# lat:  65,  70,  65,  50,  50,  65

IR_LAB_coords = np.array([
    [-65, 65],
    [-30, 70],
    [-20, 65],
    [-42, 50],
    [-65, 50],
    [-65, 65]
])

polygon_lab = regionmask.Regions(
    [IR_LAB_coords],
    names=['Labrador Sea'],
    abbrevs=['IR_LAB'],
    name='IR_LAB'
)

# Optional: plot mask and polygon --------------------------------------------
if show_plots:
    ax = plt.subplot(111, projection=ccrs.PlateCarree())
    mask_bathy.plot(
        transform=ccrs.PlateCarree(),
        cmap=cm.haline,
        add_colorbar=False
    )
    ax.coastlines()
    IR_LAB_shape.plot(ax=ax, transform=ccrs.PlateCarree(),
                      color='blue', edgecolor='blue', alpha=0.7)
    polygon_lab.plot_regions(ax=ax, add_label=False)

# Convert polygon to xarray mask ---------------------------------------------
mask_region = polygon_lab.mask(
    mask_bathy.longitude.data,
    mask_bathy.latitude.data
).rename({"lat": "latitude", "lon": "longitude"})

# Final IR_LAB mask (polygon + bathymetry) -----------------------------------
IR_LAB_mask = xr.where(mask_region == 0, mask_bathy.squeeze(), np.nan).squeeze()

###############################################################################
# 3. Save Output
###############################################################################

ds = xr.DataArray(
    data=IR_LAB_mask.data,
    dims=["latitude", "longitude"],
    coords=dict(
        longitude=IR_LAB_mask.longitude.data,
        latitude=IR_LAB_mask.latitude.data,
    ),
    attrs=dict(description="Labrador Sea mask"),
    name="IR_LAB"
)

# Metadata -------------------------------------------------------------------
ds['longitude'].attrs.update({
    'standard_name': 'longitude',
    'long_name': 'longitude',
    'units': 'degrees_east',
    'axis': 'X',
    '_CoordinateAxisType': 'Lon'
})

ds['latitude'].attrs.update({
    'standard_name': 'latitude',
    'long_name': 'latitude',
    'units': 'degrees_north',
    'axis': 'Y',
    '_CoordinateAxisType': 'Lat'
})

# Save file ------------------------------------------------------------------
ds.to_netcdf(path_save + "mask_IR_LAB.nc", 'w', format='NETCDF4_CLASSIC')

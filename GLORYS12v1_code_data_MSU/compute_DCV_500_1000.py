###############################################################################
#  Deep Convection Volume (DCV) in the IR_LAB Region
###############################################################################
# ---------------------------------------------------------------------------
# 1. Load daily MLD data
# 2. Compute yearly maximum MLD
# 3. Apply IR_LAB mask and compute convected volume
# 4. Save results as CSV
# ---------------------------------------------------------------------------
# Author    : Maria Sanchez Urrea
# Email     : maria.sanchez.urrea@gmail.com
# Created   : November 2024
###############################################################################

import os
import xarray as xr
import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings('ignore')

###############################################################################
# 1. User settings
###############################################################################

year_selection = list(np.arange(1993, 2019 + 1))

# Paths
path_data   = "D:/Arctic/projects/deep_convection_NA/data/daily_data/"
path_mask   = "D:/Arctic/projects/deep_convection_NA/data/"
path_coords = "D:/Arctic/projects/deep_convection_NA/data/"
path_save   = "D:/Arctic/projects/deep_convection_NA/outputs/"

###############################################################################
# 2. Load daily MLD data (all years)
###############################################################################

all_files = os.listdir(path_data)

files = [
    f for f in all_files
    if f.startswith("mlotst_") and any(f"_{year}" in f for year in year_selection)
]

files_with_path = [os.path.join(path_data, f) for f in files]

dset = xr.open_mfdataset(files_with_path, chunks="auto")
mld = dset.mlotst
dset.close()

###############################################################################
# 3. Load IR_LAB mask
###############################################################################

mask_file = path_mask + "mask_IR_LAB.nc"
dset_mask = xr.open_dataset(mask_file)
mask = dset_mask.IR_LAB
dset_mask.close()

###############################################################################
# 4. Load GLORYS12v1 grid-cell area (m²)
###############################################################################

coords_file = path_coords + "GLO-MFC_001_030_coordinates.nc"
dset_coord = xr.open_dataset(coords_file)

cell_area = dset_coord.e1t * dset_coord.e2t   # m²
dset_coord.close()

# Subset study region
lonW, lonE = -65, -15
latS, latN = 40, 70

cell_area = cell_area.isel(
    longitude=(cell_area.longitude >= lonW) & (cell_area.longitude <= lonE),
    latitude=(cell_area.latitude >= latS) & (cell_area.latitude <= latN)
)

###############################################################################
# 5. Compute yearly maximum MLD
###############################################################################

mld_yearly_max = mld.groupby("time.year").max("time", keep_attrs=True)

###############################################################################
# 6. Apply IR_LAB mask
###############################################################################

mld_yearly_max_masked = mld_yearly_max.where(mask == 1, np.nan)
area_masked = cell_area.where(mask == 1, np.nan)

###############################################################################
# 7. Compute Deep Convection Volume (m³)
###############################################################################

# --- Full column convection (0–MLDmax) -------------------------------------
dcv_full = mld_yearly_max_masked * area_masked
dcv_full_yearly = dcv_full.sum(["latitude", "longitude"], skipna=True)

# --- Convection between 500 and 2000 m -------------------------------------
mld_cap_2000 = xr.where(mld_yearly_max_masked > 2000, 2000, mld_yearly_max_masked)
mld_500_2000 = mld_cap_2000 - 500
mld_500_2000 = xr.where(mld_500_2000 < 0, 0, mld_500_2000)

dcv_500_2000 = mld_500_2000 * area_masked
dcv_500_2000_yearly = dcv_500_2000.sum(["latitude", "longitude"], skipna=True)

dcv_500_2000_yearly.name = "DCV_500_2000"

# Convert to DataFrame
df_dcv = dcv_500_2000_yearly.to_dataframe().reset_index()

###############################################################################
# 8. Save CSV
###############################################################################

output_csv = path_save + "deep_convection_volume_annual_IR_LAB_500_2000.csv"

df_dcv.to_csv(
    output_csv,
    encoding="ISO-8859-1",
    sep=",",
    decimal=".",
    index=False
)

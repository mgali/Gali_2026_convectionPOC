###############################################################################
#  Compute Yearly Area of Deep Mixed Layer (MLD ≥ 500, 1000, 1500, 2000 m)
###############################################################################
# ---------------------------------------------------------------------------
# 1. Load grid-cell area from GLORYS12v1
# 2. For each year:
#       - Load daily MLD
#       - Compute yearly maximum MLD
#       - Compute area exceeding selected MLD thresholds
# 3. Save results as CSV
# ---------------------------------------------------------------------------
# Author    : Maria Sanchez Urrea
# Email     : mariasu@icm.csic.es
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

start_year = 1993
end_year   = 2019
year_selection = list(np.arange(start_year, end_year + 1))

# Paths
path_coords = 'D:/Arctic/projects/deep_convection_NA/data/'
path_data   = 'D:/Arctic/projects/deep_convection_NA/data/daily_data/'
path_save   = 'D:/Arctic/projects/deep_convection_NA/outputs/'

###############################################################################
# 2. Load GLORYS12v1 grid-cell area
###############################################################################

coords_file = path_coords + "GLO-MFC_001_030_coordinates.nc"
dset = xr.open_dataset(coords_file)

# e1t = grid spacing in x-direction (m)
# e2t = grid spacing in y-direction (m)
cell_area = (dset.e1t * dset.e2t) * 1e-6   # convert to km²
dset.close()

###############################################################################
# 3. Loop over years and compute MLD area above thresholds
###############################################################################

df_mld_area = pd.DataFrame()

for yy in year_selection:

    # -----------------------------------------------------------------------
    # Load daily MLD for the given year
    # -----------------------------------------------------------------------
    mld_file = path_data + f"mlotst_{yy}.nc"
    dset = xr.open_dataset(mld_file)
    mld = dset.mlotst
    dset.close()

    # -----------------------------------------------------------------------
    # Compute yearly maximum MLD
    # -----------------------------------------------------------------------
    mld_year_max = mld.groupby('time.year').max('time', keep_attrs=True)
    del mld

    # -----------------------------------------------------------------------
    # Compute area where MLD exceeds thresholds
    # -----------------------------------------------------------------------
    thresholds = [500, 1000, 1500, 2000]
    area_results = {}

    for th in thresholds:
        mld_mask = mld_year_max.where(mld_year_max >= th)
        area_mask = xr.where(~np.isnan(mld_mask), 1, np.nan) * cell_area
        area_sum = area_mask.sum(['latitude', 'longitude'], skipna=True)
        area_results[f"MLDmax_{th}"] = float(area_sum.values)

    # -----------------------------------------------------------------------
    # Store results for this year
    # -----------------------------------------------------------------------
    df_mld_area = pd.concat([
        df_mld_area,
        pd.DataFrame({
            "year": [yy],
            "MLDmax_500":  [area_results["MLDmax_500"]],
            "MLDmax_1000": [area_results["MLDmax_1000"]],
            "MLDmax_1500": [area_results["MLDmax_1500"]],
            "MLDmax_2000": [area_results["MLDmax_2000"]],
        })
    ], ignore_index=True)

###############################################################################
# 4. Save CSV
###############################################################################

output_csv = path_save + "mlotst_yearly_max_area_500_1000_1500_2000_rectangule.csv"

df_mld_area.to_csv(
    output_csv,
    encoding="ISO-8859-1",
    sep=',',
    decimal='.',
    index=False
)

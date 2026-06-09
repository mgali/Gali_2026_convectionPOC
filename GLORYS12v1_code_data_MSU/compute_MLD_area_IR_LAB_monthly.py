###############################################################################
#  Monthly Mean and Maximum MLD over the IR_LAB Region
###############################################################################
# ---------------------------------------------------------------------------
# 1. Load IR_LAB mask
# 2. For each year:
#       - Load daily MLD
#       - Compute monthly mean and max
#       - Apply IR_LAB mask
#       - Compute spatial mean over the region
# 3. Save results as CSV
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

start_year = 1993
end_year   = 2019
year_selection = list(np.arange(start_year, end_year + 1))

# Paths
path_data = "D:/Arctic/projects/deep_convection_NA/data/daily_data/"
path_mask = "D:/Arctic/projects/deep_convection_NA/data/"
path_save = "D:/Arctic/projects/deep_convection_NA/outputs/"

###############################################################################
# 2. Load IR_LAB mask
###############################################################################

mask_file = path_mask + "mask_IR_LAB.nc"
dset = xr.open_dataset(mask_file)
mask = dset.IR_LAB
dset.close()

###############################################################################
# 3. Loop over years and compute monthly mean/max MLD
###############################################################################

df_mean_max = pd.DataFrame()

for yy in year_selection:

    # -----------------------------------------------------------------------
    # Load daily MLD for the given year
    # -----------------------------------------------------------------------
    mld_file = path_data + f"mlotst_{yy}.nc"
    dset = xr.open_dataset(mld_file)
    mld = dset.mlotst
    dset.close()

    # -----------------------------------------------------------------------
    # Compute monthly mean and max
    # -----------------------------------------------------------------------
    mld_monthly_mean = mld.groupby('time.month').mean('time', keep_attrs=True)
    mld_monthly_max  = mld.groupby('time.month').max('time', keep_attrs=True)
    del mld

    # -----------------------------------------------------------------------
    # Apply IR_LAB mask
    # -----------------------------------------------------------------------
    mld_mean_masked = mld_monthly_mean.where(mask == 1)
    mld_max_masked  = mld_monthly_max.where(mask == 1)

    # -----------------------------------------------------------------------
    # Compute spatial mean over the region
    # -----------------------------------------------------------------------
    mean_1D = mld_mean_masked.mean(['latitude', 'longitude'], skipna=True)
    max_1D  = mld_max_masked.mean(['latitude', 'longitude'], skipna=True)

    # Convert to DataFrame
    df_mean = mean_1D.to_dataframe().reset_index().rename(columns={'mlotst': 'mlotst_mean'})
    df_max  = max_1D.to_dataframe().reset_index().rename(columns={'mlotst': 'mlotst_max'})

    # Merge monthly mean and max
    df_year = df_mean.merge(df_max, on='month', how='left')
    df_year['year'] = yy
    df_year['date'] = pd.to_datetime(df_year[['year', 'month']].assign(day=1))

    df_year = df_year[['date', 'year', 'month', 'mlotst_mean', 'mlotst_max']]

    # Append to full dataset
    df_mean_max = pd.concat([df_mean_max, df_year], ignore_index=True)

###############################################################################
# 4. Save CSV
###############################################################################

output_csv = path_save + f"mlotst_monthly_mean_max_IR_LAB_{start_year}_{end_year}.csv"

df_mean_max.to_csv(
    output_csv,
    encoding="ISO-8859-1",
    sep=',',
    decimal='.',
    index=False
)

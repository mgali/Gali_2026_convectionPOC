###############################################################################
#  Compute Yearly Maximum MLD (1993–2019)
###############################################################################
# ---------------------------------------------------------------------------
# 1. Load daily MLD data
# 2. Compute yearly maximum MLD over the selected region
# 3. Save results as NetCDF
# ---------------------------------------------------------------------------
# Author    : Maria Sanchez Urrea
# Email     : mariasu@icm.csic.es
# Created   : November 2024
###############################################################################

import os
import glob
import xarray as xr
import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings('ignore')

###############################################################################
# 1. User settings
###############################################################################

# Years to include in the analysis
year_selection = list(np.arange(1993, 2019 + 1))

# Paths
path_data = 'D:/Arctic/projects/deep_convection_NA/data/daily_data/'
path_save = 'D:/Arctic/projects/deep_convection_NA/outputs/'

###############################################################################
# 2. Load Data
###############################################################################

# List all files in the directory
all_files = os.listdir(path_data)

# Select only files corresponding to the chosen years
files = [
    f for f in all_files
    if f.startswith('mlotst_') and any(f'_{year}' in f for year in year_selection)
]

# Add full path
files_with_path = [os.path.join(path_data, f) for f in files]

# Open all selected files as a single dataset
dset = xr.open_mfdataset(files_with_path, chunks='auto')
mld = dset.mlotst
dset.close()

###############################################################################
# 3. Compute Yearly Maximum MLD
###############################################################################

# Group by year and compute the maximum value for each year
mld_yearly_max = mld.groupby('time.year').max('time', keep_attrs=True)

###############################################################################
# 4. Save Output
###############################################################################

output_file = (
    path_save +
    f"mlotst_max_yearly_{year_selection[0]}_{year_selection[-1]}.nc"
)

mld_yearly_max.to_netcdf(output_file, 'w', format='NETCDF4_CLASSIC')

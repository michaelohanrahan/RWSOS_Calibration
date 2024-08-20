"""
This script is used for preparing observation data in a specific format (the same as the obs dataset used in Meuse model)
The script is created by Jing Deng.
Date: 2024-08-19
"""

import xarray as xr
import numpy as np
from pathlib import Path

work_dir = Path(r'c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\rhine\data\1-external')
fn_obs = work_dir / 'discharge_obs_hr_appended.nc'  # original obs file

# load original obs
obs = xr.open_dataset(fn_obs)
obs = obs.assign_coords(stations=('stations', obs.stations.values))

# create standard obs nc file format (reference: Meuse)
variables = ['Q', 'x', 'y', 'z', 'lon', 'lat', 'station_id', 'station_names']
time_rng = obs.time.values
wflow_id = obs.stations.values
runs = ['Obs.']

# Initialize arrays
S1 = np.full((len(time_rng), len(wflow_id), len(runs)), np.nan, dtype=np.float64)  # Float array for 'Q'
S2_float = np.full((len(wflow_id), len(runs)), np.nan, dtype=np.float64)  # Float array for other numeric variables
S2_str_id = np.full((len(wflow_id), len(runs)), '', dtype='|S64')  # Byte string array for 'station_id'
S2_str_names = np.full((len(wflow_id), len(runs)), '', dtype='|S255')  # Byte string array for 'station_names'

# Create the empty Dataset with the correct types
h = {}
for var in variables:
    if var == 'Q':
        h[var] = (('time', 'wflow_id', 'runs'), S1)
    elif var == 'station_id':
        h[var] = (('wflow_id', 'runs'), S2_str_id)
    elif var == 'station_names':
        h[var] = (('wflow_id', 'runs'), S2_str_names)
    else:
        h[var] = (('wflow_id', 'runs'), S2_float)

ds = xr.Dataset(
    data_vars=h,
    coords={'time': time_rng,
            'wflow_id': wflow_id,
            'runs': runs})

# copy data from original obs to new ds
for var in variables:
    if var == 'Q':
        ds[var].loc[:, :, 'Obs.'] = obs['Qm'].loc[:, :].values
    else:
        ds[var].loc[:, 'Obs.'] = obs[var].loc[:].values


# check if obs Qm are correctly copied to ds Q
ds_values = ds.sel(runs='Obs.', wflow_id=709).Q.values
obs_values = obs.sel(stations=709).Qm.values
# Compare for NaNs separately
nan_indices_ds = np.isnan(ds_values)
nan_indices_obs = np.isnan(obs_values)
nan_equal = np.all(nan_indices_ds == nan_indices_obs)  # Check if NaNs are in the same positions in both arrays
# Check if all non-NaN values are equal
values_equal = np.all(ds_values[nan_indices_ds == False] == obs_values[nan_indices_obs == False])
# Combine the checks
overall_equal = nan_equal and values_equal
print(f'Are arrays equal considering NaNs? {overall_equal}')


# save to file
fn_ds = work_dir / 'discharge_obs_hr_FORMAT_allvars.nc'
ds.to_netcdf(fn_ds)
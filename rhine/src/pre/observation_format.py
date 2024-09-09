"""
This script is used for preparing observation data in a specific format (the same as the obs dataset used in Meuse model)
The script is created by Jing Deng.
Date: 2024-08-19

Important modification on Sep 4, 2024:
In the original observation data, the stations (id) start from 0, which is not allowed in wflow model to be used as wflow_id. 
Therefore, the id of the stations 0 (b'Livange', [b'HLUX_Livange', [6.1149567, 49.52647733]) was manually changed to 727.
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

# Initialize the data dictionary
h = {}
# Create 3D array for 'Q'
h['Q'] = (('time', 'wflow_id', 'runs'), np.full((len(time_rng), len(wflow_id), len(runs)), np.nan, dtype=np.float64))
# Loop to create 2D arrays for numeric variables
for var in ['x', 'y', 'z', 'lon', 'lat']:
    h[var] = (('wflow_id', 'runs'), np.full((len(wflow_id), len(runs)), np.nan, dtype=np.float64))
# Initialize arrays for string variables with appropriate data types
h['station_id'] = (('wflow_id', 'runs'), np.full((len(wflow_id), len(runs)), '', dtype='|S64'))
h['station_names'] = (('wflow_id', 'runs'), np.full((len(wflow_id), len(runs)), '', dtype='|S255'))


ds = xr.Dataset(
    data_vars=h,
    coords={'time': time_rng,
            'wflow_id': wflow_id,
            'runs': runs})

# copy data from original obs to new ds
for var in variables:
    if var == 'Q':
        ds[var].loc[:, :, 'Obs.'] = obs['Qm'].values
    else:
        ds[var].loc[:, 'Obs.'] = obs[var].values


# check if obs vars are correctly copied to ds vars
for var in variables:
    if var == 'Q':
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

    else:
        ds_values = ds.sel(runs='Obs.', wflow_id=709)[var].values
        obs_values = obs.sel(stations=709)[var].values
        values_equal = np.all(ds_values == obs_values)
        overall_equal = values_equal
        print(f'Are arrays equal? {overall_equal}')



# save to file
fn_ds = work_dir / 'discharge_obs_hr_FORMAT_allvars.nc'
ds.to_netcdf(fn_ds)

# modify the wflow_id of station 0 to 727
fn_ds = work_dir / 'discharge_obs_hr_FORMAT_allvars.nc'
ds = xr.open_dataset(fn_ds)

# Change wflow_id 0 to 727
ds['wflow_id'] = ds['wflow_id'].where(ds['wflow_id'] != 0, 727)
# Save the modified dataset
ds.to_netcdf(work_dir / 'discharge_obs_hr_FORMAT_allvars_wflowid_0_to_727.nc')




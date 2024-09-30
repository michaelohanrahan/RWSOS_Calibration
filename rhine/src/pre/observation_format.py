"""
This script is used for preparing observation data in a specific format (the same as the obs dataset used in Meuse model)
The script is created by Jing Deng.
Date: 2024-08-19

Important modification on Sep 4, 2024:
In the original observation data, the stations (id) start from 0, which is not allowed in wflow model to be used as wflow_id. 
Therefore, the id of the stations 0 (b'Livange', [b'HLUX_Livange', [6.1149567, 49.52647733]) was manually changed to 727.

Important modification on Sep 30, 2024:
New data for Basel and Maxau stations

"""

import xarray as xr
import numpy as np
from pathlib import Path
import pandas as pd

#%%
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


#%% New data for Basel and Maxau stations
# Function to detect file encoding
import chardet
def detect_encoding(file_path):
    with open(file_path, 'rb') as file:
        raw_data = file.read()
    return chardet.detect(raw_data)['encoding']

# Detect the file encoding
def load_and_preprocess_csv(file_path):
    detected_encoding = detect_encoding(file_path)
    df = pd.read_csv(file_path, 
                     sep=';',  # Semicolon separator
                     decimal=',',  # Comma as decimal separator
                     encoding=detected_encoding, # UTF-8 encoding for special characters
                     skiprows=14)

    df['time'] = pd.to_datetime(df['Datum'] + ' ' + df['Uhrzeit'], 
                                     format='%d.%m.%Y %H:%M:%S', 
                                     errors='coerce')
    
    df['time'] = df['time'].astype('datetime64[ns]')
    df.set_index('time', inplace=True)

    df = df.rename(columns={'Durchfluss [m³/s]': 'Q'})
    df['Q'] = df['Q'].astype(float)

    df = df.drop(columns=['Datum', 'Uhrzeit'])
    return df

# update Q values for new data
def update_q_values(obs, df, wflow_id):
    # Convert df to an xarray DataArray
    da = xr.DataArray(df['Q'].values, coords={'time': df.index}, dims=['time'])

    # Create a boolean mask for the rows in obs where wflow_id matches
    mask = obs['wflow_id'] == wflow_id

    # Reindex da to match the time coordinates of obs
    da_reindexed = da.reindex(time=obs.sel(runs='Obs.')['time'], method='nearest')

    # Update Q values based on da
    new_q = obs.sel(runs='Obs.')['Q'].copy(deep=True)
    new_q = xr.where(mask, da_reindexed, new_q)

    # Set Q to NaN for all times not in da
    time_in_da = obs.sel(runs='Obs.')['time'].isin(da.time)
    new_q = xr.where(mask & ~time_in_da, np.nan, new_q)

    # Update the original obs Dataset
    obs_updated = obs.copy(deep=True)
    obs_updated['Q'].loc[{'runs': 'Obs.'}] = new_q
    return obs_updated


work_dir = Path(r'c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\rhine\data\1-external')
fn_obs = work_dir / 'discharge_obs_hr_FORMAT_allvars_wflowid_0_to_727.nc'
new_maxau = work_dir / 'new_data_bafg/Rhein-Q60/23700200-Rhein-Maxau-Q60.csv'
new_basel = work_dir / 'new_data_bafg/Rhein-Q60/2310010-Rhein-Basel,Rheinhalle-Q60.csv'

# load original obs
obs = xr.open_dataset(fn_obs)

# load new data
df_maxau = load_and_preprocess_csv(new_maxau)
df_basel = load_and_preprocess_csv(new_basel)

# update Q values for new data
obs_updated = update_q_values(obs, df_maxau, 688)
obs_updated = update_q_values(obs_updated, df_basel, 705)

# save to file
fn_ds = work_dir / 'discharge_obs_hr_FORMAT_allvars_wflowid_0_to_727_with_new_data.nc'
obs_updated.to_netcdf(fn_ds)



# # compare old and new data (maxau_688, basel_705)
# import plotly.graph_objects as go
# # Maxau
# old_da = obs.sel(wflow_id=688, time=slice('2005-01-01', '2016-12-31'))
# old_df = old_da.to_dataframe().reset_index()
# old_df['time'] = pd.to_datetime(old_df['time'])

# df_maxau.index = pd.to_datetime(df_maxau.index)

# fig = go.Figure()
# fig.add_trace(go.Scatter(x=old_df['time'], y=old_df['Q'], mode='lines', name='Old Data'))
# fig.add_trace(go.Scatter(x=df_maxau.index, y=df_maxau['Q'], mode='lines', name='New Data'))
# fig.update_layout(
#     title='Maxau_688: Old vs New Data',
#     xaxis_title='Time',
#     yaxis_title='Discharge (m3/s)',
#     legend_title='Data Source'
# )
# fig.write_html(work_dir / 'new_data_bafg/Rhein-Q60' / 'compare_Maxau_688.html')

# # Basel
# old_da = obs.sel(wflow_id=705, time=slice('2005-01-01', '2016-12-31'))
# old_df = old_da.to_dataframe().reset_index()
# old_df['time'] = pd.to_datetime(old_df['time'])

# df_basel.index = pd.to_datetime(df_basel.index)

# fig = go.Figure()
# fig.add_trace(go.Scatter(x=old_df['time'], y=old_df['Q'], mode='lines', name='Old Data'))
# fig.add_trace(go.Scatter(x=df_basel.index, y=df_basel['Q'], mode='lines', name='New Data'))
# fig.update_layout(
#     title='Basel_705: Old vs New Data',
#     xaxis_title='Time',
#     yaxis_title='Discharge (m3/s)',
#     legend_title='Data Source'
# )
# fig.write_html(work_dir / 'new_data_bafg/Rhein-Q60' / 'compare_Basel_705.html')
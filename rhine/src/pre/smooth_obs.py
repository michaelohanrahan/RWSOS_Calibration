import xarray as xr
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import numpy as np
from pathlib import Path

root_dir = Path(r'c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\rhine')
observations_path = root_dir / 'data/1-external/discharge_obs_hr_FORMAT_allvars_wflowid_0_to_727.nc'
ds = xr.open_dataset(observations_path)
time_range=('1996', None)

sigma=5
wflow_ids = [425, 576, 539, 704, 423, 476]

da = ds.sel(wflow_id=425, runs='Obs.', time=slice(*time_range)).Q
count1 = np.count_nonzero(~np.isnan(da.values))
da_smooth = xr.DataArray(
    gaussian_filter1d(da.values, sigma=sigma, mode='nearest'),
    dims=da.dims,
    coords=da.coords,
    name='Q',
    attrs=da.attrs
)
count2 = np.count_nonzero(~np.isnan(da_smooth.values))
print(f"{'*'*10} {425} {'*'*10}")
print(f"data points before: {count1}, after: {count2}")
print(f"data matches time dimension: {len(da.time.values) == len(da_smooth.time.values)}")

# Plotting da and da_smooth
import plotly.graph_objects as go

fig = go.Figure()
fig.add_trace(go.Scatter(x=da.time, y=da, mode='lines', name='Original Data'))
fig.add_trace(go.Scatter(x=da_smooth.time, y=da_smooth, mode='lines', name='Smoothed Data'))
fig.update_layout(title='Original vs Smoothed Data', xaxis_title='Time', yaxis_title='Discharge (m3/s)')
fig.show()







# t0='2017-06-25'
# t1='2017-07-30'
sigma=5
das = []
for wflow_id in ds.wflow_id.values:
    da = ds.sel(
        # time=slice(t0,t1),
        wflow_id=wflow_id, 
        runs='Obs.'
        ).Q
    #count non nan values
    count1 = np.count_nonzero(~np.isnan(da.values))
    da_smooth = gaussian_filter1d(da, sigma=sigma, mode='nearest')
    count2 = np.count_nonzero(~np.isnan(da_smooth))
    print(f"{'*'*10} {wflow_id} {'*'*10}")
    print(f"data points before: {count1}, after: {count2}")
    # assert count1 == count2, "Data points before and after smoothing do not match"
    print(f"data matches time dimension: {len(da.time.values) == len(da_smooth)}")
    
    da_smooth = xr.DataArray(
        np.float64(da_smooth),
        dims=da.dims,
        coords=da.coords,
        name='Q',
        attrs=da.attrs
        )
    # print(da_smooth)
    das.append(da_smooth)
    # a
ds_smooth = xr.Dataset({'Q':xr.concat(das, dim='wflow_id')})
print(ds_smooth)

ds_combined = xr.concat([ds.sel(runs='HBV'), ds_smooth], dim='runs')
print(ds_combined)

ds_combined.to_netcdf('discharge_hourlyobs_smoothed.nc')
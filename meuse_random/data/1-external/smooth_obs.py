import xarray as xr
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import numpy as np
import os 

os.chdir(r'p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse_random\data\1-external')
print(os.getcwd())
ds = xr.open_dataset('discharge_hourlyobs_HBV_combined.nc')
print(ds)

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
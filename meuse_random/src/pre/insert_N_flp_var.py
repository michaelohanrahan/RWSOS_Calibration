import xarray as xr 
import numpy as np
import matplotlib.pyplot as plt

# Open the dataset
ds = xr.open_dataset(r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\3-input\staticmaps\staticmaps_old.nc")

# Create a constant array with the same shape and dimensions as 'wflow_dem'
constant_value = 0.072
var = np.full_like(ds['wflow_dem'], constant_value)

# Apply the NaN mask from 'wflow_dem' to the constant array
var = np.where(np.isnan(ds['wflow_dem']), np.nan, var)

# Create a DataArray with the same dims, coords, and attrs as 'wflow_dem'
da = xr.DataArray(var, dims=ds['wflow_dem'].dims, coords=ds['wflow_dem'].coords, attrs=ds['wflow_dem'].attrs)

# now we insert N_Floodplain into the dataset
ds = ds.assign(N_Floodplain=da)
ds['N_Floodplain'].attrs = {'long_name': 'N_Floodplain', 'units': '-'}
ds.to_netcdf(r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\3-input\staticmaps\staticmaps.nc")
import xarray as xr
import geopandas as gpd
import os 

os.chdir(r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\1-external")

ex = xr.open_dataset('discharge_obs_combined_EXAMPLE.nc')

'''
xarray.Dataset
Dimensions:
    Q_gauges_obs:   67    
    time:           105937

Coordinates:
    Q_gauges_obs    (Q_gauges_obs)  <U8             '12105900' ... '12089500'
    time            (time)          datetime64[ns]  1986-10-01T06:00:00 ... 2023-01-...
    
Data variables:
    Q               (time, Q_gauges_obs)  float64
    
Indexes:
    Q_gauges_obs    PandasIndex
    time            PandasIndex
Attributes: (0)

'''
ds = xr.open_dataset('ds_obs_model_combined.nc')

ds.sel(runs='Obs.')

for id in ds.wflow_id.values:
    #check what proportion of the whole time series is missing in variable Q
    # print(ds)
    HBVnull = ds.sel(wflow_id=id, runs='HBV').Q.isnull().sum().values/ds.sel(wflow_id=id, runs='Obs.').Q.size
    OBSnull = ds.sel(wflow_id=id, runs='Obs.').Q.isnull().sum().values/ds.sel(wflow_id=id, runs='Obs.').Q.size
    
    print(id, HBVnull, OBSnull)
    
    if HBVnull < 0.75:
        print(f'keep {id}')
        with open('HBV_in_ds', 'a') as f:
            f.write(f'{id}\n')
            
    if OBSnull < 0.75:
        print(f'keep {id}')
        with open('OBS_in_ds', 'a') as f:
            f.write(f'{id}\n')

ds = ds.sel(runs=['Obs.', 'HBV'])

ds.to_netcdf('discharge_obs_HBV_combined.nc')

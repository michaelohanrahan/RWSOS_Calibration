import xarray as xr
import geopandas as gpd
import pandas as pd 
import os 

os.chdir(r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\1-external")

ex = xr.open_dataset('discharge_obs_combined_EXAMPLE.nc')

def clip_nan(da):
    # Convert to pandas Series
    s = da.to_series()
    
    # Remove NaNs from the beginning and end
    s = s.loc[s.first_valid_index():s.last_valid_index()]
    
    #keep id as coordinate
    r_da = xr.DataArray(s)
    r_da.name = da.name
    
    #create two new indexes
    
    r_da.coords['wflow_id'] = 9
    r_da = r_da.expand_dims('wflow_id')
    r_da.coords['runs'] = 'Obs.'
    r_da = r_da.expand_dims('runs')
    r_da.coords['time'] = r_da.indexes['time']
    
    r_da = r_da.set_index(wflow_id='wflow_id', runs='runs', time='time')
    
    return r_da
    
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

ds2 = xr.open_dataset(r"P:\archivedprojects\11208719-interreg\data\spw\Discharge\c_final\hourly_spw.nc")
Sal = ds2.sel(id=7319) #Salzinnes - Ronet
Sal = clip_nan(Sal.Q)
Sal_ds = xr.Dataset({'Q': Sal})

ds = ds.sel(runs=['Obs.', 'HBV'])
ds = xr.merge([ds, Sal], join='outer')

for id in ds.wflow_id.values:
    #check what proportion of the whole time series is missing in variable Q
    # print(ds)
    HBVnull = ds.sel(wflow_id=id, runs='HBV').Q.isnull().sum().values/ds.sel(wflow_id=id, runs='Obs.').Q.size
    OBSnull = ds.sel(wflow_id=id, runs='Obs.').Q.isnull().sum().values/ds.sel(wflow_id=id, runs='Obs.').Q.size
    
    print(id, HBVnull, OBSnull)
    
    if HBVnull < 0.75:
        print(f'keep {id}')
        with open('HBV_in_ds.txt', 'a') as f:
            f.write(f'{id}\n')
            
    if OBSnull < 0.75:
        print(f'keep {id}')
        with open('OBS_in_ds.txt', 'a') as f:
            f.write(f'{id}\n')



ds.to_netcdf('discharge_obs_HBV_combined.nc')

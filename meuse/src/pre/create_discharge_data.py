import xarray as xr
import geopandas as gpd
import pandas as pd 
import os 
from hydromt import DataCatalog as DC
from argparse import ArgumentParser as AP
import traceback
import sys
import numpy as np
import dask
from shapely.geometry import Point
import matplotlib.pyplot as plt
from icecream import ic
import shutil
import numpy as np
from dask.diagnostics import ProgressBar

def cleanup(id, da):
    # Convert to pandas Series
    s = da.to_series()
    
    # Remove NaNs from the beginning and end
    s = s.loc[s.first_valid_index():s.last_valid_index()]
    
    #keep id as coordinate
    r_da = xr.DataArray(s)
    r_da.name = da.name
    
    #create two new indexes
    r_da.coords['wflow_id'] = id
    r_da = r_da.expand_dims('wflow_id')
    r_da.coords['runs'] = 'Obs.'
    r_da = r_da.expand_dims('runs')
    r_da.coords['time'] = r_da.indexes['time']
    
    r_da = r_da.set_index(wflow_id='wflow_id', runs='runs', time='time')
    
    return r_da


if __name__ == "__main__":
    os.chdir(r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\1-external")
    #open dataset with example structure
    if os.path.exists('discharge_obs_combined.nc') and not os.path.exists('discharge_obs_combined_BACKUP.nc'):
        shutil.move('discharge_obs_combined.nc', 'discharge_obs_combined_BACKUP.nc')
    
    df = xr.open_dataset('discharge_obs_HBV_combined.nc')
    
    min_date = df.time.min().values
    max_date = df.time.max().values
    
    df_ids = df.wflow_id.values
    id_ref = gpd.read_file('../2-interim/QGIS/to_wflow/hourly_gauges_all.geojson', fid='wflow_id')
   
    #the ids chosen not to include in the gaugemap
    with open('../2-interim/ignore_list.txt', 'r') as f:
        ignore_list = [int(line.strip()) for line in f if line.strip().isdigit()]
    
    #concat the ids with ignores 
    ignore = np.concatenate([ignore_list, df_ids])
    
    #remove the ignore ids from the id_ref
    id_ref = id_ref.loc[~id_ref.wflow_id.isin(ignore)]
    
    #dont use wflow_id for these values because they are the old ones
    wanted_ids = id_ref.id.values
    
    #load the example dataset
    datacatalog = DC(f'observed_discharge/data_meuse.yml')
    flong = 'hourly'
    
    #Find the sources that make hourly obs data
    sources = [f for f in datacatalog.keys if flong in f and 'stats' not in f]
    
    #these are the datasets that contain any of the wanted ids
    datasets_to_integrate = {}
    
    for source in sources:
        ds = datacatalog.get_dataset(source)
        
        if isinstance(ds, xr.Dataset):
            ds = ds.Q
        
        ds = ds.sel(time=slice(min_date, max_date))
        
        if np.any([id in wanted_ids for id in ds.wflow_id.values]):
            # print(ds)
            ids = np.intersect1d(wanted_ids, ds.wflow_id.values)
            # print(ds, '\n', ids)
            for id in ids:
                datasets_to_integrate[id] = ds.sel(wflow_id=id)
    
    for id, da in datasets_to_integrate.items():
        #insert into the df
        ds = xr.Dataset({'Q':cleanup(id, da)})
        df = xr.merge([df, ds], join='outer')
    
    if np.all(id_ref.wflow_id.isin(df.wflow_id.values)):
        print('All ids are in the dataset')
        df = df.chunk({'wflow_id': 1, 'time':-1, 'runs':1})
        delayed = df.to_netcdf('discharge_hourlyobs_HBV_combined.nc', compute=False)
        with ProgressBar():
            delayed.compute()
    
    
    
    
    
# ds = xr.open_dataset('ds_obs_model_combined.nc')

# ds2 = xr.open_dataset(r"P:\archivedprojects\11208719-interreg\data\spw\Discharge\c_final\hourly_spw.nc")
# Sal = ds2.sel(id=7319) #Salzinnes - Ronet
# Sal = clip_nan(Sal.Q)
# Sal_ds = xr.Dataset({'Q': Sal})

# ds = ds.sel(runs=['Obs.', 'HBV'])
# ds = xr.merge([ds, Sal], join='outer')

# for id in ds.wflow_id.values:
#     #check what proportion of the whole time series is missing in variable Q
#     # print(ds)
#     HBVnull = ds.sel(wflow_id=id, runs='HBV').Q.isnull().sum().values/ds.sel(wflow_id=id, runs='Obs.').Q.size
#     OBSnull = ds.sel(wflow_id=id, runs='Obs.').Q.isnull().sum().values/ds.sel(wflow_id=id, runs='Obs.').Q.size
    
#     print(id, HBVnull, OBSnull)
    
#     if HBVnull < 0.75:
#         print(f'keep {id}')
#         with open('HBV_in_ds.txt', 'a') as f:
#             f.write(f'{id}\n')
            
#     if OBSnull < 0.75:
#         print(f'keep {id}')
#         with open('OBS_in_ds.txt', 'a') as f:
#             f.write(f'{id}\n')



# ds.to_netcdf('discharge_obs_HBV_combined.nc')

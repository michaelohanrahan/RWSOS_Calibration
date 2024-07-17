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

def health_check(ds: xr.Dataset, health_check_path: str):
    if not isinstance(ds, xr.Dataset):
        ds = xr.Dataset(data_vars={'Q': ds})
        print(f'\n {"="*10} dataset converted to xr.Dataset: {ds}')
    
    for id in ds.wflow_id.values:
        q_data = ds.sel(wflow_id=id).Q
        nans = np.isnan(q_data).sum().compute()  # Compute the sum of NaNs
        non_nan_indices = np.where(~np.isnan(q_data.compute()))[0]  # Compute q_data before finding non-NaN indices
        
        if non_nan_indices.size > 0:
            min_nonnan_date = q_data.time[non_nan_indices[0]].values
            max_nonnan_date = q_data.time[non_nan_indices[-1]].values
        else:
            min_nonnan_date = 'all NaN'
            max_nonnan_date = 'all NaN'
        
        zeros = (q_data == 0).sum().compute()  # Compute the sum of zeros
        
        try:
            station_name = ds.sel(wflow_id=id).station_name.data.compute()  # Attempt to compute station_name
        except:
            station_name = 'name not found'
        
        with open(health_check_path, 'a') as f:
            f.write(f'WF_id: {id}, Station: {station_name}, NaN: {nans}, first valid date: {min_nonnan_date}, last valid date: {max_nonnan_date}, n zeros: {zeros}\n')

def get_dc_data(ls:list, model:str, cwd:str, plat:str, freq:str='H'):
    '''
    The main function is passed the argument of a list of wflow ids to combine into a subcatchment set
    args:
        ls: list of wflow ids to combine
        model: the model name
        cwd: the current working directory
        plat: the platform to use
    creates:
        a netcdf file with the combined subcatchment data
        
    returns:

    '''
    os.chdir(cwd)
    print(f'Working in {os.getcwd()}')
    
    #load the example dataset
    datacatalog = DC(f'observed_discharge/data_{model}{plat}.yml')
    
    if freq == 'H':
        flong = 'hourly'
    elif freq == 'D':
        flong = 'daily'
        
    #Find the sources that make hourly obs data
    sources = [f for f in datacatalog.keys if flong in f and 'stats' not in f]
    print(f'{flong} keys: {sources}')
    
    health_check_path = os.path.join(cwd, f'{flong}_health_check.txt')
    
    if os.path.exists(health_check_path):
        os.remove(health_check_path)
    
    for key in hourly:
        ds = datacatalog.get_dataset(key)
        src = datacatalog.get_source(key)
        print(f'\n {"*"*20} \nWorking on {key}')
        print(f'dataset: {ds}')
        health_check(ds, health_check_path)  # Updated variable name here
        print(f'Finished {key}\n {"*"*20}')
  
if __name__ == '__main__':
    parser = AP()
    parser.add_argument('--wflow_ids', type=list, default=[9])
    parser.add_argument('--model', type=str, default='meuse')
    parser.add_argument('--cwd', type=str, default=r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\1-external")
    parser.add_argument('--freq', type=str, default='H')
    args = parser.parse_args()
    
    if sys.platform == 'win32':
        plat=''
    else:
        plat='_linux'
    try:
        print('')
        get_hourly_dc(args.wflow_ids, args.model, args.cwd, plat, args.freq)
    except Exception as e:
        print(f'Error: {e}')
        traceback.print_exc()
    

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

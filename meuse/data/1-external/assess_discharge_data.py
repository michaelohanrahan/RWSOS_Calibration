import xarray as xr
import geopandas as gpd
import pandas as pd 
import os 
from hydromt import DataCatalog as DC
from argparse import ArgumentParser as AP
import traceback
import sys
import numpy as np
<<<<<<< HEAD:meuse/data/1-external/assess_discharge_data.py
import geopandas as gpd
=======
import dask
>>>>>>> 13294849e6b5fc61d23333e82c83108544645b68:meuse/data/1-external/create_discharge_data_v2.py
from shapely.geometry import Point

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
<<<<<<< HEAD:meuse/data/1-external/assess_discharge_data.py
    
'''
EXAMPLE DATASET COMPATIBLE WITH ORIGINAL WORKFLOW

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

KEY DIFFERENCES TO EMPLOY ARE THE TWO PREVIOUS MODEL RUNS... Those will be HBV and FL1D in addition to Q
'''
=======

>>>>>>> 13294849e6b5fc61d23333e82c83108544645b68:meuse/data/1-external/create_discharge_data_v2.py
def health_check(ds: xr.Dataset, health_check_path: str):
    valid_stations = []  # Step 1: Initialize an empty list

    if not isinstance(ds, xr.Dataset):
        ds = xr.Dataset(data_vars={'Q': ds})
        
    for id in ds.wflow_id.values:
        q_data = ds.sel(wflow_id=id).Q
        nans = int(np.isnan(q_data).sum().compute())
        non_nan_indices = np.where(~np.isnan(q_data.compute()))[0]
        
        if non_nan_indices.size > 0:
            min_nonnan_date = q_data.time[non_nan_indices[0]].values
            max_nonnan_date = q_data.time[non_nan_indices[-1]].values
            valid_count = len(non_nan_indices)
            time_diffs = np.diff(q_data.time[non_nan_indices].values.astype('datetime64[h]'))
            most_common_diff = pd.Series(time_diffs).mode()[0]
            frequency = f"{most_common_diff / np.timedelta64(1, 'h')} hours"
        else:
            min_nonnan_date = 'all NaN'
            max_nonnan_date = 'all NaN'
            valid_count = 0
            frequency = 'N/A'
        
        zeros = int((q_data == 0).sum().compute())
        
        try:
            station_name = ds.sel(wflow_id=id).station_name.data.compute()
        except:
            try:
                station_name = ds.sel(wflow_id=id)['Libellé'].data.compute()
            except:
                station_name = 'name not found'
        
        nan_proportion = nans / valid_count if valid_count > 0 else 0.
        
        if nan_proportion <= 0.25 and non_nan_indices.size > 0:  # Step 2: Check if proportion 
            valid_stations.append({
                'id': int(id),
                'name': str(station_name),
                'invalid_proportion': float(nan_proportion),
                'valid_count': int(valid_count),
                'min_date': str(min_nonnan_date),
                'max_date': str(max_nonnan_date),
                'frequency': str(frequency),
                'zeros': int(zeros),
                'x': float(ds.sel(wflow_id=id).x.values),
                'y': float(ds.sel(wflow_id=id).y.values)
            })
        
            with open(health_check_path, 'a') as f:
                f.write(f'{"."*50}\nWF_id: {id}, len valid: {valid_count}, proportion NaN: {nan_proportion:.3f}, frequency: {frequency}, Station: {station_name}, \n\
                        first valid date: {min_nonnan_date}, last valid date: {max_nonnan_date}, \n\
                        n zeros: {zeros}\n {"."*50}\n')
    return valid_stations
    
<<<<<<< HEAD:meuse/data/1-external/assess_discharge_data.py
  
=======
>>>>>>> 13294849e6b5fc61d23333e82c83108544645b68:meuse/data/1-external/create_discharge_data_v2.py
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
<<<<<<< HEAD:meuse/data/1-external/assess_discharge_data.py
        
    all_stations = []
=======
    
    all_stations = []  # Initialize an empty list to store all valid stations
    
>>>>>>> 13294849e6b5fc61d23333e82c83108544645b68:meuse/data/1-external/create_discharge_data_v2.py
    for key in sources:
        ds = datacatalog.get_dataset(key)
        src = datacatalog.get_source(key)
        with open(health_check_path, 'a') as f:
            f.write(f'\n {"="*100} \nWorking on {key}\n')
            f.write(f'dataset: {ds}\n')
            
        print(f'\n {"*"*20} \nWorking on {key}')
        print(f'dataset: {ds}')
        valid_stations = health_check(ds, health_check_path)  # Updated variable name here
<<<<<<< HEAD:meuse/data/1-external/assess_discharge_data.py
=======
        #create geojson 
>>>>>>> 13294849e6b5fc61d23333e82c83108544645b68:meuse/data/1-external/create_discharge_data_v2.py
        all_stations.extend(valid_stations)
        print(f'Finished {key}\n {"*"*20}')
    # Step 3: Convert the list of dictionaries into a DataFrame
    df_valid_stations = pd.DataFrame(all_stations)
    print(df_valid_stations)
    # Step 4: Convert DataFrame to GeoDataFrame and save as GeoJSON
    gdf = gpd.GeoDataFrame(df_valid_stations, geometry=[Point(xy) for xy in zip(df_valid_stations['x'], df_valid_stations['y'])])
    gdf.to_file("valid_stations.geojson", driver='GeoJSON')
    
    # Step 3: Convert the list of dictionaries into a DataFrame
    df_valid_stations = pd.DataFrame(all_stations)
    print(df_valid_stations)
    # Step 4: Convert DataFrame to GeoDataFrame and save as GeoJSON
    gdf = gpd.GeoDataFrame(df_valid_stations, geometry=[Point(xy) for xy in zip(df_valid_stations['x'], df_valid_stations['y'])])
    gdf.to_file(f"{flong}_valid_stations.geojson", driver='GeoJSON')
    
  
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
        get_dc_data(args.wflow_ids, args.model, args.cwd, plat, args.freq)
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

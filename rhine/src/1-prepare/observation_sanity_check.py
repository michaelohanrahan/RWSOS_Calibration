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

def health_check(ds: xr.Dataset, health_check_path: str, flong):
    valid_stations = []  # Step 1: Initialize an empty list

    if not isinstance(ds, xr.Dataset):
        ds = xr.Dataset(data_vars={'Q': ds})
        
    for id in ds.wflow_id.values:
        q_data = ds.sel(wflow_id=id).Q
        non_nan_indices = np.where(~np.isnan(q_data.compute()))[0]
        
        # Filter out zeros from non-NaN indices
        non_zero_non_nan_indices = non_nan_indices[q_data[non_nan_indices] != 0]
        
        if non_zero_non_nan_indices.size > 0:
            min_nonnan_date = q_data.time[non_zero_non_nan_indices[0]].values
            min_doi = pd.to_datetime(np.datetime64('2005-01-01T00:00:00.000000000')).strftime('%Y-%m-%d %H:%M:%S')
            max_nonnan_date = q_data.time[non_zero_non_nan_indices[-1]].values
            nans = int(np.isnan(q_data.sel(time=slice(min_doi, max_nonnan_date))).sum().compute())
            valid_count = len(non_zero_non_nan_indices)
            time_diffs = np.diff(q_data.time[non_zero_non_nan_indices].values.astype('datetime64[h]'))
            most_common_diff = pd.Series(time_diffs).mode()[0]
            frequency = f"{most_common_diff / np.timedelta64(1, 'h')} hours"
            total_len = int(len(pd.date_range(min_doi, max_nonnan_date, freq='H')))
        else:
            min_nonnan_date = 'all NaN or zero'
            max_nonnan_date = 'all NaN or zero'
            valid_count = 0
            frequency = 'N/A'
            total_len = 0
            
        
        zeros = int((q_data == 0).sum().compute())
        
        try:
            station_name = ds.sel(wflow_id=id).station_name.data.compute()
            if ',' in station_name:
                station_name = station_name.split(',')[0]
        except:
            try:
                station_name = ds.sel(wflow_id=id)['Libellé'].data.compute()
            except:
                station_name = None
        ic(station_name)
        # Sanitize station_name to remove or replace invalid characters for file names
        sanitized_station_name = str(station_name).replace('/', '_').replace('\\', '_')  # Example of replacing slashes with underscores

        plotted = False
        if not os.path.exists(f'observed_discharge/visualisation_{flong}/{station_name}.png') and station_name not in ['None', None] and os.path.exists(r"P:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\2-interim\QGIS\distance_filtered_10pct_keephighzeros_hourly_gauges.geojson"):
            try:
                plot(ds.sel(wflow_id=id).Q, sanitized_station_name, flong)
                plotted = True
            except TypeError as e:
                print(f'Error: {e}')
                traceback.print_exc()
                continue
                
        nan_proportion = nans / valid_count if valid_count > 0 else 0.
        
        # if nan_proportion <= 0.25 and non_nan_indices.size > 0:  # Step 2: Check if proportion 
        valid_stations.append({
            'id': int(id),
            'name': str(station_name),
            'invalid_proportion_2005-max': float(nan_proportion),
            'total_valid_2005-max': int(valid_count),
            'total_len_2005-max': total_len,
            'min_date': str(min_nonnan_date),
            'min_doi': str(min_doi),
            'max_date': str(max_nonnan_date),
            'frequency': str(frequency),
            'zeros': int(zeros),
            'x': float(ds.sel(wflow_id=id).x.values),
            'y': float(ds.sel(wflow_id=id).y.values),
            'plotted':plotted
        })
        
        with open(health_check_path, 'a') as f:
            f.write(f'{"."*50}\nWF_id: {id}, len valid: {valid_count}, proportion NaN: {nan_proportion:.3f}, frequency: {frequency}, Station: {station_name}, \n\
                    first valid date: {min_nonnan_date}, last valid date: {max_nonnan_date}, \n\
                    n zeros: {zeros}\n {"."*50}\n')
    return valid_stations

def plot(series, name, flong):
    fig, axs = plt.subplots(3, 1, figsize=(15, 10))
    plt.suptitle(name)
    series.plot(ax=axs[0])
    
    series.sel(time=slice('2005-01-01', '2005-12-31')).plot(ax=axs[1])
    axs[1].set_title('2005')
    series.sel(time=slice('2011-01-01', '2011-12-31')).plot(ax=axs[2])
    axs[2].set_title('2011')
    os.makedirs(f'observed_discharge/visualisation_{flong}', exist_ok=True)
    plt.savefig(f'observed_discharge/visualisation_{flong}/{name}.png', dpi=300)

def get_dc_data(model:str, cwd:str, plat:str, freq:str='H'):
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
    
    health_check_path = os.path.join(cwd, f'observed_discharge/{flong}_health_check.txt')
    
    if os.path.exists(health_check_path):
        os.remove(health_check_path)
    
    all_stations = []  # Initialize an empty list to store all valid stations
    
    for key in sources:
        ds = datacatalog.get_dataset(key)
        src = datacatalog.get_source(key)
        print(f'\n {"*"*20} \nWorking on {key}')
        print(f'dataset: {ds}')
        valid_stations = health_check(ds, health_check_path, flong)  # Updated variable name here
        #create geojson 
        all_stations.extend(valid_stations)
        print(f'Finished {key}\n {"*"*20}')
    
    # Step 3: Convert the list of dictionaries into a DataFrame
    df_valid_stations = pd.DataFrame(all_stations)
    print(df_valid_stations)
    # Step 4: Convert DataFrame to GeoDataFrame and save as GeoJSON
    gdf = gpd.GeoDataFrame(df_valid_stations, geometry=[Point(xy) for xy in zip(df_valid_stations['x'], df_valid_stations['y'])])
    gdf.to_file(f"observed_discharge/{flong}_stations.geojson", driver='GeoJSON')
    
  
if __name__ == '__main__':
    
    from pathlib import Path
    work_dir = Path(r'p:\11209265-grade2023\wflow\wflow_rhine_julia\measurements\vanBart')
    fn_obs = work_dir / 'discharge_obs_hr_appended.nc'
    
    # load original obs
    obs = xr.open_dataset(fn_obs)
    # create standard obs nc file format (reference: meuse: create_discharge_data.py, rhine:)
    variables = ['Q']
    S = np.zeros((len(rng), len(obs_keys), len(runs)))\
    v = (('time', 'wflow_id', 'runs'), S)
    h = {k:v for k in variables}
    
    ds = xr.Dataset(
        data_vars=h,
        coords={'time': rng,
                'stations': [station.split("_")[1] for station in obs_keys],
                'runs': runs})
    ds = ds * np.nan
    
    
    
    
    parser = AP()
    parser.add_argument('--model', type=str, default='meuse')
    parser.add_argument('--cwd', type=str, default=r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\1-external")
    parser.add_argument('--freq', type=str, default='H')
    args = parser.parse_args()
    
    if sys.platform == 'win32':
        plat=''
    else:
        plat='_linux'
        
    try:
        get_dc_data(model=args.model, 
                    cwd=args.cwd, 
                    plat=plat, 
                    freq=args.freq)
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



# ds.to_netcdf('discharge_obs_HBV_combined.nc'
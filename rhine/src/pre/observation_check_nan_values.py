"""
This script is borrowed from Meuse project, used for 
1) sanity check on observation data
2) create ignore_list.txt for gauges with high nan percentage

Created by Jing Deng
Date: 2024-08-19
"""

import xarray as xr
import geopandas as gpd
import pandas as pd 
import os
from pathlib import Path 
from hydromt import DataCatalog as DC
from argparse import ArgumentParser as AP
import traceback
import sys
import numpy as np
import dask
from shapely.geometry import Point
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from icecream import ic
import folium

def compute_nan_percentage_and_valid_period(obs: xr.Dataset,
                                            time_range: tuple = (None, None)):
    
    obs_temp = obs.sel(time=slice(*time_range))
    nan_percentage_list = []
    valid_start_list = []
    valid_end_list = []
    nan_percentage_within_valid_period_list = []

    for id in obs_temp.wflow_id.values:
        q_data = obs_temp.sel(wflow_id=id, runs='Obs.').Q
        
        # Compute nan percentage
        nan_percentage = q_data.isnull().values.mean() * 100
        nan_percentage_list.append(nan_percentage)
        
        # Find valid data period
        valid_data = q_data.dropna('time')
        if not valid_data.isnull().all():
            valid_start = valid_data.time[0].values
            valid_end = valid_data.time[-1].values
            valid_start_list.append(valid_start)
            valid_end_list.append(valid_end)
            nan_percentage_within_valid_period = q_data.sel(time=slice(valid_start, valid_end)).isnull().values.mean() * 100
            nan_percentage_within_valid_period_list.append(nan_percentage_within_valid_period)
        else:
            valid_start_list.append(None)
            valid_end_list.append(None)
            nan_percentage_within_valid_period_list.append(None)

    # Extract latitude and longitude for each wflow_id
    latitudes = obs_temp['lat'].sel(runs='Obs.').values
    longitudes = obs_temp['lon'].sel(runs='Obs.').values

    # Create a DataFrame
    data = pd.DataFrame({
        'wflow_id': obs_temp['wflow_id'].values,
        'nan_percentage': nan_percentage_list,
        'lat': latitudes,
        'lon': longitudes,
        'valid_start': valid_start_list,
        'valid_end': valid_end_list,
        'nan_percentage_within_valid_period': nan_percentage_within_valid_period_list
    })
    
    return data


def plot_gauge_map_nan_percentage(data: pd.DataFrame,
                                  work_dir: Path,
                                  time_range: tuple = (None, None)
                                  ):
    
    # Define the center of the Rhine River Basin
    map_center = [50.0, 8.0]  # Example center latitude and longitude

    # Define a color map and normalization for the legend
    cmap = plt.get_cmap('viridis')
    norm = mcolors.Normalize(vmin=min(data['nan_percentage']), vmax=max(data['nan_percentage']))

    # Create a Folium map object
    m = folium.Map(location=map_center, zoom_start=7, tiles='OpenStreetMap')

    # Add points to the map
    for _, row in data.iterrows():
        color = mcolors.to_hex(cmap(norm(row['nan_percentage'])))
        folium.CircleMarker(
            location=[row['lat'], row['lon']],
            radius=7,
            color=color,
            fill=True,
            fill_color=color,
            fill_opacity=0.7,
            popup=f"wflow_id: {row['wflow_id']}<br>Missing Percentage: {row['nan_percentage']:.2f}%"
        ).add_to(m)

    # Create a color legend with fixed size and location
    legend_html = '''
        <div style="position: fixed; 
                    bottom: 50px; left: 50px; width: 200px; height: 300px; 
                    border:2px solid grey; background-color:white; z-index:9999; 
                    font-size:14px; padding: 10px;">
        <b>Missing Percentage</b><br>
        <div style="height: 250px; width: 20px; background: linear-gradient(to top, 
                    #440154 0%, 
                    #482576 10%, 
                    #3f4a7d 20%, 
                    #31688e 30%, 
                    #26828e 40%, 
                    #1f9e89 50%, 
                    #35b779 60%, 
                    #6fef6c 70%, 
                    #aafc4f 80%, 
                    #d8e15f 90%, 
                    #f7f7f7 100%);
                    margin-right: 10px; display: inline-block; float: left;"></div>
        <div style="display: flex; flex-direction: column; margin-left: 10px;">
            <div>100%</div>
            <div>80%</div>
            <div>60%</div>
            <div>40%</div>
            <div>20%</div>
            <div>0%</div>
        </div>
        </div>
    '''

    # Save the map to an HTML file
    m.save(work_dir/f'missing_value_gauge_map_{time_range[0]}_{time_range[1]}.html')


def remove_gauge_from_obs(obs: xr.Dataset,
                          data: pd.DataFrame,
                          threshold: float
                          ):
    """Remove gauges with high nan percentage from observation dataset

    Args:
        obs (xr.Dataset): original observation dataset that contains all gauges
        data (pd.DataFrame): dataframe that contains nan percentage for each gauge
        threshold (float): if nan percentage is higher than threshold, remove the gauge. (in 100%)

    Returns:
        obs_filtered: observation dataset with gauges removed
    """
    
    high_nan_percentage_ids = data[data['nan_percentage'] > threshold]['wflow_id']
    print(f"Number of wflow_id with nan_percentage > {threshold}%: {len(high_nan_percentage_ids)}")
    print(f'Gauges to be removed: {high_nan_percentage_ids}')
    
    obs_temp = obs.copy()
    obs_filtered = obs_temp.where(~obs_temp.wflow_id.isin(high_nan_percentage_ids), drop=True)
    
    return obs_filtered


def high_nan_gauge_ignore_list(data: pd.DataFrame,
                               threshold: float,
                               wflow_id_exception: list, # gauge like Basel (705) has higher than 60% nan percentage, but we still want to keep it
                               work_dir: Path
                               ):
    
    
    high_nan_percentage_ids = data[(data['nan_percentage'] > threshold) & (~data['wflow_id'].isin(wflow_id_exception))]['wflow_id']
    print(f"Number of wflow_id with nan_percentage > {threshold}%: {len(high_nan_percentage_ids)}")
    print(f'Gauges to be added to ignore_list.txt: {high_nan_percentage_ids}')   
    
    with open(work_dir / 'ignore_list.txt', 'w') as f:
        for id in high_nan_percentage_ids:
            f.write(f"{id}\n")
    
    
if __name__ == '__main__':
    
    work_dir = Path(r'c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\rhine\data\1-external')
    obs = xr.open_dataset(work_dir / 'discharge_obs_hr_FORMAT_allvars_wflowid_0_to_727.nc')
    
    # 1996-2016
    data = compute_nan_percentage_and_valid_period(obs, time_range=('1996-01-01', None))
    data.to_csv(work_dir / 'nan_percentage_1996_2016_valid_period.csv', index=False)
    
    # test
    data = pd.read_csv(work_dir / 'nan_percentage_1996_2016_valid_period.csv')
    
    # add gauges that have nan percentage > 60 to ignore_list.txt
    work_dir2 = Path(r'c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\rhine\data\2-interim')
    high_nan_gauge_ignore_list(data, threshold=60, wflow_id_exception=[705], work_dir=work_dir2)
    
    # plot
    plot_gauge_map_nan_percentage(data, work_dir, time_range=('2005', '2016'))
    
    # Calculate statistics of nan_percentage
    nan_percentage_hist = data['nan_percentage_within_valid_period'].plot.hist(bins=20)
    plt.xlabel('NAN Percentage, 100%')
    plt.ylabel('Frequency')
    plt.title('Histogram of NaN Percentage')
    plt.show()
    
    # # filter out high nan percentage data
    # data_filtered = data[data['nan_percentage'] <= 60]
    # plot_gauge_map_nan_percentage(data_filtered, work_dir, time_range=('1996', '2016'))
    
    
    # modify obs dataset to remove filtered gauges
    # data = pd.read_csv(work_dir / 'nan_percentage_1996_2016.csv')
    # obs_filtered = remove_gauge_from_obs(obs, data, threshold=60)
    # obs_filtered.to_netcdf(work_dir / 'discharge_obs_hr_FORMAT_allvars_filtered60.nc')
    
    # check the nan value in Main
    # gauges selected
    temp = pd.read_csv(r'c:\Users\deng_jg\work\temp.csv')
    numbers = temp.values.flatten()
    numbers = [int(num) for num in numbers if pd.notna(num)]
    
    # Check for missing data in obs xarray dataset for wflow_id = 172
    missing_year = np.unique(obs.sel(wflow_id=172).Q.isnull().time.dt.year.values)
    
    for gauge in numbers:
        for year in [1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 
                 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 
                 2014, 2015]:
            percent = obs.sel(wflow_id=gauge, time=slice(str(year),str(year+1))).Q.isnull().values.mean()*100
            print(f'gauge {gauge} year {year}: {percent}% missing')
        
    
   
     
